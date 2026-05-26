//! Per-record scan orchestrator: raw scan + consensus-dimer scan +
//! method selection (raw_fallback / consensus_single / consensus_multi).

use super::consensus::{build_consensus_dimer, extract_consensus_monomers, ConsensusEntry};
use super::find_ssrs::{find_ssrs, Hit};
use super::{Config, MotifSpec};
use ahash::AHashMap;

/// Per-record summary in the prototype's column order.
#[derive(Debug, Clone)]
pub struct SummaryRow {
    pub record_id: String,
    pub length_bp: usize,
    pub ssr_flag: String, // "yes" / "no"
    /// "NA" or motif (upper canonical).
    pub dominant_motif: String,
    /// "NA" or motif length.
    pub dominant_motif_length: String,
    pub dominant_motif_repeats: u64,
    /// Rounded to 2 decimals.
    pub dominant_motif_coverage_pct: f64,
    pub total_ssr_coverage_pct: f64,
    /// "NA" when no motif; else `m:XX.X%;m2:YY.Y%` (one decimal).
    pub top_motifs: String,
    pub ssr_method: String, // raw_fallback / consensus_single / consensus_multi
    /// "NA" or integer.
    pub consensus_period_bp: String,
    /// "NA" or `KMER(count);KMER(count);…` uppercase.
    pub consensus_monomer: String,
    pub raw_dominant_motif: String,
    pub raw_dominant_motif_coverage_pct: f64,
    pub raw_total_coverage_pct: f64,
    pub raw_n_regions: u64,
    pub raw_top_motifs: String,
}

#[derive(Debug, Clone)]
pub struct RegionRow {
    pub record_id: String,
    pub ssr_number: u32,
    pub motif_length: usize,
    pub motif_sequence: String,
    pub repeats: usize,
    pub start: usize,
    pub end: usize,
    pub normalized_motif: String,
}

/// One per-canonical aggregation row in `by_motif`.
#[derive(Debug, Clone)]
struct ByMotif {
    motif_length: usize,
    total_repeats: u64,
    total_coverage_bp: u64,
    n_regions: u64,
}

/// Output of `scan_sequence` — hits + summary + per-canonical aggregation.
#[derive(Debug, Clone)]
struct ScanOutput {
    hits: Vec<Hit>,
    summary: PartialSummary,
    by_motif: Vec<(String, ByMotif)>,
}

#[derive(Debug, Clone)]
struct PartialSummary {
    ssr_flag: String,
    n_ssr_regions: u64,
    dominant_motif: String,
    dominant_motif_length: String,
    dominant_motif_repeats: u64,
    dominant_motif_coverage_pct: f64,
    total_ssr_coverage_pct: f64,
    top_motifs: String,
}

fn empty_summary(n_hits: u64) -> PartialSummary {
    PartialSummary {
        ssr_flag: "no".into(),
        n_ssr_regions: n_hits,
        dominant_motif: "NA".into(),
        dominant_motif_length: "NA".into(),
        dominant_motif_repeats: 0,
        dominant_motif_coverage_pct: 0.0,
        total_ssr_coverage_pct: 0.0,
        top_motifs: "NA".into(),
    }
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// Aggregate `find_ssrs` hits over `length_bp`. Returns the summary +
/// the per-canonical map (the multi-motif branch consumes the map).
fn summarise_hits(hits: &[Hit], length_bp: usize, threshold_pct: f64) -> ScanOutput {
    if hits.is_empty() || length_bp == 0 {
        return ScanOutput {
            hits: hits.to_vec(),
            summary: empty_summary(hits.len() as u64),
            by_motif: Vec::new(),
        };
    }
    let mut map: AHashMap<String, ByMotif> = AHashMap::new();
    // Preserve first-appearance order for stable Vec output below.
    let mut order: Vec<String> = Vec::new();
    for h in hits {
        let m = &h.normalized_motif;
        let entry = map.entry(m.clone()).or_insert_with(|| {
            order.push(m.clone());
            ByMotif {
                motif_length: h.motif_length,
                total_repeats: 0,
                total_coverage_bp: 0,
                n_regions: 0,
            }
        });
        let ssr_len = (h.end - h.start + 1) as u64;
        entry.total_repeats += h.repeats as u64;
        entry.total_coverage_bp += ssr_len;
        entry.n_regions += 1;
    }
    // Drop motifs that are integer multiples of a shorter one.
    let keep_set: std::collections::HashSet<String> = {
        let names: Vec<String> = order.clone();
        super::find_ssrs::get_unique_motifs(&names)
            .into_iter()
            .collect()
    };
    let by_motif: Vec<(String, ByMotif)> = order
        .into_iter()
        .filter(|n| keep_set.contains(n))
        .map(|n| {
            let v = map.remove(&n).unwrap();
            (n, v)
        })
        .collect();
    if by_motif.is_empty() {
        return ScanOutput {
            hits: hits.to_vec(),
            summary: empty_summary(hits.len() as u64),
            by_motif: Vec::new(),
        };
    }
    // dominant = max total_coverage_bp; tie-break: first inserted.
    let (dom_idx, _) = by_motif
        .iter()
        .enumerate()
        .max_by_key(|(_, (_, d))| d.total_coverage_bp)
        .unwrap();
    let (dom_m, dom_d) = by_motif[dom_idx].clone();
    let total_cov_bp: u64 = by_motif.iter().map(|(_, d)| d.total_coverage_bp).sum();
    let dom_pct = 100.0 * (dom_d.total_coverage_bp as f64) / (length_bp as f64);
    let total_pct = 100.0 * (total_cov_bp as f64) / (length_bp as f64);
    let mut sorted_top = by_motif.clone();
    sorted_top.sort_by_key(|(_, d)| std::cmp::Reverse(d.total_coverage_bp));
    let top_str = sorted_top
        .iter()
        .take(3)
        .map(|(m, d)| {
            format!(
                "{}:{:.1}%",
                m,
                100.0 * (d.total_coverage_bp as f64) / (length_bp as f64)
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    let summary = PartialSummary {
        ssr_flag: if dom_pct >= threshold_pct {
            "yes"
        } else {
            "no"
        }
        .into(),
        n_ssr_regions: hits.len() as u64,
        dominant_motif: dom_m,
        dominant_motif_length: dom_d.motif_length.to_string(),
        dominant_motif_repeats: dom_d.total_repeats,
        dominant_motif_coverage_pct: round2(dom_pct),
        total_ssr_coverage_pct: round2(total_pct),
        top_motifs: top_str,
    };
    ScanOutput {
        hits: hits.to_vec(),
        summary,
        by_motif,
    }
}

/// Multi-motif summary using RAW per-motif coverage. Used in
/// `consensus_multi` branch.
fn build_multimotif_summary(
    raw_by_motif: &[(String, ByMotif)],
    validated_canonicals: &[String],
    length_bp: usize,
    threshold_pct: f64,
) -> PartialSummary {
    let map: AHashMap<&str, &ByMotif> = raw_by_motif.iter().map(|(k, v)| (k.as_str(), v)).collect();
    let mut rows: Vec<(String, ByMotif)> = validated_canonicals
        .iter()
        .filter_map(|c| map.get(c.as_str()).map(|d| (c.clone(), (*d).clone())))
        .collect();
    if rows.is_empty() {
        return empty_summary(0);
    }
    rows.sort_by_key(|(_, d)| std::cmp::Reverse(d.total_coverage_bp));
    let (dom_m, dom_d) = rows[0].clone();
    let total_cov_bp: u64 = rows.iter().map(|(_, d)| d.total_coverage_bp).sum();
    let dom_pct = 100.0 * (dom_d.total_coverage_bp as f64) / (length_bp as f64);
    let total_pct_raw = 100.0 * (total_cov_bp as f64) / (length_bp as f64);
    let total_pct = total_pct_raw.min(100.0);
    let top_str = rows
        .iter()
        .take(3)
        .map(|(m, d)| {
            format!(
                "{}:{:.1}%",
                m,
                100.0 * (d.total_coverage_bp as f64) / (length_bp as f64)
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    let n_regions = rows.iter().map(|(_, d)| d.n_regions).sum();
    PartialSummary {
        ssr_flag: if dom_pct >= threshold_pct {
            "yes"
        } else {
            "no"
        }
        .into(),
        n_ssr_regions: n_regions,
        dominant_motif: dom_m,
        dominant_motif_length: dom_d.motif_length.to_string(),
        dominant_motif_repeats: dom_d.total_repeats,
        dominant_motif_coverage_pct: round2(dom_pct),
        total_ssr_coverage_pct: round2(total_pct),
        top_motifs: top_str,
    }
}

fn scan_sequence(seq: &[u8], specs: &[MotifSpec], threshold_pct: f64) -> ScanOutput {
    let hits = find_ssrs(seq, specs);
    summarise_hits(&hits, seq.len(), threshold_pct)
}

/// Top-level per-record scan. Returns the summary row + the raw-region
/// rows for the regions TSV.
pub fn scan_record(
    rec_id: &str,
    seq: &[u8],
    top_period: Option<usize>,
    cfg: &Config,
) -> (SummaryRow, Vec<RegionRow>) {
    let raw = scan_sequence(seq, &cfg.specs, cfg.ssr_flag_threshold_pct);
    let region_rows: Vec<RegionRow> = raw
        .hits
        .iter()
        .map(|h| RegionRow {
            record_id: rec_id.to_string(),
            ssr_number: h.ssr_number,
            motif_length: h.motif_length,
            motif_sequence: h.motif_sequence.clone(),
            repeats: h.repeats,
            start: h.start,
            end: h.end,
            normalized_motif: h.normalized_motif.clone(),
        })
        .collect();

    let mut method = "raw_fallback".to_string();
    let mut consensus_period_str = "NA".to_string();
    let mut consensus_monomer_str = "NA".to_string();
    let mut auth: PartialSummary = raw.summary.clone();

    if let Some(period) = top_period {
        consensus_period_str = period.to_string();
        let monomers: Vec<ConsensusEntry> = extract_consensus_monomers(
            seq,
            period,
            cfg.consensus_max_monomers,
            cfg.consensus_freq_ratio_min,
        );
        let mut validated: Vec<(ConsensusEntry, PartialSummary)> = Vec::new();
        for m in &monomers {
            let dimer = build_consensus_dimer(
                &m.kmer,
                cfg.consensus_dimer_copies,
                cfg.consensus_dimer_min_bp,
            );
            let dimer_scan =
                scan_sequence(dimer.as_bytes(), &cfg.specs, cfg.ssr_flag_threshold_pct);
            if dimer_scan.summary.ssr_flag == "yes" {
                validated.push((m.clone(), dimer_scan.summary));
            }
        }
        if !validated.is_empty() {
            consensus_monomer_str = validated
                .iter()
                .map(|(m, _)| format!("{}({})", m.kmer.to_ascii_uppercase(), m.count))
                .collect::<Vec<_>>()
                .join(";");
            let mut unique_found: Vec<String> = Vec::new();
            for (_, ds) in &validated {
                if !unique_found.iter().any(|u| u == &ds.dominant_motif) {
                    unique_found.push(ds.dominant_motif.clone());
                }
            }
            if unique_found.len() == 1 {
                auth = validated[0].1.clone();
                method = "consensus_single".into();
            } else {
                auth = build_multimotif_summary(
                    &raw.by_motif,
                    &unique_found,
                    seq.len(),
                    cfg.ssr_flag_threshold_pct,
                );
                method = "consensus_multi".into();
            }
        }
    }

    // v0.11: ssr_flag is recomputed from the array-scale raw total so
    // it tracks what fraction of the array is actually SSR. The
    // `auth.ssr_flag` from the consensus_single branch was derived
    // from the *dimer's* coverage (≈100% by construction) and would
    // fire `yes` for any array whose kite top period happened to
    // contain a known SSR motif — see the v0.10→v0.11 fix notes.
    let array_scale_ssr_flag = if raw.summary.total_ssr_coverage_pct >= cfg.ssr_flag_threshold_pct {
        "yes".to_string()
    } else {
        "no".to_string()
    };

    let summary = SummaryRow {
        record_id: rec_id.to_string(),
        length_bp: seq.len(),
        ssr_flag: array_scale_ssr_flag,
        dominant_motif: auth.dominant_motif,
        dominant_motif_length: auth.dominant_motif_length,
        dominant_motif_repeats: auth.dominant_motif_repeats,
        dominant_motif_coverage_pct: auth.dominant_motif_coverage_pct,
        total_ssr_coverage_pct: auth.total_ssr_coverage_pct,
        top_motifs: auth.top_motifs,
        ssr_method: method,
        consensus_period_bp: consensus_period_str,
        consensus_monomer: consensus_monomer_str,
        raw_dominant_motif: raw.summary.dominant_motif,
        raw_dominant_motif_coverage_pct: raw.summary.dominant_motif_coverage_pct,
        raw_total_coverage_pct: raw.summary.total_ssr_coverage_pct,
        raw_n_regions: raw.summary.n_ssr_regions,
        raw_top_motifs: raw.summary.top_motifs,
    };
    (summary, region_rows)
}

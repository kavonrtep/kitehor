//! `kitehor report` — defensive observation-only output.
//!
//! Reads a FASTA, runs kite + (peak clustering) + ssr-scan +
//! irregularity scan internally, and writes a single TSV with one
//! row per record. Cells contain raw measurements (peaks, scores,
//! coverages, event counts); no categorical verdicts, no cascade
//! decisions, no `combined_class`. Designed as a clean-slate
//! alternative to the `analyze`/`summary` pipeline when reliable
//! upstream calls aren't available.
//!
//! Output format: tab-separated columns; lists within a cell are
//! `;`-separated; fields within a list entry are `:`-separated.
//!
//! Stage outputs consumed:
//!   * kite::analyze → top peaks (all kept)
//!   * rule_classify::cluster_peaks → clustered peaks at the same tol
//!   * ssr::scan::scan_record → motif counts + coverage
//!   * irregularity::analyse_record → indel-event metrics
//!
//! Tandem-validate is NOT in this pipeline — its verdict column is a
//! label, not a measurement, and we deliberately exclude it.

pub mod io;

use anyhow::Result;
use rayon::prelude::*;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Config {
    pub kite: crate::kite::KiteConfig,
    pub ssr: crate::ssr::Config,
    pub irregularity: crate::irregularity::Config,
    /// Cluster-peaks relative-period tolerance. Default 0.015 = same as
    /// rule_classify's default.
    pub cluster_tol: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            kite: crate::kite::KiteConfig::default(),
            ssr: crate::ssr::Config::default(),
            irregularity: crate::irregularity::Config::default(),
            cluster_tol: 0.015,
        }
    }
}

/// One row of the report TSV — flat measurement bag.
#[derive(Debug, Clone)]
pub struct ReportRow {
    pub record_id: String,
    pub array_length: usize,

    pub kite_n_peaks: usize,
    /// `period:score2_norm;…` sorted by score2_norm desc.
    pub kite_peaks: String,
    pub kite_n_clusters: usize,
    /// `rep_period:total_score:n_peaks;…` sorted by total_score desc.
    /// `n_peaks` is the number of raw peaks merged into the cluster
    /// (singleton clusters carry `:1`).
    pub kite_clusters: String,

    pub ssr_total_coverage_pct: f64,
    pub ssr_dominant_motif: String,
    pub ssr_dominant_motif_length: String,
    pub ssr_dominant_motif_repeats: u64,
    pub ssr_dominant_coverage_pct: f64,
    pub ssr_top_motifs: String,

    pub irreg_flag: String,
    pub irreg_n_kmer_groups: usize,
    pub irreg_indel_event_count: Option<usize>,
    pub irreg_indel_burden_pct: Option<f64>,
    pub irreg_indel_max_shift_bp: Option<f64>,
    pub irreg_indel_drift_bp_per_kb: Option<f64>,
    pub irreg_dropout_rate_per_pair: Option<f64>,
    pub irreg_baseline_jitter_bp: Option<f64>,
}

/// Build the report row for one record by running kite + cluster +
/// ssr + irregularity.
pub fn build_row(rec_id: &str, seq: &[u8], cfg: &Config) -> ReportRow {
    let record = crate::sequence::ArrayRecord::from_raw(rec_id.to_string(), seq);

    // 1. Kite.
    let kite_res = crate::kite::analyze(&record, &cfg.kite);

    // Top period for downstream stages: rank-1 by score2_norm (matches
    // analyze.rs's choice).
    let top_period: Option<usize> = kite_res
        .peaks
        .iter()
        .max_by(|a, b| {
            a.score2_norm
                .partial_cmp(&b.score2_norm)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|p| p.period);

    // 2. Cluster the kite peaks.
    let peak_rows: Vec<crate::rule_classify::decide::PeakRow> = kite_res
        .peaks
        .iter()
        .enumerate()
        .map(|(i, p)| crate::rule_classify::decide::PeakRow {
            rank: (i + 1) as u32,
            period: p.period,
            score2_norm: p.score2_norm,
        })
        .collect();
    let mut clusters = crate::rule_classify::cluster_peaks(&peak_rows, cfg.cluster_tol);
    clusters.sort_by(|a, b| {
        b.total_score
            .partial_cmp(&a.total_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Encode peaks (already sorted by score desc when kite returns them).
    let kite_peaks_str = encode_peaks(&kite_res.peaks);
    let kite_clusters_str = encode_clusters(&clusters);

    // 3. SSR scan.
    let (ssr_row, _regions) = crate::ssr::scan::scan_record(rec_id, seq, top_period, &cfg.ssr);

    // 4. Irregularity scan.
    let irr = crate::irregularity::analyse_record(
        rec_id,
        seq,
        top_period.map(|p| p as f64),
        &cfg.irregularity,
    );

    ReportRow {
        record_id: rec_id.to_string(),
        array_length: kite_res.length_bp,

        kite_n_peaks: kite_res.peaks.len(),
        kite_peaks: kite_peaks_str,
        kite_n_clusters: clusters.len(),
        kite_clusters: kite_clusters_str,

        ssr_total_coverage_pct: ssr_row.raw_total_coverage_pct,
        ssr_dominant_motif: ssr_row.raw_dominant_motif.clone(),
        ssr_dominant_motif_length: ssr_row.dominant_motif_length.clone(),
        ssr_dominant_motif_repeats: ssr_row.dominant_motif_repeats,
        ssr_dominant_coverage_pct: ssr_row.raw_dominant_motif_coverage_pct,
        ssr_top_motifs: ssr_row.raw_top_motifs.clone(),

        irreg_flag: irr.flag.to_string(),
        irreg_n_kmer_groups: irr.n_kmer_groups,
        irreg_indel_event_count: irr.indel_event_count,
        irreg_indel_burden_pct: irr.indel_burden_pct,
        irreg_indel_max_shift_bp: irr.indel_max_shift_bp,
        irreg_indel_drift_bp_per_kb: irr.indel_drift_bp_per_kb,
        irreg_dropout_rate_per_pair: irr.dropout_rate_per_pair,
        irreg_baseline_jitter_bp: irr.baseline_jitter_bp,
    }
}

/// Encode kite peaks as `period:score2_norm;…` sorted by score desc.
fn encode_peaks(peaks: &[crate::kite::KitePeak]) -> String {
    let mut sorted: Vec<&crate::kite::KitePeak> = peaks.iter().collect();
    sorted.sort_by(|a, b| {
        b.score2_norm
            .partial_cmp(&a.score2_norm)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted
        .iter()
        .map(|p| {
            format!(
                "{}:{}",
                p.period,
                trim_zeros(&format!("{:.6}", p.score2_norm)),
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// Encode clusters as `rep_period:total_score:n_peaks;…`.
fn encode_clusters(clusters: &[crate::rule_classify::Cluster]) -> String {
    clusters
        .iter()
        .map(|c| {
            let rep = (c.rep_period * 100.0).round() / 100.0;
            format!(
                "{}:{}:{}",
                trim_zeros(&format!("{:.2}", rep)),
                trim_zeros(&format!("{:.6}", c.total_score)),
                c.n_peaks,
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn trim_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Entry point — read FASTA, parallel scan, write `<prefix>.report.tsv`.
pub fn run_subcommand(fasta: &Path, out_prefix: &Path, cfg: &Config) -> Result<usize> {
    let records = crate::ssr::io::read_fasta_ordered(fasta)?;
    let rows: Vec<ReportRow> = records
        .par_iter()
        .map(|(rid, seq)| build_row(rid, seq, cfg))
        .collect();
    let path = io::report_path(out_prefix);
    io::write_report(&path, &rows)?;
    Ok(rows.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_zeros_basics() {
        assert_eq!(trim_zeros("1.000000"), "1");
        assert_eq!(trim_zeros("1.500000"), "1.5");
        assert_eq!(trim_zeros("1.5"), "1.5");
        assert_eq!(trim_zeros("123"), "123");
        assert_eq!(trim_zeros("0.000000"), "0");
    }

    #[test]
    fn encode_clusters_format() {
        let cs = vec![
            crate::rule_classify::Cluster {
                rep_period: 100.5,
                total_score: 0.5,
                n_peaks: 2,
                min_rank: 1,
                periods: vec![100, 101],
            },
            crate::rule_classify::Cluster {
                rep_period: 200.0,
                total_score: 0.1,
                n_peaks: 1,
                min_rank: 3,
                periods: vec![200],
            },
        ];
        assert_eq!(encode_clusters(&cs), "100.5:0.5:2;200:0.1:1");
    }

    #[test]
    fn empty_peaks_yields_empty_strings() {
        let row = build_row("empty", &[b'N'; 50], &Config::default());
        assert_eq!(row.kite_n_peaks, 0);
        assert_eq!(row.kite_peaks, "");
        assert_eq!(row.kite_n_clusters, 0);
        // irregularity has no period when kite has no peaks.
        assert_eq!(row.irreg_flag, "no_period");
    }

    #[test]
    fn synthetic_tandem_pure_has_dominant_peak() {
        // 60-bp monomer with no obvious internal sub-period, tiled
        // 100× = 6 kb. We don't pin the exact top peak (kite legitimately
        // picks up integer-divisor sub-periods too); we just check that
        // the row was assembled, peaks are non-empty, irregularity is
        // clean.
        let monomer = b"ACGTAGCTACTAGCATCGATCGGCATCGCATAGCTAGCATCGATCGTAGCATCGATCGTC";
        assert_eq!(monomer.len(), 60);
        let seq: Vec<u8> = monomer.iter().cycle().take(6000).copied().collect();
        let row = build_row("synthetic", &seq, &Config::default());
        assert!(row.kite_n_peaks >= 1);
        assert!(!row.kite_peaks.is_empty());
        assert!(row.kite_n_clusters >= 1);
        // 60 should appear somewhere in the peaks string (top-1 or
        // top-2; kite may rank a divisor higher by score2 weighting).
        assert!(
            row.kite_peaks.contains("60:") || row.kite_clusters.contains("60:"),
            "expected period 60 somewhere; got peaks={} clusters={}",
            row.kite_peaks,
            row.kite_clusters,
        );
        // Clean tandem → irregularity is `ok` and step_count is 0.
        assert_eq!(row.irreg_flag, "ok");
        assert_eq!(row.irreg_indel_event_count, Some(0));
    }
}

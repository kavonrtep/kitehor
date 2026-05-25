//! End-to-end orchestrator — `kitehor analyze <fasta> -o <prefix>`.
//!
//! Runs five stages in sequence, sharing intermediate state in memory
//! while still writing every per-stage TSV to disk for debugging and
//! downstream analysis (TSV-per-stage is a hard contract — see
//! `docs/new/rule_proto_impl_plan.md` §0).
//!
//! Stage map:
//!
//! - kite-periodicity → `<prefix>.kite.tsv` + `.kite.peaks.tsv`
//! - rule-classify → `<prefix>.verdicts.tsv`
//! - tandem-validate → `<prefix>.tandem_validate.tsv`
//! - ssr-scan → `<prefix>.ssr.tsv` + `.ssr.regions.tsv`
//! - summary-merge → `<prefix>.summary.tsv`

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Aggregated config for the orchestrator.
#[derive(Debug, Clone, Default)]
pub struct Config {
    pub kite: crate::kite::KiteConfig,
    pub rule: crate::rule_classify::Config,
    pub tandem_validate: crate::tandem_validate::Config,
    pub ssr: crate::ssr::Config,
    pub summary: crate::summary::Config,
}

#[derive(Debug, Clone)]
pub struct Report {
    pub n_records: usize,
    pub n_hor: usize,
    pub n_tr: usize,
    pub n_unresolved: usize,
}

/// Run the full pipeline.
pub fn run(fasta: &Path, out_prefix: &Path, cfg: &Config) -> Result<Report> {
    run_with(fasta, out_prefix, cfg, None)
}

/// Run the full pipeline, optionally writing a FASTA-like periodogram
/// bundle to `periodogram_path`.
pub fn run_with(
    fasta: &Path,
    out_prefix: &Path,
    cfg: &Config,
    periodogram_path: Option<&Path>,
) -> Result<Report> {
    use crate::sequence::ArrayRecord;

    // 1. FASTA → records.
    let records: Vec<(String, Vec<u8>)> = crate::ssr::io::read_fasta_ordered(fasta)?;
    log::info!(
        "analyze: loaded {} record(s) from {:?}",
        records.len(),
        fasta
    );
    let array_records: Vec<ArrayRecord> = records
        .iter()
        .map(|(id, seq)| ArrayRecord::from_raw(id.clone(), seq))
        .collect();

    // 2. Kite (parallel over records).
    let kite_results: Vec<crate::kite::KiteResult> = array_records
        .par_iter()
        .map(|r| crate::kite::analyze(r, &cfg.kite))
        .collect();
    write_kite_outputs(out_prefix, &kite_results)?;
    if let Some(path) = periodogram_path {
        let n = crate::periodogram::write_periodogram_bundle(path, &kite_results, &cfg.kite)
            .with_context(|| format!("writing periodogram bundle to {:?}", path))?;
        log::info!(
            "analyze: periodogram bundle wrote {} record(s) to {:?}",
            n,
            path
        );
    }

    // 3. Rule-classify.
    let verdicts: Vec<crate::rule_classify::Verdict> = kite_results
        .iter()
        .map(|kr| crate::rule_classify::classify(kr, &cfg.rule))
        .collect();
    let verdicts_path = crate::rule_classify::io::verdicts_path(out_prefix);
    crate::rule_classify::io::write_verdicts(&verdicts_path, &verdicts)?;

    // 4a–4b. Tandem-validate + SSR. Independent — run in parallel.
    // Tandem-validate consumes the kite peaks; SSR consumes the kite
    // rank-1 period as a hint.
    let kite_peaks_by_id: ahash::AHashMap<String, Vec<crate::tandem_validate::scan::PeakRow>> =
        kite_results
            .iter()
            .map(|kr| {
                let rows = kr
                    .peaks
                    .iter()
                    .enumerate()
                    .map(|(i, p)| crate::tandem_validate::scan::PeakRow {
                        rank: (i + 1) as u32,
                        period: p.period,
                        score2_norm: p.score2_norm,
                    })
                    .collect();
                (kr.array_id.clone(), rows)
            })
            .collect();
    let top_periods: ahash::AHashMap<String, usize> = kite_results
        .iter()
        .filter_map(|kr| {
            kr.peaks
                .iter()
                .max_by(|a, b| {
                    a.score2_norm
                        .partial_cmp(&b.score2_norm)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|p| (kr.array_id.clone(), p.period))
        })
        .collect();

    // Convert verdicts → tandem_validate input shape (all verdicts kept,
    // including simple_tr and unresolved).
    let tv_verdicts: Vec<crate::tandem_validate::scan::VerdictRow> = verdicts
        .iter()
        .map(|v| crate::tandem_validate::scan::VerdictRow {
            case_id: v.case_id.clone(),
            verdict: match v.kind {
                crate::rule_classify::VerdictKind::Hor => "hor".to_string(),
                crate::rule_classify::VerdictKind::SimpleTr => "simple_tr".to_string(),
                crate::rule_classify::VerdictKind::Unresolved => "unresolved".to_string(),
            },
            founder: v.founder,
            tile: v.tile,
            multiplicity: v.multiplicity.map(|m| m as i64),
        })
        .collect();

    let (tv_res, ssr_res) = rayon::join(
        || -> Result<()> {
            let rows = crate::tandem_validate::scan::scan_records(
                &records,
                &tv_verdicts,
                &kite_peaks_by_id,
                &cfg.tandem_validate,
            );
            crate::tandem_validate::io::write_rows(
                &crate::tandem_validate::io::out_path(out_prefix),
                &rows,
            )?;
            Ok(())
        },
        || -> Result<()> {
            let mut sum_rows: Vec<crate::ssr::scan::SummaryRow> = Vec::new();
            let mut reg_rows: Vec<crate::ssr::scan::RegionRow> = Vec::new();
            for (rec_id, seq) in &records {
                let (s, r) = crate::ssr::scan::scan_record(
                    rec_id,
                    seq,
                    top_periods.get(rec_id).copied(),
                    &cfg.ssr,
                );
                sum_rows.push(s);
                reg_rows.extend(r);
            }
            crate::ssr::io::write_summary(&crate::ssr::io::summary_path(out_prefix), &sum_rows)?;
            crate::ssr::io::write_regions(&crate::ssr::io::regions_path(out_prefix), &reg_rows)?;
            Ok(())
        },
    );
    tv_res?;
    ssr_res?;

    // 5. Summary-merge — runs against the freshly-written TSVs so the
    // merge logic is exercised through the same code path as the
    // standalone subcommand.
    let verdicts_p = crate::rule_classify::io::verdicts_path(out_prefix);
    let tv_p = crate::tandem_validate::io::out_path(out_prefix);
    let ssr_p = crate::ssr::io::summary_path(out_prefix);
    let _ = crate::summary::run_subcommand(&verdicts_p, &tv_p, &ssr_p, out_prefix, &cfg.summary)?;

    // Final tally.
    let n_records = records.len();
    let mut n_hor = 0;
    let mut n_tr = 0;
    let mut n_unresolved = 0;
    for v in &verdicts {
        match v.kind {
            crate::rule_classify::VerdictKind::Hor => n_hor += 1,
            crate::rule_classify::VerdictKind::SimpleTr => n_tr += 1,
            crate::rule_classify::VerdictKind::Unresolved => n_unresolved += 1,
        }
    }
    Ok(Report {
        n_records,
        n_hor,
        n_tr,
        n_unresolved,
    })
}

/// Mirror `kite-periodicity`'s default output schema so the analyze
/// bundle includes the same `<prefix>.kite.tsv` + `.kite.peaks.tsv`
/// pair the standalone subcommand emits.
fn write_kite_outputs(out_prefix: &Path, results: &[crate::kite::KiteResult]) -> Result<()> {
    let primary = primary_path(out_prefix);
    let long = long_peaks_path(out_prefix);
    if let Some(parent) = primary.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f =
        std::fs::File::create(&primary).with_context(|| format!("creating {:?}", &primary))?;
    writeln!(
        f,
        "case_id\tarray_length\tn_peaks_kept\tmonomer_size\tscore\tmonomer_size_2\tscore_2\tmonomer_size_3\tscore_3"
    )?;
    let na = "NA";
    for r in results {
        let p1 = r.peaks.first();
        let p2 = r.peaks.get(1);
        let p3 = r.peaks.get(2);
        let fmt_p = |p: Option<&crate::kite::KitePeak>| {
            p.map(|x| x.period.to_string()).unwrap_or_else(|| na.into())
        };
        let fmt_s = |p: Option<&crate::kite::KitePeak>| {
            p.map(|x| format!("{:.10}", x.score))
                .unwrap_or_else(|| na.into())
        };
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.array_id,
            r.length_bp,
            r.peaks.len(),
            fmt_p(p1),
            fmt_s(p1),
            fmt_p(p2),
            fmt_s(p2),
            fmt_p(p3),
            fmt_s(p3),
        )?;
    }

    let mut g = std::fs::File::create(&long).with_context(|| format!("creating {:?}", &long))?;
    writeln!(
        g,
        "case_id\tarray_length\trank\tperiod\tpeak_height\tscore\tscore2\tscore2_norm\tbackground"
    )?;
    for r in results {
        for (i, p) in r.peaks.iter().enumerate() {
            writeln!(
                g,
                "{}\t{}\t{}\t{}\t{:.4}\t{:.10}\t{:.10}\t{:.10}\t{:.4}",
                r.array_id,
                r.length_bp,
                i + 1,
                p.period,
                p.peak_height,
                p.score,
                p.score2,
                p.score2_norm,
                p.background
            )?;
        }
    }
    Ok(())
}

fn primary_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".kite.tsv");
    PathBuf::from(p)
}

fn long_peaks_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".kite.peaks.tsv");
    PathBuf::from(p)
}

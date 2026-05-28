//! TSV writer for the `kitehor report` output.
//!
//! 19 tab-separated columns, one row per record. Cells that hold
//! lists (kite_peaks, kite_clusters, ssr_top_motifs) use `;`
//! between entries and `:` between fields within an entry. Empty
//! cells (never the string "NA") indicate missing/non-applicable
//! values — e.g. irregularity columns when the scan flag is
//! `no_period` or `too_short`.

use super::ReportRow;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

const COLUMNS: &[&str] = &[
    "record_id",
    "array_length",
    "kite_n_peaks",
    "kite_peaks",
    "kite_n_clusters",
    "kite_clusters",
    "ssr_total_coverage_pct",
    "ssr_dominant_motif",
    "ssr_dominant_motif_length",
    "ssr_dominant_motif_repeats",
    "ssr_dominant_coverage_pct",
    "ssr_top_motifs",
    "irreg_flag",
    "irreg_n_kmer_groups",
    "irreg_indel_event_count",
    "irreg_indel_burden_pct",
    "irreg_indel_max_shift_bp",
    "irreg_indel_drift_bp_per_kb",
    "irreg_dropout_rate_per_pair",
    "irreg_baseline_jitter_bp",
];

/// Output path `<prefix>.report.tsv`.
pub fn report_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".report.tsv");
    PathBuf::from(p)
}

pub fn write_report(path: &Path, rows: &[ReportRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", COLUMNS.join("\t"))?;
    for r in rows {
        let line = format!(
            "{rid}\t{l}\t{kn}\t{kp}\t{kc}\t{kcs}\t{stc}\t{sdm}\t{sml}\t{smr}\t{sdc}\t{stm}\t{ifl}\t{ing}\t{iec}\t{ibp}\t{ims}\t{idr}\t{idp}\t{ibj}",
            rid = r.record_id,
            l = r.array_length,
            kn = r.kite_n_peaks,
            kp = r.kite_peaks,
            kc = r.kite_n_clusters,
            kcs = r.kite_clusters,
            stc = fmt_pct(r.ssr_total_coverage_pct),
            sdm = ssr_str(&r.ssr_dominant_motif),
            sml = ssr_str(&r.ssr_dominant_motif_length),
            smr = r.ssr_dominant_motif_repeats,
            sdc = fmt_pct(r.ssr_dominant_coverage_pct),
            stm = ssr_str(&r.ssr_top_motifs),
            ifl = r.irreg_flag,
            ing = r.irreg_n_kmer_groups,
            iec = fmt_opt_usize(r.irreg_indel_event_count),
            ibp = fmt_opt_f64(r.irreg_indel_burden_pct),
            ims = fmt_opt_f64(r.irreg_indel_max_shift_bp),
            idr = fmt_opt_f64(r.irreg_indel_drift_bp_per_kb),
            idp = fmt_opt_f64(r.irreg_dropout_rate_per_pair),
            ibj = fmt_opt_f64(r.irreg_baseline_jitter_bp),
        );
        writeln!(f, "{}", line)?;
    }
    Ok(())
}

/// Treat the SSR module's "NA" sentinel as an empty cell (the report
/// uses empty rather than NA to indicate "no value").
fn ssr_str(s: &str) -> &str {
    if s == "NA" {
        ""
    } else {
        s
    }
}

fn fmt_pct(v: f64) -> String {
    if !v.is_finite() {
        return String::new();
    }
    trim_zeros(&format!("{:.4}", v))
}

fn fmt_opt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => trim_zeros(&format!("{:.6}", x)),
        Some(_) => String::new(),
        None => String::new(),
    }
}

fn fmt_opt_usize(v: Option<usize>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => String::new(),
    }
}

fn trim_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ssr_str_maps_na_to_empty() {
        assert_eq!(ssr_str("NA"), "");
        assert_eq!(ssr_str("AT"), "AT");
    }

    #[test]
    fn fmt_opt_f64_handles_none() {
        assert_eq!(fmt_opt_f64(None), "");
        assert_eq!(fmt_opt_f64(Some(0.5)), "0.5");
    }
}

//! TSV I/O for irregularity scan — 14-column schema matching the
//! Python prototype.

use super::scan::RecordResult;
use anyhow::{anyhow, Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

const COLUMNS: &[&str] = &[
    "record_id",
    "array_length",
    "period_P",
    "n_kmer_groups",
    "n_pairs_total",
    "baseline_jitter_bp",
    "indel_event_count",
    "indel_burden_pct",
    "indel_max_shift_bp",
    "indel_drift_bp_per_kb",
    "dropout_event_count",
    "dropout_rate_per_pair",
    "flag",
    "notes",
];

/// Output path for the irregularity TSV (`<prefix>.irregularity.tsv`).
pub fn irregularity_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".irregularity.tsv");
    PathBuf::from(p)
}

/// Read kite top periods (`monomer_size` column of kite summary TSV)
/// keyed by `case_id`. Returns an empty map on missing column. Used
/// by the standalone subcommand to pair FASTA records with their
/// kite-detected period.
pub fn read_kite_top_periods(path: &Path) -> Result<ahash::AHashMap<String, f64>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening kite summary at {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let cid_idx = headers
        .iter()
        .position(|h| h == "case_id" || h == "record_id")
        .ok_or_else(|| anyhow!("kite summary missing case_id/record_id column"))?;
    let ms_idx = headers
        .iter()
        .position(|h| h == "monomer_size")
        .ok_or_else(|| anyhow!("kite summary missing monomer_size column"))?;
    let mut out: ahash::AHashMap<String, f64> = ahash::AHashMap::new();
    for rec in rdr.records() {
        let rec = rec.context("reading kite summary row")?;
        let cid = rec.get(cid_idx).unwrap_or("").to_string();
        if cid.is_empty() {
            continue;
        }
        let s = rec.get(ms_idx).unwrap_or("").trim();
        if s.is_empty() || s.eq_ignore_ascii_case("na") || s.eq_ignore_ascii_case("nan") {
            continue;
        }
        if let Ok(p) = s.parse::<f64>() {
            if p > 0.0 {
                out.insert(cid, p);
            }
        }
    }
    Ok(out)
}

/// Write the irregularity TSV. 14 columns, `%.6g` float format, empty
/// cells for `None`. Matches `tools/rule_proto/irregularity_v2.py`.
pub fn write_irregularity(path: &Path, rows: &[RecordResult]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", COLUMNS.join("\t"))?;
    for r in rows {
        let line = format!(
            "{rid}\t{l}\t{p}\t{ng}\t{np}\t{bj}\t{ec}\t{bp}\t{ms}\t{dr}\t{dc}\t{drp}\t{flag}\t{notes}",
            rid = r.record_id,
            l = r.array_length,
            p = fmt_opt_f64(r.period_p),
            ng = r.n_kmer_groups,
            np = r.n_pairs_total,
            bj = fmt_opt_f64(r.baseline_jitter_bp),
            ec = fmt_opt_usize(r.indel_event_count),
            bp = fmt_opt_f64(r.indel_burden_pct),
            ms = fmt_opt_f64(r.indel_max_shift_bp),
            dr = fmt_opt_f64(r.indel_drift_bp_per_kb),
            dc = fmt_opt_usize(r.dropout_event_count),
            drp = fmt_opt_f64(r.dropout_rate_per_pair),
            flag = r.flag,
            notes = r.notes,
        );
        writeln!(f, "{}", line)?;
    }
    Ok(())
}

fn fmt_opt_f64(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{:.6}", x)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string(),
        Some(_) => "nan".into(),
        None => String::new(),
    }
}

fn fmt_opt_usize(v: Option<usize>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmt_opt_f64_drops_trailing_zeros() {
        assert_eq!(fmt_opt_f64(Some(1.0)), "1");
        assert_eq!(fmt_opt_f64(Some(1.5)), "1.5");
        assert_eq!(fmt_opt_f64(Some(0.0)), "0");
        assert_eq!(fmt_opt_f64(None), "");
    }
}

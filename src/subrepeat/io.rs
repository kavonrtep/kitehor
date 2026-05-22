//! TSV I/O for `subrepeat-scan`.

use super::scan::{PeakRow, SummaryRow, WindowRow};
use ahash::AHashMap;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn summary_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".subrepeat.tsv");
    PathBuf::from(p)
}

pub fn windows_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".windows.tsv");
    PathBuf::from(p)
}

/// Read a kite peaks TSV and return `{record_id: Vec<PeakRow>}`. Preserves
/// per-record peak order from the file.
pub fn read_kite_peaks_grouped(path: &Path) -> Result<AHashMap<String, Vec<PeakRow>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(hdr_line) = lines.next() else {
        return Ok(AHashMap::new());
    };
    let cols: Vec<&str> = hdr_line.split('\t').collect();
    let cid = cols
        .iter()
        .position(|c| *c == "case_id")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing case_id"))?;
    let rank = cols
        .iter()
        .position(|c| *c == "rank")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing rank"))?;
    let period = cols
        .iter()
        .position(|c| *c == "period")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing period"))?;
    let s2n = cols
        .iter()
        .position(|c| *c == "score2_norm")
        .ok_or_else(|| anyhow::anyhow!("kite peaks missing score2_norm"))?;
    let mut out: AHashMap<String, Vec<PeakRow>> = AHashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let id = cells[cid].to_string();
        let row = PeakRow {
            rank: cells[rank].parse().unwrap_or(0),
            period: cells[period].parse().unwrap_or(0),
            score2_norm: cells[s2n].parse().unwrap_or(0.0),
        };
        out.entry(id).or_default().push(row);
    }
    Ok(out)
}

/// Format a float the way pandas's `to_csv` does without `float_format`.
fn fmt_pandas_default(x: f64) -> String {
    if x.is_nan() {
        return "".to_string();
    }
    if x.is_infinite() {
        return if x.is_sign_negative() { "-inf" } else { "inf" }.to_string();
    }
    if x == x.trunc() && x.abs() < 1e16 {
        return format!("{}.0", x as i64);
    }
    let s = format!("{}", x);
    if !s.contains('.') && !s.contains('e') && !s.contains('E') {
        format!("{s}.0")
    } else {
        s
    }
}

pub fn write_summary(path: &Path, rows: &[SummaryRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(
        f,
        "record_id\tlength_bp\thost_period_bp\tsubrepeat_period_bp\tsubrepeat_flag\treason\tn_windows_total\tn_windows_sub\tn_windows_non_sub\tn_subrepeat_blocks\tsubrepeat_coverage_bp\tsubrepeat_coverage_pct\tblocks"
    )?;
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.record_id,
            r.length_bp,
            r.host_period_bp,
            r.subrepeat_period_bp,
            r.subrepeat_flag,
            r.reason,
            r.n_windows_total,
            r.n_windows_sub,
            r.n_windows_non_sub,
            r.n_subrepeat_blocks,
            r.subrepeat_coverage_bp,
            fmt_pandas_default(r.subrepeat_coverage_pct),
            r.blocks,
        )?;
    }
    Ok(())
}

pub fn write_windows(path: &Path, rows: &[WindowRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    if rows.is_empty() {
        writeln!(f)?;
        return Ok(());
    }
    writeln!(
        f,
        "record_id\twindow_start\twindow_end\ttop_period\ttop_score2_norm\tclass_raw\tclass"
    )?;
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.record_id,
            r.window_start,
            r.window_end,
            r.top_period,
            fmt_pandas_default(r.top_score2_norm),
            r.class_raw,
            r.class_,
        )?;
    }
    Ok(())
}

//! TSV I/O for `tandem-validate`. Output mirrors the Python prototype's
//! 16-column layout at `%.6g` precision so per-record `decision_hint`
//! diffs against the Python reference reduce to a string compare.

use super::scan::{CandidateKind, Decision, PeakRow, Row, VerdictRow};
use crate::rule_classify::io::fmt_g;
use ahash::AHashMap;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn out_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".tandem_validate.tsv");
    PathBuf::from(p)
}

/// Read a kite peaks TSV grouped by `case_id`, preserving per-record
/// peak order. Self-contained — does not depend on the `subrepeat`
/// module so this stage stands alone once the prior stages are retired.
pub fn read_kite_peaks_grouped(path: &Path) -> Result<AHashMap<String, Vec<PeakRow>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(hdr) = lines.next() else {
        return Ok(AHashMap::new());
    };
    let cols: Vec<&str> = hdr.split('\t').collect();
    let find = |name: &str| -> Result<usize> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| anyhow::anyhow!("kite peaks missing {name}"))
    };
    let cid = find("case_id")?;
    let rank = find("rank")?;
    let period = find("period")?;
    let s2n = find("score2_norm")?;
    let mut out: AHashMap<String, Vec<PeakRow>> = AHashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() <= cid.max(rank).max(period).max(s2n) {
            continue;
        }
        out.entry(cells[cid].to_string())
            .or_default()
            .push(PeakRow {
                rank: cells[rank].parse().unwrap_or(0),
                period: cells[period].parse().unwrap_or(0),
                score2_norm: cells[s2n].parse().unwrap_or(0.0),
            });
    }
    Ok(out)
}

/// Read `verdicts.tsv` and return one [`VerdictRow`] per record. Unlike
/// `hor_validate::io::read_verdicts`, we keep **every** verdict — the
/// detector runs on `hor`, `simple_tr`, and `unresolved` rows.
pub fn read_verdicts(path: &Path) -> Result<Vec<VerdictRow>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(hdr) = lines.next() else {
        return Ok(Vec::new());
    };
    let cols: Vec<&str> = hdr.split('\t').collect();
    let find = |name: &str| -> Result<usize> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing {name}"))
    };
    let cid = find("case_id")?;
    let vidx = find("verdict")?;
    let fidx = find("founder")?;
    let midx = find("multiplicity")?;
    let tidx = find("tile")?;
    let mut out: Vec<VerdictRow> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() <= cid.max(vidx).max(fidx).max(midx).max(tidx) {
            continue;
        }
        out.push(VerdictRow {
            case_id: cells[cid].to_string(),
            verdict: cells[vidx].to_string(),
            founder: cells[fidx].parse::<f64>().ok(),
            tile: cells[tidx].parse::<f64>().ok(),
            multiplicity: cells[midx].parse::<f64>().ok().map(|x| x.round() as i64),
        });
    }
    Ok(out)
}

pub fn write_rows(path: &Path, rows: &[Row]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    writeln!(
        f,
        "record_id\tverdict\thost_period\tmultiplicity\twindow_bp\tn_candidates\tcandidates\tbest_candidate_period\tbest_candidate_kind\tdensity\tspatial_contrast\tphase_contrast\tn_windows_total\tn_windows_present\tdecision_hint\treason"
    )?;
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.record_id,
            r.verdict,
            fmt_opt_f64(r.host_period, 6),
            fmt_opt_i64(r.multiplicity),
            fmt_opt_usize(r.window_bp),
            r.n_candidates,
            r.candidates,
            fmt_opt_f64(r.best_candidate_period, 6),
            fmt_opt_kind(r.best_candidate_kind),
            fmt_opt_f64(r.density, 6),
            fmt_opt_f64(r.spatial_contrast, 6),
            fmt_phase_contrast(r.phase_contrast),
            r.n_windows_total,
            r.n_windows_present,
            decision_label(r.decision_hint),
            r.reason,
        )?;
    }
    Ok(())
}

fn fmt_opt_f64(v: Option<f64>, precision: usize) -> String {
    match v {
        Some(x) => fmt_g(precision, x),
        None => String::new(),
    }
}

fn fmt_opt_i64(v: Option<i64>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => String::new(),
    }
}

fn fmt_opt_usize(v: Option<usize>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => String::new(),
    }
}

fn fmt_opt_kind(v: Option<CandidateKind>) -> String {
    match v {
        Some(CandidateKind::Founder) => "founder".to_string(),
        Some(CandidateKind::Other) => "other".to_string(),
        None => String::new(),
    }
}

/// Phase contrast renders as `"NA"` when not computed (window covers a
/// full host cycle so within-cycle bins are degenerate), matching the
/// Python prototype's literal string.
fn fmt_phase_contrast(v: Option<f64>) -> String {
    match v {
        Some(x) => fmt_g(6, x),
        None => "NA".to_string(),
    }
}

/// Stringify a [`Decision`] for both the `decision_hint` and `reason`
/// columns. See `Decision` for the variant → label mapping.
pub fn decision_label(d: Decision) -> &'static str {
    match d {
        Decision::LocalizedSubrepeat => "localized_subrepeat",
        Decision::ConfirmsHost => "confirms_host",
        Decision::Ambiguous => "ambiguous",
        Decision::NoSignal => "no_signal",
        Decision::NoCandidates => "no_candidates",
        Decision::NoWindows => "no_windows",
        Decision::SkipK2 => "skip_k2",
        Decision::NoHost => "no_host",
        Decision::NoVerdict => "no_verdict",
    }
}

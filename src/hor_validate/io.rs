//! TSV I/O for `hor-validate`. Format uses `%.6g` per the prototype.

use super::scan::{HorVerdict, ValidationRow};
use crate::rule_classify::io::fmt_g;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn out_path(out_prefix: &Path) -> PathBuf {
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".hor_within_tile.tsv");
    PathBuf::from(p)
}

/// Read `verdicts.tsv` and return one [`HorVerdict`] per HOR row.
pub fn read_verdicts(path: &Path) -> Result<Vec<HorVerdict>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let mut lines = text.lines();
    let Some(hdr) = lines.next() else {
        return Ok(Vec::new());
    };
    let cols: Vec<&str> = hdr.split('\t').collect();
    let cid = cols
        .iter()
        .position(|c| *c == "case_id")
        .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing case_id"))?;
    let vidx = cols
        .iter()
        .position(|c| *c == "verdict")
        .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing verdict"))?;
    let fidx = cols
        .iter()
        .position(|c| *c == "founder")
        .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing founder"))?;
    let midx = cols
        .iter()
        .position(|c| *c == "multiplicity")
        .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing multiplicity"))?;
    let tidx = cols
        .iter()
        .position(|c| *c == "tile")
        .ok_or_else(|| anyhow::anyhow!("verdicts.tsv missing tile"))?;
    let mut out: Vec<HorVerdict> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells[vidx] != "hor" {
            continue;
        }
        // Drop rows where `tile` doesn't parse — mirrors
        // `pd.to_numeric(..., errors="coerce")` + `.notna()`.
        let Ok(tile) = cells[tidx].parse::<f64>() else {
            continue;
        };
        let Ok(founder) = cells[fidx].parse::<f64>() else {
            continue;
        };
        let mult: Option<u32> = cells[midx].parse().ok();
        out.push(HorVerdict {
            case_id: cells[cid].to_string(),
            founder,
            tile,
            multiplicity: mult,
        });
    }
    Ok(out)
}

pub fn write_validation(path: &Path, rows: &[ValidationRow]) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    // Always emit the header — summary-merge reads this file by column
    // name and would fail on a header-less empty file. (The Python
    // prototype crashes when there are no HOR calls; the Rust port
    // turns that into "header only, no data rows".)
    writeln!(
        f,
        "record_id\tglobal_founder_bp\tglobal_tile_bp\tglobal_founder_score\tglobal_tile_score\tglobal_founder_tile_ratio\twithin_top_period\twithin_top_score\twithin_founder_score\twithin_founder_top_ratio\tdecision_hint\tfounder_density\tphase_contrast\tdensity_n_windows\tdensity_hint\tskip_reason"
    )?;
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.record_id,
            fmt_g(6, r.global_founder_bp),
            r.global_tile_bp,
            fmt_g(6, r.global_founder_score),
            fmt_g(6, r.global_tile_score),
            fmt_g(6, r.global_founder_tile_ratio),
            r.within_top_period,
            fmt_g(6, r.within_top_score),
            fmt_g(6, r.within_founder_score),
            fmt_g(6, r.within_founder_top_ratio),
            r.decision_hint,
            r.founder_density,
            r.phase_contrast,
            r.density_n_windows,
            r.density_hint,
            r.skip_reason,
        )?;
    }
    Ok(())
}

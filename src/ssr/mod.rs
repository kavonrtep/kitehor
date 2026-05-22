//! SSR (short-motif tandem repeat) scanner — Rust port of
//! `tools/rule_proto/ssr_scan.py`.
//!
//! TideCluster-style `find_ssrs` over motif lengths 1..=14 plus a
//! kite-driven consensus-dimer correction that handles interrupted
//! SSRs (e.g., `(GT)n` with mutations where raw coverage is ~9% but
//! the dimer is pure (GT)n).

pub mod consensus;
pub mod find_ssrs;
pub mod io;
pub mod scan;

use anyhow::Result;
use std::path::Path;

/// Per-motif-length minimum repeat count. Default per TideCluster:
/// `[(1, 20), (2, 9), (3, 6), (4..=14, 5)]`.
#[derive(Debug, Clone)]
pub struct MotifSpec {
    pub motif_length: usize,
    pub min_repeats: usize,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub ssr_flag_threshold_pct: f64,
    pub specs: Vec<MotifSpec>,
    pub consensus_dimer_copies: usize,
    pub consensus_dimer_min_bp: usize,
    pub consensus_max_monomers: usize,
    pub consensus_freq_ratio_min: f64,
}

impl Default for Config {
    fn default() -> Self {
        let mut specs = vec![
            MotifSpec {
                motif_length: 1,
                min_repeats: 20,
            },
            MotifSpec {
                motif_length: 2,
                min_repeats: 9,
            },
            MotifSpec {
                motif_length: 3,
                min_repeats: 6,
            },
        ];
        for n in 4..=14 {
            specs.push(MotifSpec {
                motif_length: n,
                min_repeats: 5,
            });
        }
        Self {
            ssr_flag_threshold_pct: 30.0,
            specs,
            consensus_dimer_copies: 4,
            consensus_dimer_min_bp: 30,
            consensus_max_monomers: 3,
            consensus_freq_ratio_min: 0.3,
        }
    }
}

/// Parse a `--motif-min-reps "1:20,2:9,3:6,…,14:5"` string into a
/// `Vec<MotifSpec>`. Used by the CLI handler.
pub fn parse_motif_min_reps(s: &str) -> anyhow::Result<Vec<MotifSpec>> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (l, r) = part
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("expected `len:min_reps`, got {:?}", part))?;
        let motif_length: usize = l
            .trim()
            .parse()
            .map_err(|e| anyhow::anyhow!("invalid motif length {:?}: {}", l, e))?;
        let min_repeats: usize = r
            .trim()
            .parse()
            .map_err(|e| anyhow::anyhow!("invalid min_repeats {:?}: {}", r, e))?;
        out.push(MotifSpec {
            motif_length,
            min_repeats,
        });
    }
    Ok(out)
}

/// Subcommand entry point.
pub fn run_subcommand(
    fasta: &Path,
    out_prefix: &Path,
    kite_peaks: Option<&Path>,
    cfg: &Config,
) -> Result<usize> {
    let records = io::read_fasta_ordered(fasta)?;
    let top_periods = match kite_peaks {
        Some(p) => io::read_kite_top_periods(p)?,
        None => Default::default(),
    };
    let mut summary_rows: Vec<scan::SummaryRow> = Vec::with_capacity(records.len());
    let mut region_rows: Vec<scan::RegionRow> = Vec::new();
    for (rec_id, seq) in &records {
        let (sum, regs) = scan::scan_record(rec_id, seq, top_periods.get(rec_id).copied(), cfg);
        summary_rows.push(sum);
        region_rows.extend(regs);
    }
    let summary_path = io::summary_path(out_prefix);
    let regions_path = io::regions_path(out_prefix);
    io::write_summary(&summary_path, &summary_rows)?;
    io::write_regions(&regions_path, &region_rows)?;
    Ok(summary_rows.len())
}

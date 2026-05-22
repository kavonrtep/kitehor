//! Subrepeat scanner — port of `tools/rule_proto/subrepeat_scan.py`.
//!
//! Detects arrays where a long host repeat unit alternates in space with
//! a short embedded sub-repeat (e.g., rDNA: ~10 kb host monomer
//! containing ~15× a 200 bp IGS sub-repeat). The global kite spectrum
//! shows BOTH periodicities, which can mislead the rule classifier.
//! This scanner slides per-record windows and asks: which of the two
//! candidate periods dominates *in this window?*

pub mod io;
pub mod scan;

use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub tol: f64,
    pub window_mult_sub: usize,
    pub step_frac: usize,
    pub top_n_sub: usize,
    pub top_n_host: usize,
    pub sub_floor: f64,
    pub window_score_floor: f64,
    pub min_run: usize,
    pub host_sub_ratio_min: usize,
    pub min_window_bp: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            tol: 0.05,
            window_mult_sub: 5,
            step_frac: 4,
            top_n_sub: 3,
            top_n_host: 10,
            sub_floor: 0.05,
            window_score_floor: 0.3,
            min_run: 3,
            host_sub_ratio_min: 3,
            min_window_bp: 1000,
        }
    }
}

/// Subcommand entry point.
pub fn run_subcommand(
    fasta: &Path,
    out_prefix: &Path,
    kite_peaks: &Path,
    cfg: &Config,
) -> Result<usize> {
    let records = crate::ssr::io::read_fasta_ordered(fasta)?;
    let global_peaks = io::read_kite_peaks_grouped(kite_peaks)?;
    let mut summary_rows: Vec<scan::SummaryRow> = Vec::with_capacity(records.len());
    let mut window_rows: Vec<scan::WindowRow> = Vec::new();
    for (rec_id, seq) in &records {
        let empty: Vec<scan::PeakRow> = Vec::new();
        let peaks = global_peaks.get(rec_id).unwrap_or(&empty);
        let (sum, wrows) = scan::scan_record(rec_id, seq, peaks, cfg);
        summary_rows.push(sum);
        window_rows.extend(wrows);
    }
    let sum_path = io::summary_path(out_prefix);
    let win_path = io::windows_path(out_prefix);
    io::write_summary(&sum_path, &summary_rows)?;
    io::write_windows(&win_path, &window_rows)?;
    Ok(summary_rows.len())
}

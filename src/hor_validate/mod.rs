//! Within-tile HOR validator — port of
//! `tools/rule_proto/hor_within_tile_check.py`.
//!
//! For each rule-classify HOR call, run kite on (a) the first tile and
//! (b) sliding windows across the whole array, then compute the
//! founder density + phase contrast. Output classifies each HOR call as
//! one of `{spatially_confirms_hor, localized_duplication, ambiguous,
//! insufficient_phase_bins, k_too_low_for_test(k=N), NA}`.

pub mod io;
pub mod scan;

use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub period_match_tol: f64,
    pub min_k_for_density: u32,
    pub density_window_tile_frac: usize,
    pub min_founder_mult: usize,
    pub min_density_window_bp: usize,
    pub max_density_windows: usize,
    pub density_rel_floor: f64,
    pub phase_fold_bins: usize,
    pub density_dup_max: f64,
    pub density_hor_min: f64,
    pub phase_contrast_dup_min: f64,
    pub phase_contrast_hor_max: f64,
    pub max_tile_bp: usize,
    pub min_window_bp: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            period_match_tol: 0.02,
            min_k_for_density: 4,
            density_window_tile_frac: 3,
            min_founder_mult: 3,
            min_density_window_bp: 200,
            max_density_windows: 1000,
            density_rel_floor: 0.2,
            phase_fold_bins: 10,
            density_dup_max: 0.35,
            density_hor_min: 0.7,
            phase_contrast_dup_min: 0.4,
            phase_contrast_hor_max: 0.15,
            max_tile_bp: 200_000,
            min_window_bp: 200,
        }
    }
}

/// Subcommand entry point.
pub fn run_subcommand(
    fasta: &Path,
    verdicts_path: &Path,
    global_peaks_path: &Path,
    out_prefix: &Path,
    cfg: &Config,
) -> Result<usize> {
    let records = crate::ssr::io::read_fasta_ordered(fasta)?;
    let verdicts = io::read_verdicts(verdicts_path)?;
    let global = crate::subrepeat::io::read_kite_peaks_grouped(global_peaks_path)?;
    let rows = scan::run(&records, &verdicts, &global, cfg);
    let out_path = io::out_path(out_prefix);
    io::write_validation(&out_path, &rows)?;
    Ok(rows.len())
}

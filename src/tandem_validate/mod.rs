//! Unified spatial-localization subrepeat detector — port of
//! `tools/rule_proto/tandem_validate.py` (spec v5).
//!
//! Replaces the older `hor_validate` (within-tile founder density+phase
//! test) and `subrepeat` (array-position nested-TR scan) stages with a
//! single detector. The two prior stages tested the same underlying
//! property — spatial heterogeneity of a sub-host periodicity — at two
//! different scales with two different in-window criteria. This module
//! unifies them under one scan loop.
//!
//! See `docs/new/tandem_validate_spec.md` for the algorithm contract
//! and `docs/new/tandem_validate_port_plan.md` for the port plan this
//! module implements.

pub mod io;
pub mod scan;

use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    // --- Candidate selection ---
    /// Reject candidates with `period < cand_min_period`. Default 20
    /// (rule_classify's floor; anything ≤ kite's `k=6` is noise).
    pub cand_min_period: usize,
    /// Absolute `score2_norm` floor. Default 0 — admit all kite peaks
    /// that pass kite's own peak-above-background filter.
    pub cand_score_floor: f64,
    /// Relative score floor: candidate's `score2_norm` must be at least
    /// this fraction of the record's max `score2_norm`. Adaptive
    /// per-record. Default 0.03.
    pub cand_rel_score_floor: f64,
    /// At most this many non-founder candidates per record. Default 5.
    pub cand_top_n: usize,
    /// Candidate period must be inside `host * host_inside_ratio`
    /// (strict for `Other`; +1% slack for `Founder`). Default 1/3.
    pub host_inside_ratio: f64,
    /// Tolerance for matching candidates to founder-harmonic explained
    /// rungs. Default 0.05.
    pub founder_tol: f64,

    // --- Window sizing ---
    /// Floor on window size = `host * window_host_frac`. Default 1/3.
    pub window_host_frac: f64,
    /// Window must hold `window_cand_mult × max_candidate` bp for kite
    /// to detect the candidate. Default 3.
    pub window_cand_mult: f64,
    /// Hard lower bound on window size (bp). Default 200.
    pub min_window_bp: usize,

    // --- In-window presence (two criteria, gated by Candidate::kind) ---
    /// Period-match tolerance (`|p_top − p_cand| / p_cand ≤ tol`).
    /// Default 0.02.
    pub period_match_tol: f64,
    /// STRICT (kind = Other): window's top score must be ≥ this for
    /// "present". Discriminates legit nested (heterogeneous top scores)
    /// from plain tandem (uniformly strong tops). Default 0.3.
    pub window_score_floor: f64,
    /// LOOSE (kind = Founder): in-window founder score is "present" if
    /// `sum(scores within ±tol of founder) ≥ presence_rel_floor × top`.
    /// The founder is rank-2/3 of kite in a clean HOR (the tile is
    /// rank-1), so strict would never fire. Default 0.2.
    pub presence_rel_floor: f64,

    // --- Spatial / phase binning ---
    /// Number of bins for spatial + phase contrast. Default 10.
    pub n_bins: usize,

    // --- Decision thresholds ---
    /// `density ≤ this` → localized. Default 0.35.
    pub density_dup_max: f64,
    /// `density ≥ this AND contrasts low` → uniform. Default 0.7.
    pub density_hor_min: f64,
    /// `contrast ≥ this` → localized. Default 0.4.
    pub contrast_dup_min: f64,
    /// `contrast ≤ this` (both spatial and phase) required for uniform.
    /// Default 0.15.
    pub contrast_hor_max: f64,
    /// `n_present < this` short-circuits to `no_signal`. Default 3.
    pub min_present_windows: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            cand_min_period: 20,
            cand_score_floor: 0.0,
            cand_rel_score_floor: 0.03,
            cand_top_n: 5,
            host_inside_ratio: 1.0 / 3.0,
            founder_tol: 0.05,
            window_host_frac: 1.0 / 3.0,
            window_cand_mult: 3.0,
            min_window_bp: 200,
            period_match_tol: 0.02,
            window_score_floor: 0.3,
            presence_rel_floor: 0.2,
            n_bins: 10,
            density_dup_max: 0.35,
            density_hor_min: 0.7,
            contrast_dup_min: 0.4,
            contrast_hor_max: 0.15,
            min_present_windows: 3,
        }
    }
}

/// Subcommand entry point: read inputs, scan, write `<prefix>.tandem_validate.tsv`.
pub fn run_subcommand(
    fasta: &Path,
    verdicts_path: &Path,
    peaks_path: &Path,
    out_prefix: &Path,
    cfg: &Config,
) -> Result<usize> {
    let records = crate::ssr::io::read_fasta_ordered(fasta)?;
    let verdicts = io::read_verdicts(verdicts_path)?;
    let peaks = io::read_kite_peaks_grouped(peaks_path)?;
    let rows = scan::scan_records(&records, &verdicts, &peaks, cfg);
    let out = io::out_path(out_prefix);
    io::write_rows(&out, &rows)?;
    Ok(rows.len())
}

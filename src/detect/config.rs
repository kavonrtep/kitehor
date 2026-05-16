//! Detector configuration (`detect_impl_plan.md §6.0`, A11).
//!
//! Every numeric threshold the detector uses lives here. Defaults
//! match the upstream spec §7–§9. Override via `--config detect.toml`.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::path::Path;

#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct DetectorConfig {
    // Widths
    pub min_width: usize,
    pub max_width: usize,
    pub max_widths_per_array: usize,
    pub neighborhood_n: usize,
    pub max_hor_k: usize,
    pub divisor_top_n: usize,

    // Wrap / IC
    pub min_rows_per_width: usize,
    pub ic_threshold_min: f64,
    pub ic_threshold_hor_base: f64,
    pub ic_threshold_hor_unit: f64,
    pub ic_threshold_simple_tr: f64,
    /// Simple-TR rescue threshold on R(1). When base IC is below
    /// `ic_threshold_simple_tr` but R(1) ≥ this value and phase_sep
    /// is low, fire simple_TR. Catches indel-drifted simple TRs
    /// where column IC is diluted to ~0.14 but adjacent-row
    /// similarity remains ~0.9.
    pub simple_tr_r1_rescue: f64,
    /// Lower IC floor used only by the rescue path and the
    /// "any-width-supported?" early-exit check. Distinct from
    /// `ic_threshold_min` so we can admit indel-drifted widths
    /// (IC ≈ 0.10–0.15) for rescue without polluting the general
    /// candidate path.
    pub ic_threshold_rescue: f64,

    // Embeddings
    pub embedding_k: usize,
    pub embedding_dim_hash: Option<usize>,

    // Phase / multiplicity
    pub phase_separation_threshold: f64,
    pub primitive_correction_delta: f64,
    pub min_hor_units: usize,
    /// R(1) floor at base width for regime-B HOR. Below this, the
    /// HOR is in regime C — slot consensuses too diverged to call
    /// the base period statistically valid; fall to simple_TR at
    /// HOR-unit width.
    pub regime_c_r1_threshold: f64,
    /// Input-period-score threshold for the "coexisting repeats →
    /// mixed" fallback. Two+ input periods with `period_score >=
    /// strong_period_score` after multiplicity dedup → mixed.
    pub strong_period_score: f64,
    /// Regime-A tag fires when the simple_TR candidate's
    /// `|R(best_lag) - R(1)| < regime_a_r_curve_flatness` AND
    /// `R(1) > regime_a_r1_floor`. Reason field then mentions
    /// "regime A".
    pub regime_a_r_curve_flatness: f64,
    pub regime_a_r1_floor: f64,

    // Shift
    pub shift_local_range_bp: i32,
    pub shift_breakpoint_threshold: i32,
    pub shift_breakpoint_window_frac: f64,

    // Irregularity / segmentation
    pub block_size_rows_min: usize,
    pub stratification_same_threshold: f64,
    pub stratification_diff_threshold: f64,
    /// Demote `HOR` → `irregular_HOR` when block-level IC variance
    /// (irregularity_score) exceeds this. Default 0.30 (DH3).
    pub irregularity_demote_threshold: f64,

    // M7 — analysis blocks for same-width mixed detection
    // (`docs/new/detect_m7_plan.md` Q1 + Q8).
    /// Upper bound on the number of internal analysis blocks used
    /// per array for the mixed-detection consensus-identity test.
    /// `block_rows = max(min_segment_rows, ceil(n_rows / max_segments_per_array))`.
    pub max_segments_per_array: usize,
    /// Minimum rows in an analysis block. Smaller blocks are
    /// merged into their neighbour or dropped before the
    /// identity test.
    pub min_segment_rows: usize,
    /// Minimum non-N positions (as a fraction of the consensus
    /// length) required for a pairwise identity comparison to be
    /// admitted into the mixed-override test. Pairs below this
    /// floor return `None` (uninformative).
    pub min_identity_coverage: f64,
    /// For HOR / irregular_HOR arrays, the comparison consensus
    /// is built at `hor_length_bp` only when each analysis block
    /// contains at least this many complete HOR units. Otherwise
    /// fall back to `base_width_bp`.
    pub min_complete_units_per_block: usize,

    // Confidence
    pub confidence_weights: ConfidenceWeights,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ConfidenceWeights {
    pub alpha: f64,   // phase_separation
    pub beta:  f64,   // R(best_k) − R(unrelated_k)
    pub gamma: f64,   // log10(n_complete_copies + 1)
    pub delta: f64,   // mean_column_IC
    pub epsilon: f64, // irregularity_score
    pub zeta: f64,    // wobble_amplitude / w
    pub eta: f64,    // |mean_shift| / w
}

impl Default for ConfidenceWeights {
    fn default() -> Self {
        // Starting values from `detect_spec.md` §9 — calibrate in M6.
        Self {
            alpha: 3.0,
            beta: 2.0,
            gamma: 0.5,
            delta: 2.0,
            epsilon: 4.0,
            zeta: 1.0,
            eta: 2.0,
        }
    }
}

impl Default for DetectorConfig {
    fn default() -> Self {
        Self {
            min_width: 20,
            max_width: 5_000,
            max_widths_per_array: 40,
            neighborhood_n: 3,
            max_hor_k: 30,
            divisor_top_n: 5,

            min_rows_per_width: 8,
            // M6 calibration on ground_truth_v2. Lowered IC
            // thresholds catch (a) high-divergence HOR units whose
            // column IC drops because slot consensuses look random,
            // and (b) indel-affected simple TRs whose IC is diluted
            // by row drift but R(1) stays high.
            ic_threshold_min: 0.20,
            ic_threshold_hor_base: 0.30,
            ic_threshold_hor_unit: 0.25,
            ic_threshold_simple_tr: 0.30,
            simple_tr_r1_rescue: 0.40,
            ic_threshold_rescue: 0.01,

            embedding_k: 4,
            embedding_dim_hash: None,

            // M6 calibration: progressively lowered to 0.01 so
            // heavy-wobble HORs (where row misalignment weakens
            // R(k)) still cross the threshold while clean simple
            // TRs (phase_sep ≈ 0.0001) stay below.
            phase_separation_threshold: 0.01,
            primitive_correction_delta: 0.05,
            // M6: lowered from 5 → 3 to catch HOR fixtures with
            // n_copies=50 and k=16 (only 3 complete units).
            min_hor_units: 3,
            regime_c_r1_threshold: 0.5,
            strong_period_score: 0.85,
            regime_a_r_curve_flatness: 0.05,
            regime_a_r1_floor: 0.85,

            shift_local_range_bp: 5,
            shift_breakpoint_threshold: 3,
            shift_breakpoint_window_frac: 0.5,

            block_size_rows_min: 100,
            stratification_same_threshold: 0.90,
            stratification_diff_threshold: 0.80,
            irregularity_demote_threshold: 0.50,

            // M7 defaults (`docs/new/detect_m7_plan.md` Q8).
            max_segments_per_array: 32,
            min_segment_rows: 20,
            min_identity_coverage: 0.70,
            min_complete_units_per_block: 3,

            confidence_weights: ConfidenceWeights::default(),
        }
    }
}

impl DetectorConfig {
    /// Validate ranges that aren't enforced by the type system.
    pub fn validate(&self) -> Result<()> {
        if self.min_width == 0 || self.max_width <= self.min_width {
            anyhow::bail!(
                "min_width ({}) must be > 0 and < max_width ({})",
                self.min_width,
                self.max_width
            );
        }
        if self.max_widths_per_array == 0 {
            anyhow::bail!("max_widths_per_array must be >= 1");
        }
        if self.max_hor_k < 2 {
            anyhow::bail!("max_hor_k must be >= 2");
        }
        if self.embedding_k < 2 || self.embedding_k > 8 {
            anyhow::bail!(
                "embedding_k must be in [2, 8] (default 4); got {}",
                self.embedding_k
            );
        }
        if !(0.0..=1.0).contains(&self.phase_separation_threshold) {
            anyhow::bail!(
                "phase_separation_threshold must be in [0,1]; got {}",
                self.phase_separation_threshold
            );
        }
        if self.shift_local_range_bp < 1 {
            anyhow::bail!("shift_local_range_bp must be >= 1");
        }
        if self.stratification_diff_threshold > self.stratification_same_threshold {
            anyhow::bail!(
                "stratification_diff_threshold ({}) must be <= stratification_same_threshold ({})",
                self.stratification_diff_threshold,
                self.stratification_same_threshold
            );
        }
        if self.max_segments_per_array < 2 {
            anyhow::bail!(
                "max_segments_per_array ({}) must be >= 2",
                self.max_segments_per_array
            );
        }
        if self.min_segment_rows < 2 {
            anyhow::bail!(
                "min_segment_rows ({}) must be >= 2",
                self.min_segment_rows
            );
        }
        if !(0.0..=1.0).contains(&self.min_identity_coverage) {
            anyhow::bail!(
                "min_identity_coverage ({}) must be in [0, 1]",
                self.min_identity_coverage
            );
        }
        if self.min_complete_units_per_block < 1 {
            anyhow::bail!(
                "min_complete_units_per_block ({}) must be >= 1",
                self.min_complete_units_per_block
            );
        }
        Ok(())
    }

    /// Load from a TOML file. Missing fields fall back to `Default`.
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading detector config {:?}", path))?;
        let cfg: DetectorConfig = toml::from_str(&text)
            .with_context(|| format!("parsing detector config {:?}", path))?;
        cfg.validate()?;
        Ok(cfg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_passes_validation() {
        DetectorConfig::default().validate().unwrap();
    }

    #[test]
    fn rejects_bad_widths() {
        let mut c = DetectorConfig::default();
        c.min_width = 0;
        assert!(c.validate().is_err());
        c = DetectorConfig::default();
        c.max_width = c.min_width;
        assert!(c.validate().is_err());
    }

    #[test]
    fn rejects_bad_embedding_k() {
        let mut c = DetectorConfig::default();
        c.embedding_k = 1;
        assert!(c.validate().is_err());
        c.embedding_k = 9;
        assert!(c.validate().is_err());
    }

    #[test]
    fn rejects_inverted_stratification_thresholds() {
        let mut c = DetectorConfig::default();
        c.stratification_diff_threshold = 0.95;
        c.stratification_same_threshold = 0.90;
        assert!(c.validate().is_err());
    }

    #[test]
    fn loads_from_partial_toml() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("detect.toml");
        std::fs::write(&p, "min_width = 10\nmax_width = 1000\n").unwrap();
        let c = DetectorConfig::load(&p).unwrap();
        assert_eq!(c.min_width, 10);
        assert_eq!(c.max_width, 1000);
        // Defaults preserved for everything else.
        assert_eq!(c.max_hor_k, 30);
    }
}

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

    // Embeddings
    pub embedding_k: usize,
    pub embedding_dim_hash: Option<usize>,

    // Phase / multiplicity
    pub phase_separation_threshold: f64,
    pub primitive_correction_delta: f64,
    pub min_hor_units: usize,

    // Shift
    pub shift_local_range_bp: i32,
    pub shift_breakpoint_threshold: i32,
    pub shift_breakpoint_window_frac: f64,

    // Irregularity / segmentation
    pub block_size_rows_min: usize,
    pub stratification_same_threshold: f64,
    pub stratification_diff_threshold: f64,

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
            ic_threshold_min: 0.5,
            ic_threshold_hor_base: 0.4,
            ic_threshold_hor_unit: 0.7,
            ic_threshold_simple_tr: 0.7,

            embedding_k: 4,
            embedding_dim_hash: None,

            phase_separation_threshold: 0.15,
            primitive_correction_delta: 0.05,
            min_hor_units: 5,

            shift_local_range_bp: 5,
            shift_breakpoint_threshold: 3,
            shift_breakpoint_window_frac: 0.5,

            block_size_rows_min: 100,
            stratification_same_threshold: 0.90,
            stratification_diff_threshold: 0.80,

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

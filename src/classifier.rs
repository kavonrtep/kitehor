//! Kite-first probabilistic HOR classifier.
//!
//! Ports the calibrated R pipeline (predict_verdict.R + fit_platt.R)
//! into Rust. Three pieces:
//!   1. [`RandomForest`] — loads ranger trees from the JSON dump produced
//!      by `eval/training_data/export_ranger.R` and predicts (probability
//!      or regression mean) by averaging over trees.
//!   2. [`PlattScaler`] — sigmoid(intercept + slope · logit(s_raw)).
//!   3. [`ClassifierConfig`] — TOML-loadable bundle of thresholds, Platt
//!      coefficients, imputation medians, recovery tolerances.
//!
//! The 4-verdict logic + family demotion + k-recovery live in
//! [`crate::classify`] (the orchestrator); this module is pure model.
//!
//! Inequality direction (verified against ranger 0.16): for numeric
//! features, `x <= split_val` goes to the **left** child, otherwise right.
//! All features in our two models are numeric (no categorical splits).

use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Random forest (load + predict)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeType {
    /// Probability estimation (binary): leaf prediction is P(class = "1").
    Probability,
    /// Regression: leaf prediction is the mean target.
    Regression,
}

#[derive(Debug, Deserialize)]
struct RawTree {
    split_var: Vec<i32>,
    split_val: Vec<f64>,
    left: Vec<u32>,
    right: Vec<u32>,
    prediction: Vec<f64>,
}

#[derive(Debug, Deserialize)]
struct RawForest {
    treetype: String,
    num_trees: usize,
    feature_names: Vec<String>,
    #[serde(default)]
    class_levels: Option<Vec<String>>,
    trees: Vec<RawTree>,
}

/// Flat-array tree (one struct-of-arrays per tree). Internal nodes have
/// `split_var[i] >= 0`; leaves have `split_var[i] = -1` and the predicted
/// value sits in `prediction[i]`.
#[derive(Debug, Clone)]
pub struct Tree {
    pub split_var: Vec<i32>,
    pub split_val: Vec<f64>,
    pub left: Vec<u32>,
    pub right: Vec<u32>,
    pub prediction: Vec<f64>,
}

impl Tree {
    /// Predict a single leaf value. `x` is a feature row in the order of
    /// the forest's `feature_names`.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let mut node: usize = 0;
        loop {
            let sv = self.split_var[node];
            if sv < 0 {
                return self.prediction[node];
            }
            let v = x[sv as usize];
            // ranger: x <= split_val -> left, else right.
            if v <= self.split_val[node] {
                node = self.left[node] as usize;
            } else {
                node = self.right[node] as usize;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct RandomForest {
    pub treetype: TreeType,
    pub feature_names: Vec<String>,
    /// For probability mode, the levels (e.g., ["0", "1"]); leaf
    /// prediction is P(level = class_levels[1]).
    pub class_levels: Option<Vec<String>>,
    pub trees: Vec<Tree>,
}

impl RandomForest {
    pub fn load_json<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("opening RF model {:?}", path))?;
        let raw: RawForest = serde_json::from_str(&text)
            .with_context(|| format!("parsing RF model {:?}", path))?;
        let treetype = match raw.treetype.as_str() {
            "probability" => TreeType::Probability,
            "regression" => TreeType::Regression,
            other => bail!("unknown treetype: {}", other),
        };
        if raw.trees.len() != raw.num_trees {
            bail!(
                "tree count mismatch: declared {}, actual {}",
                raw.num_trees,
                raw.trees.len()
            );
        }
        let trees = raw
            .trees
            .into_iter()
            .map(|t| Tree {
                split_var: t.split_var,
                split_val: t.split_val,
                left: t.left,
                right: t.right,
                prediction: t.prediction,
            })
            .collect();
        Ok(RandomForest {
            treetype,
            feature_names: raw.feature_names,
            class_levels: raw.class_levels,
            trees,
        })
    }

    /// Average the per-tree leaf predictions. For probability mode this
    /// is the calibrated/raw P(class=1); for regression it is the mean.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let n = self.trees.len();
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for tree in &self.trees {
            sum += tree.predict(x);
        }
        sum / (n as f64)
    }

    /// Return a feature-name -> index map for the forest's expected
    /// feature ordering. Used to assemble feature vectors from a
    /// named row.
    pub fn feature_index(&self) -> HashMap<&str, usize> {
        self.feature_names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect()
    }

    /// Build a feature vector for this forest from a name->value map.
    /// Missing keys → NaN (caller must impute first; we do **not**
    /// silently zero-fill).
    pub fn assemble<'a, F: Fn(&str) -> Option<f64>>(&self, lookup: F) -> Vec<f64> {
        self.feature_names
            .iter()
            .map(|n| lookup(n.as_str()).unwrap_or(f64::NAN))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Platt scaling
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct PlattScaler {
    pub intercept: f64,
    pub slope: f64,
    pub clip_eps: f64,
}

impl PlattScaler {
    /// Map a raw RF probability in [0,1] to a calibrated probability via
    /// `sigmoid(intercept + slope · logit(s))`. Matches the R inference
    /// pipeline (see `predict_verdict.R`).
    pub fn calibrate(&self, s_raw: f64) -> f64 {
        let eps = self.clip_eps;
        let s = s_raw.clamp(eps, 1.0 - eps);
        let logit = (s / (1.0 - s)).ln();
        let z = self.intercept + self.slope * logit;
        // Numerically-stable sigmoid.
        if z >= 0.0 {
            let e = (-z).exp();
            1.0 / (1.0 + e)
        } else {
            let e = z.exp();
            e / (1.0 + e)
        }
    }
}

// ---------------------------------------------------------------------------
// Config (TOML)
// ---------------------------------------------------------------------------

/// Embedded default config, baked into the binary at compile time. Loaded
/// via `ClassifierConfig::default_baked()` when no `--classifier-config`
/// path is given.
pub const DEFAULT_CONFIG_TOML: &str =
    include_str!(concat!("../config/classifier.toml"));

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Thresholds {
    pub t_low: f64,
    pub t_high: f64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PlattCfg {
    pub intercept: f64,
    pub slope: f64,
    pub clip_eps: f64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Imputation {
    pub h_d1: f64,
    pub h_founder: f64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Recovery {
    pub top_k_candidates: usize,
    pub tol_bp: i64,
    pub tol_rel: f64,
    pub k_min: i64,
    pub founder_min_bp: i64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct TrainingMeta {
    pub detectable_thr: f64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ModelsCfg {
    pub hor_score: String,
    pub k_pred: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ClassifierConfig {
    pub thresholds: Thresholds,
    pub platt: PlattCfg,
    pub imputation: Imputation,
    pub recovery: Recovery,
    pub training_meta: TrainingMeta,
    pub models: ModelsCfg,
}

impl ClassifierConfig {
    pub fn default_baked() -> Result<Self> {
        toml::from_str(DEFAULT_CONFIG_TOML)
            .context("parsing baked-in default classifier.toml")
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("opening classifier config {:?}", path))?;
        toml::from_str(&text)
            .with_context(|| format!("parsing classifier config {:?}", path))
    }

    pub fn platt(&self) -> PlattScaler {
        PlattScaler {
            intercept: self.platt.intercept,
            slope: self.platt.slope,
            clip_eps: self.platt.clip_eps,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_config_parses() {
        let cfg = ClassifierConfig::default_baked().expect("default config");
        assert!(cfg.thresholds.t_low > 0.0 && cfg.thresholds.t_low < cfg.thresholds.t_high);
        assert!(cfg.thresholds.t_high < 1.0);
        assert!(cfg.platt.slope > 0.0);
    }

    #[test]
    fn platt_identity_when_slope1_intercept0() {
        let p = PlattScaler { intercept: 0.0, slope: 1.0, clip_eps: 1e-4 };
        // Sigmoid(logit(s)) = s for s away from the edges.
        for s in [0.05, 0.2, 0.5, 0.8, 0.95] {
            let out = p.calibrate(s);
            assert!((out - s).abs() < 1e-9, "s={s}, out={out}");
        }
    }

    #[test]
    fn platt_calibrates_using_baked_coefs() {
        // From classifier.toml: intercept=0.16844, slope=1.56331.
        // logit(0.5) = 0 -> z = intercept -> sigmoid(0.16844) ~= 0.5420
        let p = PlattScaler { intercept: 0.16844, slope: 1.56331, clip_eps: 1e-4 };
        let out = p.calibrate(0.5);
        assert!((out - 0.5420).abs() < 1e-3, "out={out}");
    }

    /// Build a 2-node tree by hand and check predict() obeys `<= left`.
    #[test]
    fn tiny_tree_traversal() {
        // root (split feat 0 at 1.5) -> left=leaf(0.1), right=leaf(0.9)
        let tree = Tree {
            split_var: vec![0, -1, -1],
            split_val: vec![1.5, 0.0, 0.0],
            left:      vec![1,  0,  0],
            right:     vec![2,  0,  0],
            prediction: vec![0.0, 0.1, 0.9],
        };
        assert_eq!(tree.predict(&[1.0]), 0.1);
        assert_eq!(tree.predict(&[1.5]), 0.1); // <= goes left
        assert_eq!(tree.predict(&[1.6]), 0.9);
    }
}

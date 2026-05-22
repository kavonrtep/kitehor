//! Pipeline merger — port of `tools/rule_proto/summary.py`.
//!
//! Joins verdicts.tsv + subrepeat.tsv + ssr.tsv (+ optional
//! hor_within_tile.tsv) on `record_id` (verdicts uses `case_id`, renamed
//! here). Applies the eight-rule first-match-wins decision tree to
//! produce `combined_class`. Output is one row per record with 28–32
//! columns at `%.4g` float precision.

pub mod io;

use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// `pure_ssr` fires when `ssr_dominant_motif_coverage_pct ≥` this.
    pub pure_ssr_pct_threshold: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            pure_ssr_pct_threshold: 80.0,
        }
    }
}

/// Apply the 8-rule decision tree (first-match-wins). Mirrors
/// `summary.py::combined_class`.
pub fn combined_class(
    hor_verdict: &str,
    ssr_flag: &str,
    ssr_dom_pct: f64,
    subrep_flag: &str,
    density_hint: &str,
    cfg: &Config,
) -> &'static str {
    if ssr_flag == "yes" && ssr_dom_pct >= cfg.pure_ssr_pct_threshold {
        return "pure_ssr";
    }
    if subrep_flag == "yes" {
        return "tr_with_nested_tr";
    }
    if hor_verdict == "hor" && density_hint == "localized_duplication" {
        return "tr_with_subrepeat";
    }
    if hor_verdict == "hor" && ssr_flag == "yes" {
        return "hor_with_ssr";
    }
    if hor_verdict == "hor" {
        return "hor";
    }
    if hor_verdict == "simple_tr" && ssr_flag == "yes" {
        return "tr_with_ssr";
    }
    if hor_verdict == "simple_tr" {
        return "tr";
    }
    "unresolved"
}

/// Subcommand entry point: read inputs, merge, write `<prefix>.summary.tsv`.
pub fn run_subcommand(
    verdicts: &Path,
    subrepeat: &Path,
    ssr: &Path,
    within_tile: Option<&Path>,
    out_prefix: &Path,
    cfg: &Config,
) -> Result<usize> {
    let merged = io::merge_inputs(verdicts, subrepeat, ssr, within_tile, cfg)?;
    let n = merged.rows.len();
    let path = io::summary_path(out_prefix);
    io::write_summary(&path, &merged)?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_ssr_fires_above_threshold() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("hor", "yes", 85.0, "no", "", &cfg),
            "pure_ssr"
        );
    }

    #[test]
    fn nested_tr_beats_hor() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("hor", "no", 0.0, "yes", "spatially_confirms_hor", &cfg),
            "tr_with_nested_tr"
        );
    }

    #[test]
    fn tr_with_subrepeat_via_within_tile() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("hor", "no", 0.0, "no", "localized_duplication", &cfg),
            "tr_with_subrepeat"
        );
    }

    #[test]
    fn hor_with_ssr_below_pure_threshold() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("hor", "yes", 50.0, "no", "", &cfg),
            "hor_with_ssr"
        );
    }

    #[test]
    fn plain_hor() {
        let cfg = Config::default();
        assert_eq!(combined_class("hor", "no", 0.0, "no", "", &cfg), "hor");
    }

    #[test]
    fn tr_with_ssr() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("simple_tr", "yes", 30.0, "no", "", &cfg),
            "tr_with_ssr"
        );
    }

    #[test]
    fn plain_tr() {
        let cfg = Config::default();
        assert_eq!(combined_class("simple_tr", "no", 0.0, "no", "", &cfg), "tr");
    }

    #[test]
    fn unresolved_default() {
        let cfg = Config::default();
        assert_eq!(
            combined_class("unresolved", "no", 0.0, "no", "", &cfg),
            "unresolved"
        );
    }
}

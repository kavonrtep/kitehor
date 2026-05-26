//! Pipeline merger — port of `tools/rule_proto/summary_unified.py`.
//!
//! Joins verdicts.tsv + tandem_validate.tsv + ssr.tsv on `record_id`
//! (verdicts uses `case_id`, renamed here). Applies the nine-rule
//! first-match-wins decision tree to produce `combined_class`. Output
//! is one row per record with the merged column set at `%.4g` float
//! precision.
//!
//! v0.11 cascade changes (vs. v0.10):
//! * SSR decisions are driven by the **array-scale**
//!   `ssr_raw_total_coverage_pct` (computed by the raw scanner against
//!   `seq.len()`), not by the consensus-path-overridden
//!   `dominant_motif_coverage_pct` (which under `consensus_single`
//!   reports the candidate monomer's *self*-coverage, not the array's).
//! * Two new classes added so every verdict category gets a parallel
//!   `<base>_with_ssr` partner when SSR coverage clears the
//!   "has-SSR" floor: `unresolved_with_ssr` and
//!   `tr_with_subrepeat_with_ssr`. Total: 9 classes.

pub mod io;

use anyhow::Result;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// `pure_ssr` fires when `ssr_raw_total_coverage_pct ≥` this.
    /// Default 80.0.
    pub pure_ssr_pct_threshold: f64,
    /// `_with_ssr` partner classes fire when
    /// `ssr_raw_total_coverage_pct ≥` this (but below the pure
    /// threshold). Default 30.0 — matches `ssr::Config::ssr_flag_threshold_pct`
    /// so the recomputed `ssr_flag` column tracks this cascade.
    pub ssr_has_pct_threshold: f64,
    /// Minimum `tv_density` (fraction of windows in which the candidate
    /// subrepeat motif is present) required to promote a
    /// `localized_subrepeat` decision to `tr_with_subrepeat` /
    /// `tr_with_subrepeat_with_ssr`. Below this, the record falls
    /// through to its natural class (`tr` / `tr_with_ssr` / etc.) —
    /// rejecting indel / localized-duplication artifacts that show
    /// short-period signal in only a few windows. Default 0.7.
    pub subrepeat_density_min: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            pure_ssr_pct_threshold: 80.0,
            ssr_has_pct_threshold: 30.0,
            subrepeat_density_min: 0.7,
        }
    }
}

/// Apply the 9-rule decision tree (first-match-wins). Mirrors
/// `summary_unified.py::combined_class` (v0.11 update).
///
/// * `ssr_raw_total_pct` is the array-scale total-SSR coverage
///   (`ssr.tsv::raw_total_coverage_pct`, summary-renamed to
///   `ssr_raw_total_coverage_pct`). Never use the consensus-path's
///   `ssr_dominant_motif_coverage_pct` — under `consensus_single` it
///   reflects the candidate monomer's self-coverage (≈100%), not the
///   array's, and would over-trigger `pure_ssr`.
/// * `tv_decision` is the `decision_hint` column from `tandem_validate.tsv`.
///   Only `localized_subrepeat` paired with `tv_density >=
///   subrepeat_density_min` keeps a row out of its verdict's natural
///   class. Low-density "localized" hits are presumed indel /
///   short-event artifacts and fall through to `tr` / `unresolved`.
/// * `tv_density` is the fraction of scanner windows in which the
///   candidate subrepeat motif was present (column `tv_density`).
///   `None` (missing/empty) is treated as failing the threshold.
pub fn combined_class(
    hor_verdict: &str,
    ssr_raw_total_pct: f64,
    tv_decision: &str,
    tv_density: Option<f64>,
    cfg: &Config,
) -> &'static str {
    let has_ssr = ssr_raw_total_pct >= cfg.ssr_has_pct_threshold;
    if ssr_raw_total_pct >= cfg.pure_ssr_pct_threshold {
        return "pure_ssr";
    }
    if tv_decision == "localized_subrepeat"
        && tv_density.is_some_and(|d| d >= cfg.subrepeat_density_min)
    {
        return if has_ssr {
            "tr_with_subrepeat_with_ssr"
        } else {
            "tr_with_subrepeat"
        };
    }
    if hor_verdict == "hor" {
        return if has_ssr { "hor_with_ssr" } else { "hor" };
    }
    if hor_verdict == "simple_tr" {
        return if has_ssr { "tr_with_ssr" } else { "tr" };
    }
    if has_ssr {
        "unresolved_with_ssr"
    } else {
        "unresolved"
    }
}

/// Subcommand entry point: read inputs, merge, write `<prefix>.summary.tsv`.
pub fn run_subcommand(
    verdicts: &Path,
    tandem_validate: &Path,
    ssr: &Path,
    out_prefix: &Path,
    cfg: &Config,
) -> Result<usize> {
    let merged = io::merge_inputs(verdicts, tandem_validate, ssr, cfg)?;
    let n = merged.rows.len();
    let path = io::summary_path(out_prefix);
    io::write_summary(&path, &merged)?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> Config {
        Config::default()
    }

    /// Helper: density 1.0 always clears the subrepeat-density gate.
    fn d_full() -> Option<f64> {
        Some(1.0)
    }

    // --- pure_ssr threshold ---

    #[test]
    fn pure_ssr_fires_above_threshold() {
        assert_eq!(
            combined_class("hor", 85.0, "no_candidates", None, &cfg()),
            "pure_ssr"
        );
    }

    #[test]
    fn pure_ssr_outranks_subrepeat() {
        // High raw-total SSR coverage wins over a positive
        // tandem_validate localized_subrepeat hint.
        assert_eq!(
            combined_class("hor", 90.0, "localized_subrepeat", d_full(), &cfg()),
            "pure_ssr"
        );
    }

    #[test]
    fn pure_ssr_not_triggered_by_consensus_inflated_dom_cov() {
        // Regression: under v0.10 the consensus_single path inflated
        // dominant_motif_coverage_pct to ≈100 even when the array
        // was a few % SSR. v0.11 uses raw_total_pct, so a 4.21%
        // array (e.g. TRC_18:chr7_2164485_2172903) no longer fires
        // pure_ssr.
        assert_ne!(
            combined_class("simple_tr", 4.21, "no_candidates", None, &cfg()),
            "pure_ssr"
        );
    }

    // --- _with_subrepeat family + density gate ---

    #[test]
    fn subrepeat_without_ssr_above_density() {
        assert_eq!(
            combined_class("hor", 0.0, "localized_subrepeat", d_full(), &cfg()),
            "tr_with_subrepeat"
        );
    }

    #[test]
    fn subrepeat_with_ssr_above_density() {
        assert_eq!(
            combined_class("simple_tr", 35.0, "localized_subrepeat", d_full(), &cfg()),
            "tr_with_subrepeat_with_ssr"
        );
    }

    #[test]
    fn subrepeat_below_density_falls_to_tr() {
        // v0.12: localized_subrepeat with density below the gate is
        // treated as a likely indel/duplication artifact and falls
        // through to the natural class.
        assert_eq!(
            combined_class("simple_tr", 0.0, "localized_subrepeat", Some(0.2), &cfg()),
            "tr"
        );
    }

    #[test]
    fn subrepeat_below_density_with_ssr_falls_to_tr_with_ssr() {
        assert_eq!(
            combined_class("simple_tr", 35.0, "localized_subrepeat", Some(0.2), &cfg()),
            "tr_with_ssr"
        );
    }

    #[test]
    fn subrepeat_missing_density_does_not_fire() {
        // Missing tv_density → safer default: treat as below threshold.
        assert_eq!(
            combined_class("simple_tr", 0.0, "localized_subrepeat", None, &cfg()),
            "tr"
        );
    }

    #[test]
    fn subrepeat_at_threshold_fires() {
        // Density exactly == threshold should clear the >= gate.
        let c = cfg();
        assert_eq!(
            combined_class(
                "simple_tr",
                0.0,
                "localized_subrepeat",
                Some(c.subrepeat_density_min),
                &c,
            ),
            "tr_with_subrepeat"
        );
    }

    #[test]
    fn confirms_host_falls_through_to_hor() {
        // tv = confirms_host → cascade falls through to verdict.
        assert_eq!(
            combined_class("hor", 0.0, "confirms_host", None, &cfg()),
            "hor"
        );
    }

    // --- _with_ssr partners (parallel to base) ---

    #[test]
    fn hor_with_ssr_below_pure_threshold() {
        assert_eq!(
            combined_class("hor", 50.0, "no_candidates", None, &cfg()),
            "hor_with_ssr"
        );
    }

    #[test]
    fn plain_hor() {
        assert_eq!(
            combined_class("hor", 0.0, "no_candidates", None, &cfg()),
            "hor"
        );
    }

    #[test]
    fn tr_with_ssr() {
        assert_eq!(
            combined_class("simple_tr", 35.0, "no_candidates", None, &cfg()),
            "tr_with_ssr"
        );
    }

    #[test]
    fn plain_tr() {
        assert_eq!(
            combined_class("simple_tr", 0.0, "no_candidates", None, &cfg()),
            "tr"
        );
    }

    #[test]
    fn unresolved_with_ssr() {
        // New in v0.11. Parallel to hor_with_ssr / tr_with_ssr.
        assert_eq!(
            combined_class("unresolved", 35.0, "no_candidates", None, &cfg()),
            "unresolved_with_ssr"
        );
    }

    #[test]
    fn unresolved_default() {
        assert_eq!(
            combined_class("unresolved", 0.0, "no_candidates", None, &cfg()),
            "unresolved"
        );
    }

    // --- skip / pass-through cases ---

    #[test]
    fn skip_k2_falls_through_to_hor() {
        assert_eq!(combined_class("hor", 0.0, "skip_k2", None, &cfg()), "hor");
    }

    #[test]
    fn ambiguous_falls_through() {
        assert_eq!(combined_class("hor", 0.0, "ambiguous", None, &cfg()), "hor");
        assert_eq!(
            combined_class("simple_tr", 0.0, "ambiguous", None, &cfg()),
            "tr"
        );
    }

    // --- threshold boundaries ---

    #[test]
    fn has_ssr_threshold_inclusive() {
        // Exactly at the threshold → "has SSR" fires.
        assert_eq!(
            combined_class("simple_tr", 30.0, "no_candidates", None, &cfg()),
            "tr_with_ssr"
        );
    }

    #[test]
    fn pure_ssr_threshold_inclusive() {
        assert_eq!(
            combined_class("simple_tr", 80.0, "no_candidates", None, &cfg()),
            "pure_ssr"
        );
    }

    #[test]
    fn just_below_has_ssr_threshold() {
        assert_eq!(
            combined_class("simple_tr", 29.9, "no_candidates", None, &cfg()),
            "tr"
        );
    }
}

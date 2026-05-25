//! Pipeline merger — port of `tools/rule_proto/summary_unified.py`.
//!
//! Joins verdicts.tsv + tandem_validate.tsv + ssr.tsv on `record_id`
//! (verdicts uses `case_id`, renamed here). Applies the seven-rule
//! first-match-wins decision tree to produce `combined_class`. Output
//! is one row per record with the merged column set at `%.4g` float
//! precision.
//!
//! Cascade compared to the prior 8-class merger:
//! * `tr_with_nested_tr` is retired (merged into `tr_with_subrepeat`).
//! * The two subrepeat-style triggers (subrepeat-scan's `blocks+non_sub`
//!   flag and hor-validate's `localized_duplication` hint) collapse to
//!   a single `tandem_validate decision == localized_subrepeat`.

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

/// Apply the 7-rule decision tree (first-match-wins). Mirrors
/// `summary_unified.py::combined_class`.
///
/// `tv_decision` is the `decision_hint` column from `tandem_validate.tsv`.
/// Only `localized_subrepeat` drives the `tr_with_subrepeat` class; any
/// other value (`confirms_host`, `ambiguous`, `no_signal`, `no_candidates`,
/// `no_windows`, `skip_k2`, `no_host`, `no_verdict`) lets the cascade
/// fall through to the verdict's natural class.
pub fn combined_class(
    hor_verdict: &str,
    ssr_flag: &str,
    ssr_dom_pct: f64,
    tv_decision: &str,
    cfg: &Config,
) -> &'static str {
    if ssr_flag == "yes" && ssr_dom_pct >= cfg.pure_ssr_pct_threshold {
        return "pure_ssr";
    }
    if tv_decision == "localized_subrepeat" {
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

    #[test]
    fn pure_ssr_fires_above_threshold() {
        assert_eq!(
            combined_class("hor", "yes", 85.0, "no_candidates", &cfg()),
            "pure_ssr"
        );
    }

    #[test]
    fn pure_ssr_outranks_subrepeat() {
        // A record with both high SSR coverage AND a localized subrepeat
        // signal is classified pure_ssr — the SSR rule comes first.
        assert_eq!(
            combined_class("hor", "yes", 90.0, "localized_subrepeat", &cfg()),
            "pure_ssr"
        );
    }

    #[test]
    fn subrepeat_beats_hor() {
        assert_eq!(
            combined_class("hor", "no", 0.0, "localized_subrepeat", &cfg()),
            "tr_with_subrepeat"
        );
    }

    #[test]
    fn subrepeat_beats_simple_tr() {
        assert_eq!(
            combined_class("simple_tr", "no", 0.0, "localized_subrepeat", &cfg()),
            "tr_with_subrepeat"
        );
    }

    #[test]
    fn confirms_host_falls_through_to_hor() {
        // The unified detector says the tile is confirmed — fall through
        // to the verdict's natural class.
        assert_eq!(
            combined_class("hor", "no", 0.0, "confirms_host", &cfg()),
            "hor"
        );
    }

    #[test]
    fn hor_with_ssr_below_pure_threshold() {
        assert_eq!(
            combined_class("hor", "yes", 50.0, "no_candidates", &cfg()),
            "hor_with_ssr"
        );
    }

    #[test]
    fn plain_hor() {
        assert_eq!(
            combined_class("hor", "no", 0.0, "no_candidates", &cfg()),
            "hor"
        );
    }

    #[test]
    fn tr_with_ssr() {
        assert_eq!(
            combined_class("simple_tr", "yes", 30.0, "no_candidates", &cfg()),
            "tr_with_ssr"
        );
    }

    #[test]
    fn plain_tr() {
        assert_eq!(
            combined_class("simple_tr", "no", 0.0, "no_candidates", &cfg()),
            "tr"
        );
    }

    #[test]
    fn skip_k2_falls_through_to_hor() {
        // HOR k=2 records skip the detector entirely; cascade should
        // still call them hor.
        assert_eq!(combined_class("hor", "no", 0.0, "skip_k2", &cfg()), "hor");
    }

    #[test]
    fn ambiguous_falls_through() {
        assert_eq!(combined_class("hor", "no", 0.0, "ambiguous", &cfg()), "hor");
        assert_eq!(
            combined_class("simple_tr", "no", 0.0, "ambiguous", &cfg()),
            "tr"
        );
    }

    #[test]
    fn unresolved_default() {
        assert_eq!(
            combined_class("unresolved", "no", 0.0, "no_candidates", &cfg()),
            "unresolved"
        );
    }
}

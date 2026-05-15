//! Verdict orchestrator: features + RFs + Platt scaling + 4-category
//! decision + family demotion + k-recovery.
//!
//! Port of `eval/training_data/predict_verdict.R`. Given a built
//! [`FeatureRow`], the calibrated thresholds, and both forests (HOR
//! probability + k regression), this returns a [`Verdict`] explaining
//! the model's call.
//!
//! Decision order:
//!   1. s1 == 0 || d1 == 0                  -> NoSignal
//!   2. score >= t_high                      -> Hor          (subject to demote)
//!   3. score <  t_low                       -> Tandem
//!   4. otherwise                            -> Unresolved
//!
//! Demote-to-unresolved: if RF says HOR but the kite family search
//! found no real founder (founder_d == 0 or founder_d == tile_d),
//! we have nothing to report → flip to Unresolved unless k-recovery
//! succeeds.
//!
//! k-recovery (Option B from MISTAKE_TRIAGE §16): for demoted HOR
//! cases, ask the k-predictor for an integer k. If `d1 * k` or `d1 / k`
//! matches one of `d2 .. d_{top_k+1}` within tolerance, accept the
//! recovered (founder, k, tile) tuple.

use crate::classifier::{ClassifierConfig, PlattScaler, RandomForest};
use crate::features::FeatureRow;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerdictCategory {
    Hor,
    Tandem,
    Unresolved,
    NoSignal,
}

impl VerdictCategory {
    pub fn as_str(self) -> &'static str {
        match self {
            VerdictCategory::Hor => "hor",
            VerdictCategory::Tandem => "tandem",
            VerdictCategory::Unresolved => "unresolved",
            VerdictCategory::NoSignal => "no_signal",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Verdict {
    pub category: VerdictCategory,
    /// Platt-calibrated HOR probability in [0, 1].
    pub hor_score: f64,
    /// Raw (uncalibrated) RF probability in [0, 1].
    pub hor_score_raw: f64,
    /// Founder period (bp). None when not applicable.
    pub founder: Option<usize>,
    /// Multiplicity (k). None when not applicable.
    pub multiplicity: Option<usize>,
    /// Tile period (bp). None when not applicable.
    pub tile: Option<usize>,
    /// Integer k from the k-predictor (set whenever recovery is
    /// attempted, even if recovery failed). None when k-recovery
    /// was not run.
    pub k_pred: Option<usize>,
    /// True iff this HOR call comes from k-recovery (i.e., the kite
    /// family search alone would have demoted it).
    pub recovered: bool,
}

/// Apply imputation, evaluate the forests + Platt, run the verdict
/// state machine. Mutates the feature row in-place to record the
/// imputed values (`h_d1`, `h_founder`) so callers can inspect them.
pub fn classify(
    row: &mut FeatureRow,
    cfg: &ClassifierConfig,
    platt: &PlattScaler,
    hor_model: &RandomForest,
    k_model: Option<&RandomForest>,
) -> Verdict {
    // 1. Imputation. The R pipeline imputes h_d1 / h_founder with
    //    training medians when probe-periods returned NA.
    if row.h_d1.is_nan() {
        row.h_d1 = cfg.imputation.h_d1;
    }
    if row.h_founder.is_nan() {
        row.h_founder = cfg.imputation.h_founder;
    }

    // 2. HOR-probability RF.
    let x_h = assemble(hor_model, row);
    let raw = hor_model.predict(&x_h);
    let cal = platt.calibrate(raw);

    // 3. Provisional category.
    let no_signal = row.s1 == 0.0 || row.d1 == 0;
    let mut category = if no_signal {
        VerdictCategory::NoSignal
    } else if cal >= cfg.thresholds.t_high {
        VerdictCategory::Hor
    } else if cal < cfg.thresholds.t_low {
        VerdictCategory::Tandem
    } else {
        VerdictCategory::Unresolved
    };

    // 4. Family demotion: HOR-call needs a real (founder, tile) pair
    //    from the kite family search.
    let no_family = row.family_founder_d == 0
        || row.family_tile_d == 0
        || row.family_founder_d == row.family_tile_d;
    let mut founder: Option<usize> = None;
    let mut tile: Option<usize> = None;
    let mut multiplicity: Option<usize> = None;
    let mut k_pred: Option<usize> = None;
    let mut recovered = false;

    if category == VerdictCategory::Hor && !no_family {
        founder = Some(row.family_founder_d);
        tile = Some(row.family_tile_d);
        multiplicity = Some(
            ((row.family_tile_d as f64 / row.family_founder_d as f64).round() as usize).max(1),
        );
    } else if category == VerdictCategory::Hor && no_family {
        // Try k-recovery if we have the k-model.
        if let Some(km) = k_model {
            let x_k = assemble(km, row);
            let mut k_int = km.predict(&x_k).round() as i64;
            if k_int < cfg.recovery.k_min {
                k_int = cfg.recovery.k_min;
            }
            let k_int = k_int as usize;
            k_pred = Some(k_int);

            let d1 = row.d1;
            if d1 > 0 {
                let candidates = top_k_candidates(row, cfg.recovery.top_k_candidates);
                // Hypothesis A: d1 is the HOR-unit (tile), founder = d1 / k.
                let cand_founder_a = (d1 as f64 / k_int as f64).round() as i64;
                let matched_a = match_period(
                    cand_founder_a,
                    &candidates,
                    cfg.recovery.tol_bp,
                    cfg.recovery.tol_rel,
                );
                // Hypothesis B: d1 is the founder, tile = d1 * k.
                let cand_tile_b = (d1 as i64) * (k_int as i64);
                let matched_b = match_period(
                    cand_tile_b,
                    &candidates,
                    cfg.recovery.tol_bp,
                    cfg.recovery.tol_rel,
                );

                if matched_a.is_some()
                    && cand_founder_a >= cfg.recovery.founder_min_bp
                {
                    let founder_bp = matched_a.unwrap() as usize;
                    row.family_founder_d = founder_bp;
                    row.family_tile_d = d1;
                    founder = Some(founder_bp);
                    tile = Some(d1);
                    multiplicity = Some(((d1 as f64 / founder_bp as f64).round() as usize).max(1));
                    recovered = true;
                } else if matched_b.is_some() && cand_tile_b > d1 as i64 {
                    let tile_bp = matched_b.unwrap() as usize;
                    row.family_founder_d = d1;
                    row.family_tile_d = tile_bp;
                    founder = Some(d1);
                    tile = Some(tile_bp);
                    multiplicity = Some(((tile_bp as f64 / d1 as f64).round() as usize).max(1));
                    recovered = true;
                } else {
                    // Recovery failed → demote.
                    category = VerdictCategory::Unresolved;
                }
            } else {
                category = VerdictCategory::Unresolved;
            }
        } else {
            // No k-model → demote.
            category = VerdictCategory::Unresolved;
        }
    } else if category == VerdictCategory::Tandem {
        // R convention: tandem reports only tile + multiplicity=1;
        // founder stays NA (see predict_verdict.R, the tan_idx branch).
        tile = if row.d1 > 0 { Some(row.d1) } else { None };
        multiplicity = Some(1);
    }

    Verdict {
        category,
        hor_score: cal,
        hor_score_raw: raw,
        founder,
        multiplicity,
        tile,
        k_pred,
        recovered,
    }
}

fn assemble(model: &RandomForest, row: &FeatureRow) -> Vec<f64> {
    model
        .feature_names
        .iter()
        .map(|name| row.get(name.as_str()).unwrap_or(f64::NAN))
        .collect()
}

/// Pull candidates from `d2..d_{top_k}` — i.e., the top-`top_k` peaks
/// excluding `d1` (which is the suspected founder/tile in the recovery
/// hypotheses). Matches R: `candidates <- c(d$d2[i], d$d3[i], d$d4[i],
/// d$d5[i])` for top_k=5.
fn top_k_candidates(row: &FeatureRow, top_k: usize) -> Vec<i64> {
    let all = [
        row.d2, row.d3, row.d4, row.d5, row.d6, row.d7, row.d8, row.d9, row.d10,
    ];
    let n = top_k.saturating_sub(1);
    all.iter()
        .take(n)
        .filter(|&&p| p > 0)
        .map(|&p| p as i64)
        .collect()
}

/// Return the first candidate within ±tol_bp / ±tol_rel of the target,
/// matching the R `match_period()`. Returns None if no candidate fits.
fn match_period(target: i64, candidates: &[i64], tol_bp: i64, tol_rel: f64) -> Option<i64> {
    for &c in candidates {
        if c <= 0 {
            continue;
        }
        let diff = (c - target).abs();
        let abs_tol = tol_bp.max((tol_rel * c.max(target) as f64) as i64);
        if diff <= abs_tol {
            return Some(c);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synth_row() -> FeatureRow {
        FeatureRow {
            case_id: "x".into(),
            stratum: "x".into(),
            array_length: 100000,
            s1: 0.5,
            s2: 0.2,
            s3: 0.1,
            s2_over_s1: 0.4,
            s3_over_s1: 0.2,
            family_size_best: 3,
            tile_founder_ratio: 0.4,
            tile_jitter: 1,
            d1: 178,
            d2: 356,
            d3: 534,
            log_d1_over_l: -6.3,
            d2_over_d1: 2.0,
            d3_over_d1: 3.0,
            max_d_top3_over_min_d_top3: 3.0,
            d4: 712,
            d5: 0,
            d6: 0,
            d7: 0,
            d8: 0,
            d9: 0,
            d10: 0,
            family_founder_d: 178,
            family_tile_d: 534,
            distinct_kmers_per_bp: 0.001,
            kmer_entropy: 5.0,
            singletons_ratio: 0.01,
            h_d1: 0.95,
            h_founder: 0.65,
        }
    }

    #[test]
    fn no_signal_when_kite_is_silent() {
        // Make a stub RF that returns 0.0 always (5 features).
        // For this test we use the real model via load_json; if files
        // aren't available the test is skipped. (Run from crate dir.)
        let cfg = ClassifierConfig::default_baked().unwrap();
        let h_path = "models/hor_score.rftrees.json";
        if !std::path::Path::new(h_path).exists() {
            eprintln!("skipping: model artifact missing (run from crate dir)");
            return;
        }
        let h = RandomForest::load_json(h_path).unwrap();
        let mut r = synth_row();
        r.s1 = 0.0;
        r.d1 = 0;
        let v = classify(&mut r, &cfg, &cfg.platt(), &h, None);
        assert_eq!(v.category, VerdictCategory::NoSignal);
    }

    #[test]
    fn match_period_within_relative_tol() {
        // target 534, candidates [356, 530] → 530 is within 0.02*534=10.68
        // (and tol_bp=5 gives 5, but 530 within ±5 of 534 already).
        let m = match_period(534, &[356, 530], 5, 0.02);
        assert_eq!(m, Some(530));
    }

    #[test]
    fn match_period_no_match() {
        let m = match_period(100, &[200, 300, 400], 5, 0.02);
        assert!(m.is_none());
    }
}

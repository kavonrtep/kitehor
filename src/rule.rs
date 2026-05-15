//! Rule-based HOR classifier.
//!
//! A four-condition rule applied directly to the kite peak output. Every
//! peak in `KiteResult.peaks` has already passed kite's
//! `peak > background` and `score2_norm > 0.001` filters, so we trust
//! it as a real periodicity and only ask whether the peak structure
//! forms a Higher-Order Repeat.
//!
//! ```text
//! HOR ⟺ d1 (strongest kite peak) = k × p_n
//!       for some k ∈ [2, qmax]
//!       and p_n ∈ top-N peaks by score
//!       and |d1 − k·p_n| ≤ max(tol_bp, tol_rel × k·p_n)
//!       and d1 ≥ lo_period
//!       and p_n ≥ lo_period
//! ```
//!
//! Unidirectional: `d1` must be the *tile*, `p_n` the *founder*. The
//! converse arrangement (`d1` = founder, tile at some `p_n > d1`) is
//! the typical synthetic-data shape but a poor predictor on real
//! centromeric arrays — see `docs/rule.md` for the empirical study
//! that led to this design.

use crate::kite::{KitePeak, KiteResult};

#[derive(Debug, Clone, Copy)]
pub struct RuleConfig {
    /// Founder candidate must be in the top-`top_n` kite peaks by
    /// score. Filtering deeper sub-period harmonics out of the
    /// founder-candidate pool is the main precision lever. Default 3.
    pub top_n: usize,
    /// Maximum multiplicity considered. Default 30.
    pub qmax: usize,
    /// Absolute period-match tolerance (bp). Default 5.
    pub tol_bp: usize,
    /// Relative period-match tolerance. Default 0.02.
    pub tol_rel: f64,
    /// Minimum period (bp) for either tile or founder candidate.
    /// Default 15.
    pub lo_period: usize,
}

impl Default for RuleConfig {
    fn default() -> Self {
        Self {
            top_n: 3,
            qmax: 30,
            tol_bp: 5,
            tol_rel: 0.02,
            lo_period: 15,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RuleVerdict {
    /// No kite peaks survived filtering.
    NoSignal,
    /// d1 < lo_period — too short to be a meaningful repeat.
    Unresolved,
    /// Single dominant peak with no integer-multiple partner.
    Tandem { monomer_bp: usize },
    /// Tile (= d1) is a k-fold multiple of a top-N founder peak.
    Hor {
        founder: usize,
        tile: usize,
        k: usize,
        share: f64,
    },
}

impl RuleVerdict {
    pub fn as_str(self) -> &'static str {
        match self {
            RuleVerdict::Hor { .. } => "hor",
            RuleVerdict::Tandem { .. } => "tandem",
            RuleVerdict::Unresolved => "unresolved",
            RuleVerdict::NoSignal => "no_signal",
        }
    }

    pub fn founder(&self) -> Option<usize> {
        match self {
            RuleVerdict::Hor { founder, .. } => Some(*founder),
            _ => None,
        }
    }

    pub fn tile(&self) -> Option<usize> {
        match self {
            RuleVerdict::Hor { tile, .. } => Some(*tile),
            RuleVerdict::Tandem { monomer_bp } => Some(*monomer_bp),
            _ => None,
        }
    }

    pub fn multiplicity(&self) -> Option<usize> {
        match self {
            RuleVerdict::Hor { k, .. } => Some(*k),
            RuleVerdict::Tandem { .. } => Some(1),
            _ => None,
        }
    }

    pub fn share(&self) -> Option<f64> {
        match self {
            RuleVerdict::Hor { share, .. } => Some(*share),
            _ => None,
        }
    }
}

/// Apply the kite rule to a record's peak list.
pub fn classify_peaks(peaks: &[KitePeak], cfg: &RuleConfig) -> RuleVerdict {
    let Some(d1) = peaks.first() else {
        return RuleVerdict::NoSignal;
    };
    if d1.period < cfg.lo_period {
        return RuleVerdict::Unresolved;
    }

    let end = peaks.len().min(cfg.top_n);
    let mut best: Option<(usize, usize, f64)> = None; // (founder, k, share)
    for p in peaks[1..end].iter() {
        if p.period < cfg.lo_period || p.period >= d1.period {
            continue;
        }
        let k_f = (d1.period as f64 / p.period as f64).round();
        if !(2.0..=cfg.qmax as f64).contains(&k_f) {
            continue;
        }
        let k = k_f as usize;
        let expected = k * p.period;
        let diff = d1.period.abs_diff(expected);
        let tol = cfg.tol_bp.max((cfg.tol_rel * expected as f64) as usize);
        if diff > tol {
            continue;
        }
        let share = if p.score > 0.0 && d1.score > 0.0 {
            p.score.min(d1.score) / p.score.max(d1.score)
        } else {
            0.0
        };
        if best.as_ref().map(|(_, _, s)| share > *s).unwrap_or(true) {
            best = Some((p.period, k, share));
        }
    }

    if let Some((founder, k, share)) = best {
        RuleVerdict::Hor {
            founder,
            tile: d1.period,
            k,
            share,
        }
    } else {
        RuleVerdict::Tandem {
            monomer_bp: d1.period,
        }
    }
}

/// Convenience wrapper that takes the full `KiteResult`.
pub fn classify(kite: &KiteResult, cfg: &RuleConfig) -> RuleVerdict {
    classify_peaks(&kite.peaks, cfg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kite::KitePeak;

    fn mk(period: usize, score: f64) -> KitePeak {
        KitePeak {
            period,
            peak_height: 1.0,
            score,
            score2: 0.0,
            score2_norm: 0.01,
            background: 0.0,
        }
    }

    #[test]
    fn empty_peaks_no_signal() {
        let v = classify_peaks(&[], &RuleConfig::default());
        assert_eq!(v, RuleVerdict::NoSignal);
    }

    #[test]
    fn single_peak_is_tandem() {
        let pks = vec![mk(178, 0.5)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        assert_eq!(v, RuleVerdict::Tandem { monomer_bp: 178 });
    }

    #[test]
    fn clean_k2_hor_fires() {
        // d1 = 744 (tile), d2 = 372 (founder), score share 0.19 (below
        // any historical share floor — confirms we trust kite peaks).
        let pks = vec![mk(744, 0.44), mk(372, 0.083), mk(168, 0.04)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        match v {
            RuleVerdict::Hor {
                founder, tile, k, ..
            } => {
                assert_eq!(founder, 372);
                assert_eq!(tile, 744);
                assert_eq!(k, 2);
            }
            other => panic!("expected Hor, got {other:?}"),
        }
    }

    #[test]
    fn founder_dominant_pattern_is_tandem() {
        // d1=178 (the AT178 monomer, founder-dominant), d2=356, d3=168.
        // No top-3 peak is a divisor of 178 with k>=2, so unidirectional
        // rule says tandem.
        let pks = vec![mk(178, 0.56), mk(356, 0.10), mk(168, 0.03)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        assert!(matches!(v, RuleVerdict::Tandem { monomer_bp: 178 }));
    }

    #[test]
    fn k5_hor_fires() {
        // Real TRC_1 case: d1=888, d2=178, founder=178 in top-3.
        let pks = vec![mk(888, 0.13), mk(178, 0.12), mk(356, 0.07)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        match v {
            RuleVerdict::Hor {
                founder, tile, k, ..
            } => {
                assert_eq!(founder, 178);
                assert_eq!(tile, 888);
                assert_eq!(k, 5);
            }
            other => panic!("expected Hor k=5, got {other:?}"),
        }
    }

    #[test]
    fn deep_subperiod_not_promoted() {
        // d1=310 dominant, with a deep weak peak at 51 that happens to
        // satisfy 310 = 6 * 51. With top_n=3, 51 is not a candidate.
        let pks = vec![
            mk(310, 0.40),
            mk(10687, 0.10),
            mk(10689, 0.09),
            mk(51, 0.18),
        ];
        let v = classify_peaks(&pks, &RuleConfig::default());
        // d2=10687, d3=10689 are larger than d1 -> direction A rejects.
        // 51 is at index 3, outside top_n=3, so doesn't qualify.
        assert!(matches!(v, RuleVerdict::Tandem { .. }));
    }

    #[test]
    fn lo_period_blocks_tiny_d1() {
        let pks = vec![mk(6, 0.5), mk(3, 0.2)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        assert_eq!(v, RuleVerdict::Unresolved);
    }

    #[test]
    fn tolerance_absorbs_small_drift() {
        // 744 vs 2*372 = 744, exact. Now try 741 = 2*372 (diff 3, tol 14).
        let pks = vec![mk(741, 0.44), mk(372, 0.08), mk(168, 0.04)];
        let v = classify_peaks(&pks, &RuleConfig::default());
        match v {
            RuleVerdict::Hor {
                founder, tile, k, ..
            } => {
                assert_eq!(founder, 372);
                assert_eq!(tile, 741);
                assert_eq!(k, 2);
            }
            _ => panic!("expected Hor"),
        }
    }
}

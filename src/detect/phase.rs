//! Phase separation + primitive multiplicity correction
//! (`detect_impl_plan.md §6.7`).
//!
//! ```text
//! best_k          = argmax_{k ≥ 2} R(k)
//! same_phase      = R(best_k)
//! background      = median R(k') for k' near best_k but not a multiple
//! phase_separation = same_phase − background
//! ```
//!
//! Primitive correction: if `R(d) ≥ R(best_k) − δ` for some divisor
//! `d ≥ 2` of `best_k`, prefer the smaller `d`. Avoids reporting
//! `k = 12` for arrays whose primitive multiplicity is 6 (or 3, or
//! 4) just because R(12) ≥ R(6) by a tiny margin.

use crate::detect::widths::divisors;

/// `phase_separation(R, best_k)` per the spec. `R` is indexed
/// 0-based for lags 1..=K (so `R[0]` is R(1)).
pub fn phase_separation(r_k: &[f64], best_k: usize) -> f64 {
    if best_k < 2 || best_k > r_k.len() {
        return 0.0;
    }
    let same = r_k[best_k - 1];
    // Background = median R(k') for k' near best_k, k' ≠ multiple of best_k.
    let lo = best_k.saturating_sub(best_k / 2);
    let hi = (best_k + best_k / 2).min(r_k.len());
    let mut bg: Vec<f64> = (lo.max(1)..=hi)
        .filter(|&k| k != best_k && k % best_k != 0 && best_k % k != 0)
        .map(|k| r_k[k - 1])
        .collect();
    if bg.is_empty() {
        return 0.0;
    }
    bg.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = bg[bg.len() / 2];
    same - median
}

/// Primitive multiplicity correction. If `R(d) ≥ R(best_k) − δ` for
/// some divisor `d ≥ 2` of `best_k`, return the smallest such `d`.
pub fn primitive_correct(r_k: &[f64], best_k: usize, delta: f64) -> usize {
    if best_k < 2 || best_k > r_k.len() {
        return best_k;
    }
    let best_val = r_k[best_k - 1];
    let mut k = best_k;
    for d in divisors(best_k) {
        if d < 2 || d > r_k.len() {
            continue;
        }
        if r_k[d - 1] >= best_val - delta && d < k {
            k = d;
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_separation_high_for_clean_hor() {
        // R(k) peaks sharply at k=12: R(12)=0.9, everything else 0.3
        // except multiples (24 not in range).
        let mut r = vec![0.3; 16];
        r[11] = 0.9; // R(12)
        let ps = phase_separation(&r, 12);
        assert!(ps > 0.5, "expected sharp peak → high phase_sep; got {ps}");
    }

    #[test]
    fn phase_separation_low_for_uniform_curve() {
        // All R(k) ≈ 0.8 — no preferred lag → phase_separation small.
        let r = vec![0.8; 16];
        let ps = phase_separation(&r, 12);
        assert!(ps.abs() < 0.05, "uniform curve → phase_sep ≈ 0; got {ps}");
    }

    #[test]
    fn primitive_corrects_12_down_to_6() {
        // R(6)=0.88, R(12)=0.90. δ=0.05 → 6 wins.
        let mut r = vec![0.2; 16];
        r[5] = 0.88;
        r[11] = 0.90;
        assert_eq!(primitive_correct(&r, 12, 0.05), 6);
    }

    #[test]
    fn primitive_keeps_12_when_no_divisor_close() {
        let mut r = vec![0.2; 16];
        r[5] = 0.4;
        r[11] = 0.90;
        assert_eq!(primitive_correct(&r, 12, 0.05), 12);
    }

    #[test]
    fn primitive_corrects_to_smallest_qualifying_divisor() {
        // R(3) = R(4) = R(6) = R(12) - 0.01, δ=0.05 → 3 wins.
        let mut r = vec![0.2; 16];
        r[2] = 0.89;
        r[3] = 0.89;
        r[5] = 0.89;
        r[11] = 0.90;
        assert_eq!(primitive_correct(&r, 12, 0.05), 3);
    }
}

//! Candidate width expansion (`detect_impl_plan.md §6.1`, A3).
//!
//! Four-tier prioritisation rule, applied **before** the
//! `max_widths_per_array` truncation:
//!
//! 1. **Every original input period**, sorted by `period_score` desc.
//! 2. **Every valid divisor of the top-`divisor_top_n` periods**.
//! 3. **Near-misses** (±`neighborhood_n`) around tiers 1 and 2.
//! 4. **Harmonics** (`2·p`, `3·p`) and any remaining low-score extras.
//!
//! Tiers 1 + 2 are **never** dropped by the cap. If they collectively
//! exceed it, the cap is logged and the result is returned with all
//! tier-1/2 entries intact (no tier-3/4 added).

use crate::detect::config::DetectorConfig;
use crate::detect::types::PeriodCandidate;
use std::collections::BTreeSet;

pub fn expand(periods: &[PeriodCandidate], cfg: &DetectorConfig, array_len: usize) -> Vec<usize> {
    let mut sorted: Vec<&PeriodCandidate> = periods.iter().collect();
    sorted.sort_by(|a, b| {
        b.period_score
            .partial_cmp(&a.period_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut out: Vec<usize> = Vec::new();
    let mut seen: BTreeSet<usize> = BTreeSet::new();

    let in_range = |w: usize| -> bool { w >= cfg.min_width && w <= cfg.max_width && w < array_len };
    let push = |out: &mut Vec<usize>, seen: &mut BTreeSet<usize>, w: usize| {
        if in_range(w) && seen.insert(w) {
            out.push(w);
        }
    };

    // Tier 1 — every original input period.
    for c in &sorted {
        push(&mut out, &mut seen, c.period_bp);
    }

    // Tier 2 — divisors of the top-N input periods.
    let top_n = cfg.divisor_top_n.min(sorted.len());
    for c in &sorted[..top_n] {
        for d in divisors(c.period_bp) {
            push(&mut out, &mut seen, d);
        }
    }

    let tier12_len = out.len();
    if tier12_len >= cfg.max_widths_per_array {
        if tier12_len > cfg.max_widths_per_array {
            log::warn!(
                "widths: tiers 1+2 produced {} widths; exceeds cap {} — keeping all",
                tier12_len,
                cfg.max_widths_per_array
            );
        }
        return out;
    }

    // Tier 3 — neighborhoods around tiers 1+2 items.
    let tier12: Vec<usize> = out.clone();
    'tier3: for w in tier12 {
        for n in 1..=cfg.neighborhood_n {
            for sign in [-(n as i64), n as i64] {
                let new_w = w as i64 + sign;
                if new_w > 0 {
                    push(&mut out, &mut seen, new_w as usize);
                    if out.len() >= cfg.max_widths_per_array {
                        break 'tier3;
                    }
                }
            }
        }
    }

    // Tier 4 — harmonics (low-score extras).
    'tier4: for c in &sorted {
        for mul in [2usize, 3] {
            let harm = c.period_bp.saturating_mul(mul);
            push(&mut out, &mut seen, harm);
            if out.len() >= cfg.max_widths_per_array {
                break 'tier4;
            }
        }
    }

    out
}

/// Proper divisors of `p` (i.e. `1 < d < p` with `p % d == 0`).
pub fn divisors(p: usize) -> Vec<usize> {
    let mut ds = Vec::new();
    if p < 4 {
        return ds;
    }
    let mut d = 2;
    while d * d <= p {
        if p % d == 0 {
            ds.push(d);
            let other = p / d;
            if other != d && other != p {
                ds.push(other);
            }
        }
        d += 1;
    }
    ds.sort();
    ds
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;
    use crate::detect::types::PeriodCandidate;

    fn p(period_bp: usize, score: f64) -> PeriodCandidate {
        PeriodCandidate {
            array_id: "x".into(),
            period_bp,
            period_score: score,
            source: "test".into(),
        }
    }

    #[test]
    fn empty_input_empty_output() {
        let cfg = DetectorConfig::default();
        assert!(expand(&[], &cfg, 100_000).is_empty());
    }

    #[test]
    fn tier1_always_present() {
        let cfg = DetectorConfig::default();
        let ws = expand(&[p(171, 0.94), p(2052, 0.88)], &cfg, 500_000);
        assert!(ws.contains(&171));
        assert!(ws.contains(&2052));
    }

    #[test]
    fn divisors_of_top_period_included() {
        let cfg = DetectorConfig::default();
        let ws = expand(&[p(2052, 0.88), p(171, 0.94)], &cfg, 500_000);
        // 2052 = 2*2*3*3*3*19; proper divisors >= min_width=20 include
        // 27, 36, 38, 54, 76, 108, 114, 171, 228, 342, 513, 684, 1026.
        // 171's divisors 9 and 19 fall below min_width=20 and are dropped.
        for d in [57usize, 171, 2052] {
            assert!(ws.contains(&d), "missing divisor {d}: ws={ws:?}");
        }
        assert!(!ws.contains(&9));
        assert!(!ws.contains(&19));
    }

    #[test]
    fn divisors_below_min_width_kept_when_min_width_lower() {
        let mut cfg = DetectorConfig::default();
        cfg.min_width = 5;
        let ws = expand(&[p(171, 0.94)], &cfg, 500_000);
        assert!(ws.contains(&9));
        assert!(ws.contains(&19));
        assert!(ws.contains(&57));
    }

    #[test]
    fn out_of_range_filtered() {
        let mut cfg = DetectorConfig::default();
        cfg.min_width = 100;
        cfg.max_width = 1000;
        let ws = expand(&[p(171, 0.94), p(2052, 0.88)], &cfg, 10_000);
        // 2052 > max_width → dropped.
        assert!(!ws.contains(&2052));
        // 9, 19 (divisors of 171) < min_width → dropped.
        assert!(!ws.contains(&9));
        assert!(!ws.contains(&19));
        // 171 in range.
        assert!(ws.contains(&171));
    }

    #[test]
    fn period_above_array_len_dropped() {
        let cfg = DetectorConfig::default();
        let ws = expand(&[p(171, 0.94), p(10_000, 0.50)], &cfg, 5_000);
        assert!(!ws.contains(&10_000));
        assert!(ws.contains(&171));
    }

    #[test]
    fn cap_keeps_all_of_tier12_even_if_over_cap() {
        let mut cfg = DetectorConfig::default();
        cfg.max_widths_per_array = 5;
        cfg.divisor_top_n = 1;
        // 720 has many divisors above min_width=20.
        let ws = expand(&[p(720, 0.94)], &cfg, 100_000);
        assert!(ws.contains(&720));
    }

    #[test]
    fn neighborhoods_added_around_tier12() {
        let cfg = DetectorConfig::default();
        let ws = expand(&[p(171, 0.94)], &cfg, 500_000);
        for off in [-3i64, -2, -1, 1, 2, 3] {
            let w = (171i64 + off) as usize;
            assert!(ws.contains(&w), "missing neighborhood {w}: ws={ws:?}");
        }
    }

    #[test]
    fn harmonics_added() {
        let cfg = DetectorConfig::default();
        let ws = expand(&[p(171, 0.94)], &cfg, 500_000);
        assert!(ws.contains(&342));
        assert!(ws.contains(&513));
    }

    #[test]
    fn divisors_helper_basics() {
        assert_eq!(divisors(12), vec![2, 3, 4, 6]);
        assert_eq!(divisors(171), vec![3, 9, 19, 57]);
        assert!(divisors(1).is_empty());
        assert!(divisors(2).is_empty());
        assert!(divisors(7).is_empty());
    }
}

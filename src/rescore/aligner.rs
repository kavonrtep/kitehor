//! Banded semi-global edit distance kernel for `rescore`.
//!
//! The pattern `a` is matched end-to-end; the text `b` has free start and
//! free end. The returned value is the minimum number of substitutions +
//! insertions + deletions over all alignments of `a` to some substring of
//! `b`, capped by the band geometry.
//!
//! `band` is an **indel-deviation tolerance**, not an edit-distance cap.
//! Mismatches accumulate along the main diagonal at zero band cost, so a
//! divergent monomer pair (e.g. 25 % mismatches but no net indels) is
//! reported faithfully even with `band = 20`. Only alignments whose net
//! indel count exceeds the band are excluded — for our use case (HOR
//! variants with single-digit indel rates per tile) the natural default
//! `max(20, 2·slop)` is comfortably generous.
//!
//! N handling: `b'N'` matches nothing, including another `b'N'`. The
//! sampler in `src/rescore/sample.rs` filters out windows above a small Ns
//! threshold, so the kernel rarely sees N-heavy input in practice.
//!
//! Cells use `u16`; the sentinel is `u16::MAX`. Real edit-distance values
//! are bounded by `m + n` which is well under `u16::MAX = 65535` for any
//! realistic period. `saturating_add` keeps the sentinel sticky.

/// Reusable per-call scratch storage. Hand one to each rayon worker via
/// `par_iter().map_init(Scratch::new, ...)` to avoid re-allocating the
/// two DP rows for every (record, period, sample) call.
#[derive(Debug, Default, Clone)]
pub struct Scratch {
    prev: Vec<u16>,
    curr: Vec<u16>,
}

impl Scratch {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Sentinel returned in `BandedResult::j_end` when no in-band alignment
/// reached row m (e.g. band too tight for the geometry). Callers treat
/// this as "no shift information available".
pub const J_END_NONE: u16 = u16::MAX;

/// Result of one banded semi-global alignment call.
///
/// `j_end` is the column in `b` where the optimal path exits at row `m`.
/// For our caller, the *start* column of the alignment in `b` is
/// approximated as `j_end − m`; this is exact when the alignment has no
/// net indels and off by at most `band` otherwise. Compared to the
/// caller's "natural" mapping (`j_start = slop` because `b` is the
/// adjacent tile extended by ±slop), the shift is
/// `(j_end − m) − slop = j_end − m − slop` bp — positive when the
/// alignment landed downstream of the claimed period, negative upstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BandedResult {
    pub distance: u32,
    pub j_end: u16,
}

/// Per-cell costs used by the DP recurrence. Match cost is always 0
/// (edit-distance convention). Insertions and deletions share a single
/// gap cost — affine gaps are not modeled.
///
/// When `mismatch_cost != 1` or `gap_cost != 1`, the returned value is a
/// *weighted* edit distance, not an operation count. Callers should
/// remember that `identity_from_distance` derives identity from the
/// returned value over `|A|`, so non-unit costs produce a weighted
/// identity (still in `[0, 1]`, monotone in alignment quality, but not
/// equal to "matching positions / pattern length").
#[derive(Debug, Clone, Copy)]
pub struct ScoringConfig {
    pub mismatch_cost: u16,
    pub gap_cost: u16,
}

impl Default for ScoringConfig {
    fn default() -> Self {
        Self {
            mismatch_cost: 1,
            gap_cost: 1,
        }
    }
}

/// Banded semi-global edit distance.
///
/// `a` is the pattern (length `m`); `b` is the text (length `n`). Both
/// `a`'s start and end are constrained to lie on the boundary of the
/// alignment; `b`'s start and end are free within the band envelope.
///
/// The DP envelope is:
///
/// ```text
///     j_lo(i) = max(0, i - band)
///     j_hi(i) = min(n, i + (n - m) + band)
/// ```
///
/// Cells outside the envelope are treated as +∞. If the envelope is empty
/// at any row (band too tight relative to the geometry), the function
/// returns an inf sentinel ≥ `m + n + band + 1`, which the caller treats
/// as "no in-band alignment" (identity_from_distance clamps to 0).
///
/// Cost: O(m · (2·band + (n - m) + 1)) cells per call.
pub fn semiglobal_edit_distance_banded(
    a: &[u8],
    b: &[u8],
    band: usize,
    scoring: &ScoringConfig,
    scratch: &mut Scratch,
) -> BandedResult {
    let m = a.len();
    let n = b.len();
    if m == 0 {
        return BandedResult {
            distance: 0,
            j_end: 0,
        };
    }
    if n == 0 {
        // Must delete all of A; the path exits at column 0.
        return BandedResult {
            distance: (m as u32).saturating_mul(scoring.gap_cost as u32),
            j_end: 0,
        };
    }

    let offset_max = n.saturating_sub(m);
    let inf: u16 = u16::MAX;

    // Ensure capacity and reset both rows to the sentinel. Resizing is
    // a no-op once the worker has seen a pair of this length.
    if scratch.prev.len() < n + 1 {
        scratch.prev.resize(n + 1, inf);
    } else {
        scratch.prev[..=n].fill(inf);
    }
    if scratch.curr.len() < n + 1 {
        scratch.curr.resize(n + 1, inf);
    } else {
        scratch.curr[..=n].fill(inf);
    }

    // Row 0: free start. D[0][j] = 0 for j ∈ [0, offset_max]; +∞ elsewhere.
    for j in 0..=offset_max.min(n) {
        scratch.prev[j] = 0;
    }

    for i in 1..=m {
        let j_lo = i.saturating_sub(band);
        let j_hi = (i + offset_max + band).min(n);

        if j_lo > j_hi {
            // Empty band on this row ⇒ reset the soon-to-be-prev row to
            // inf so reads in iteration i+1 propagate correctly.
            scratch.curr[..=n].fill(inf);
            std::mem::swap(&mut scratch.prev, &mut scratch.curr);
            continue;
        }

        // Sentinel the two cells immediately outside the band — they may
        // be read as curr[j_lo(i+1)-1] or (after swap) prev[j_hi(i)+1]
        // in the next iteration. Cells further outside are never read.
        if j_lo > 0 {
            scratch.curr[j_lo - 1] = inf;
        }
        if j_hi < n {
            scratch.curr[j_hi + 1] = inf;
        }

        // Boundary D[i][0] only valid in band (i.e. j_lo == 0). Cost
        // is i deletions from A, weighted by `gap_cost`.
        if j_lo == 0 {
            scratch.curr[0] = (i as u16).saturating_mul(scoring.gap_cost);
        }

        let ai = a[i - 1];
        let j_start = j_lo.max(1);
        for j in j_start..=j_hi {
            let bj = b[j - 1];
            // N matches nothing (including N).
            let cell_cost: u16 = if ai == bj && ai != b'N' {
                0
            } else {
                scoring.mismatch_cost
            };
            let diag = scratch.prev[j - 1].saturating_add(cell_cost);
            let up = scratch.prev[j].saturating_add(scoring.gap_cost);
            let left = scratch.curr[j - 1].saturating_add(scoring.gap_cost);
            scratch.curr[j] = diag.min(up).min(left);
        }
        std::mem::swap(&mut scratch.prev, &mut scratch.curr);
    }

    // Free end: argmin over D[m][j] for j ∈ [j_lo(m), j_hi(m)] ∩ [0, n].
    // Result is in `scratch.prev` after the final swap.
    let j_lo = m.saturating_sub(band);
    let j_hi = (m + offset_max + band).min(n);
    if j_lo > j_hi {
        return BandedResult {
            distance: inf as u32,
            j_end: J_END_NONE,
        };
    }
    let window = &scratch.prev[j_lo..=j_hi];
    // Argmin: first occurrence on ties (cheapest path geometrically).
    let (best_idx, best_val) = window
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| **v)
        .map(|(i, v)| (i, *v))
        .unwrap_or((0, inf));
    BandedResult {
        distance: best_val as u32,
        j_end: (j_lo + best_idx) as u16,
    }
}

/// Identity ∈ [0, 1] from edit distance and pattern length. Clamped at 0.
#[inline]
pub fn identity_from_distance(edit_distance: u32, pattern_len: usize) -> f64 {
    if pattern_len == 0 {
        return 1.0;
    }
    let id = 1.0 - (edit_distance as f64) / (pattern_len as f64);
    if id < 0.0 {
        0.0
    } else {
        id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Convenience wrapper used by existing functional tests: a band
    /// large enough that the kernel behaves identically to unbanded DP,
    /// and default unit costs. Returns just the edit distance for
    /// existing assertions.
    fn dist(a: &[u8], b: &[u8]) -> u32 {
        let band = (a.len() + b.len()).max(64);
        let mut s = Scratch::new();
        let cfg = ScoringConfig::default();
        semiglobal_edit_distance_banded(a, b, band, &cfg, &mut s).distance
    }

    #[test]
    fn exact_match_is_zero() {
        let a = b"ACGTACGT";
        assert_eq!(dist(a, a), 0);
    }

    #[test]
    fn one_mismatch() {
        let a = b"ACGTACGT";
        let b = b"ACGTAGGT";
        assert_eq!(dist(a, b), 1);
    }

    #[test]
    fn one_insertion_in_b() {
        let a = b"ACGTACGT";
        let b = b"ACGTNACGT";
        let d = dist(a, b);
        assert!(d <= 1, "expected ≤1, got {}", d);
    }

    #[test]
    fn one_deletion_in_a_vs_b() {
        let a = b"ACGTCGT";
        let b = b"ACGTACGT";
        assert_eq!(dist(a, b), 1);
    }

    #[test]
    fn substring_of_b_is_free_match() {
        let a = b"CAT";
        let b = b"GGCATTT";
        assert_eq!(dist(a, b), 0);
    }

    #[test]
    fn a_longer_than_b() {
        let a = b"ACGTACGT";
        let b = b"ACGT";
        assert_eq!(dist(a, b), 4);
    }

    #[test]
    fn empty_a() {
        assert_eq!(dist(b"", b"ACGT"), 0);
    }

    #[test]
    fn empty_b() {
        assert_eq!(dist(b"ACG", b""), 3);
    }

    #[test]
    fn n_matches_nothing() {
        assert_eq!(dist(b"NCG", b"ACG"), 1);
        assert_eq!(dist(b"ACG", b"ANG"), 1);
        assert_eq!(dist(b"NCG", b"NCG"), 1);
    }

    #[test]
    fn identity_basic() {
        assert!((identity_from_distance(0, 100) - 1.0).abs() < 1e-12);
        assert!((identity_from_distance(10, 100) - 0.9).abs() < 1e-12);
        assert!((identity_from_distance(100, 100) - 0.0).abs() < 1e-12);
        assert_eq!(identity_from_distance(200, 100), 0.0);
        assert_eq!(identity_from_distance(0, 0), 1.0);
    }

    /// Brute-force oracle. Matches the kernel's free-start semantics:
    /// `D[0][j] = 0` for `j ∈ [0, n - m]` (the natural mapping window),
    /// `+∞` elsewhere. This is *not* unrestricted semi-global — it's
    /// what `semiglobal_edit_distance_banded(_, _, ∞)` computes.
    fn dp_reference(a: &[u8], b: &[u8]) -> u32 {
        let m = a.len();
        let n = b.len();
        let offset_max = n.saturating_sub(m);
        let inf: u32 = (m + n + 1) as u32;
        let mut d = vec![vec![inf; n + 1]; m + 1];
        let zero_to = offset_max.min(n);
        d[0][..=zero_to].fill(0);
        for (i, row) in d.iter_mut().enumerate().take(m + 1).skip(1) {
            row[0] = i as u32;
        }
        for i in 1..=m {
            for j in 1..=n {
                let ai = a[i - 1];
                let bj = b[j - 1];
                let cost = if ai == bj && ai != b'N' { 0 } else { 1 };
                d[i][j] = d[i - 1][j - 1]
                    .saturating_add(cost)
                    .min(d[i - 1][j].saturating_add(1))
                    .min(d[i][j - 1].saturating_add(1));
            }
        }
        *d[m].iter().min().unwrap()
    }

    #[test]
    fn matches_reference_when_band_is_generous() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(1);
        let bases = b"ACGT";
        let mut s = Scratch::new();
        for _ in 0..200 {
            let m = rng.random_range(1..50);
            let n = rng.random_range(1..50);
            let a: Vec<u8> = (0..m).map(|_| bases[rng.random_range(0..4)]).collect();
            let b: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
            let band = m + n + 5;
            let banded =
                semiglobal_edit_distance_banded(&a, &b, band, &ScoringConfig::default(), &mut s);
            let oracle = dp_reference(&a, &b);
            assert_eq!(
                banded.distance,
                oracle,
                "mismatch a={:?} b={:?}",
                std::str::from_utf8(&a),
                std::str::from_utf8(&b)
            );
        }
    }

    #[test]
    fn band_never_underreports() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(2);
        let bases = b"ACGT";
        let mut s = Scratch::new();
        for _ in 0..200 {
            let m = rng.random_range(5..40);
            let n = rng.random_range(m..m + 10);
            let a: Vec<u8> = (0..m).map(|_| bases[rng.random_range(0..4)]).collect();
            let b: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
            let truth = dp_reference(&a, &b);
            for &band in &[0usize, 1, 2, 5, 10, m + n] {
                let banded = semiglobal_edit_distance_banded(
                    &a,
                    &b,
                    band,
                    &ScoringConfig::default(),
                    &mut s,
                );
                assert!(
                    banded.distance >= truth,
                    "band={} reported {} but true is {} (a={:?}, b={:?})",
                    band,
                    banded.distance,
                    truth,
                    std::str::from_utf8(&a),
                    std::str::from_utf8(&b)
                );
            }
        }
    }

    #[test]
    fn matches_reference_when_truth_fits_in_band() {
        let a = b"ACGTACGTACGT";
        let b = b"ACGTAGCTACAT";
        let truth = dp_reference(a, b);
        assert!(truth <= 3, "unexpected truth={}", truth);
        let mut s = Scratch::new();
        for &band in &[0usize, 1, 2, 5, 100] {
            let banded =
                semiglobal_edit_distance_banded(a, b, band, &ScoringConfig::default(), &mut s);
            assert_eq!(
                banded.distance, truth,
                "band={} truth={} got={}",
                band, truth, banded.distance
            );
        }
    }

    #[test]
    fn band_zero_caps_at_indel_count() {
        let a = b"ACGTACGT";
        let b = b"ACGTNACGT";
        let truth = dp_reference(a, b);
        assert!(truth <= 1);
        let mut s = Scratch::new();
        let banded_zero =
            semiglobal_edit_distance_banded(a, b, 0, &ScoringConfig::default(), &mut s);
        assert!(banded_zero.distance >= truth);
    }

    #[test]
    fn empty_band_with_unreachable_geometry_returns_inf_sentinel() {
        let a = b"ACGTACGT";
        let b = b"ACGT";
        let mut s = Scratch::new();
        let r = semiglobal_edit_distance_banded(a, b, 2, &ScoringConfig::default(), &mut s);
        assert!(
            r.distance > a.len() as u32,
            "expected sentinel, got {}",
            r.distance
        );
        assert_eq!(r.j_end, J_END_NONE);
    }

    #[test]
    fn scratch_reuse_matches_fresh_scratch() {
        // Property: reusing one Scratch across many heterogeneous calls
        // produces the same answers as a fresh Scratch each call.
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(3);
        let bases = b"ACGT";
        let mut reused = Scratch::new();
        for _ in 0..100 {
            let m = rng.random_range(1..60);
            let n = rng.random_range(m..m + 20);
            let a: Vec<u8> = (0..m).map(|_| bases[rng.random_range(0..4)]).collect();
            let b: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
            let band = rng.random_range(0..30);
            let reused_result = semiglobal_edit_distance_banded(
                &a,
                &b,
                band,
                &ScoringConfig::default(),
                &mut reused,
            );
            let fresh_result = semiglobal_edit_distance_banded(
                &a,
                &b,
                band,
                &ScoringConfig::default(),
                &mut Scratch::new(),
            );
            assert_eq!(
                reused_result,
                fresh_result,
                "reuse divergence band={} a={:?} b={:?}",
                band,
                std::str::from_utf8(&a),
                std::str::from_utf8(&b)
            );
        }
    }

    #[test]
    fn parallel_workers_get_same_results_as_sequential() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        use rayon::prelude::*;
        let mut rng = StdRng::seed_from_u64(4);
        let bases = b"ACGT";
        let cases: Vec<(Vec<u8>, Vec<u8>, usize)> = (0..200)
            .map(|_| {
                let m = rng.random_range(5..80);
                let n = rng.random_range(m..m + 20);
                let a: Vec<u8> = (0..m).map(|_| bases[rng.random_range(0..4)]).collect();
                let b: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
                let band = rng.random_range(0..40);
                (a, b, band)
            })
            .collect();

        // Sequential: fresh Scratch each call.
        let seq: Vec<BandedResult> = cases
            .iter()
            .map(|(a, b, band)| {
                let mut s = Scratch::new();
                semiglobal_edit_distance_banded(a, b, *band, &ScoringConfig::default(), &mut s)
            })
            .collect();

        // Parallel: one Scratch per worker via map_init.
        let par: Vec<BandedResult> = cases
            .par_iter()
            .map_init(Scratch::new, |scratch, (a, b, band)| {
                semiglobal_edit_distance_banded(a, b, *band, &ScoringConfig::default(), scratch)
            })
            .collect();

        assert_eq!(seq, par, "parallel vs sequential divergence");
    }

    #[test]
    fn mismatch_cost_scales_weighted_distance() {
        // One mismatch in a 4 bp pair, generous band so geometry is fine.
        let a = b"ACGT";
        let b = b"AGGT"; // pos 1: C→G
        let band = 16;
        let mut s = Scratch::new();
        let d1 = semiglobal_edit_distance_banded(
            a,
            b,
            band,
            &ScoringConfig {
                mismatch_cost: 1,
                gap_cost: 1,
            },
            &mut s,
        );
        let d3 = semiglobal_edit_distance_banded(
            a,
            b,
            band,
            &ScoringConfig {
                mismatch_cost: 3,
                gap_cost: 1,
            },
            &mut s,
        );
        assert_eq!(d1.distance, 1);
        // Higher mismatch cost: either pay 3 for the mismatch, or 2 gaps
        // (one in A, one in B) for cost 2. The DP picks the cheaper path.
        assert_eq!(d3.distance, 2);
    }

    #[test]
    fn gap_cost_scales_weighted_distance() {
        // Forced one-gap alignment.
        let a = b"ACGT"; // 4
        let b = b"ACGGT"; // 5 — one insertion in B
        let band = 16;
        let mut s = Scratch::new();
        let d1 = semiglobal_edit_distance_banded(
            a,
            b,
            band,
            &ScoringConfig {
                mismatch_cost: 1,
                gap_cost: 1,
            },
            &mut s,
        );
        let d4 = semiglobal_edit_distance_banded(
            a,
            b,
            band,
            &ScoringConfig {
                mismatch_cost: 1,
                gap_cost: 4,
            },
            &mut s,
        );
        assert_eq!(d1.distance, 1);
        // With gap_cost=4 and free B-start (offset_max=1), the cheapest
        // alignment is to shift B-start by 1, giving 1 mismatch (cost 1).
        assert_eq!(d4.distance, 1);
    }

    #[test]
    fn high_gap_cost_forces_mismatch_path() {
        // Even-length pair with one true mismatch; with high gap cost
        // the optimum can't escape via gaps.
        let a = b"ACGTACGT";
        let b = b"ACGTNCGT"; // pos 4: A→N (treated as mismatch)
        let band = 16;
        let mut s = Scratch::new();
        let d = semiglobal_edit_distance_banded(
            a,
            b,
            band,
            &ScoringConfig {
                mismatch_cost: 1,
                gap_cost: 10,
            },
            &mut s,
        );
        // Optimum is the single 1-cost mismatch on the diagonal; no gap
        // is cheaper. Result must be 1, independent of gap_cost.
        assert_eq!(d.distance, 1);
    }

    // --- j_end / shift coverage --------------------------------------------

    #[test]
    fn j_end_on_exact_match_equals_m_plus_offset_max() {
        // a == b ⇒ optimal path is the diagonal from (0,0) to (m,m).
        // n == m so offset_max = 0; j_end must be m.
        let a = b"ACGTACGTAC";
        let m = a.len();
        let mut s = Scratch::new();
        let r = semiglobal_edit_distance_banded(a, a, 16, &ScoringConfig::default(), &mut s);
        assert_eq!(r.distance, 0);
        assert_eq!(r.j_end as usize, m);
    }

    #[test]
    fn j_end_recovers_shift_when_b_is_shifted() {
        // Construct a "tile-pair" geometry: a is a clean tile, b is the
        // adjacent tile with ±slop slack on each side. Shift the "true
        // tile" inside b by δ bp from the natural mapping.
        //
        // Natural mapping: a starts at column `slop` in b, so j_end = m + slop.
        // Shifted by δ:    a starts at column `slop + δ`, so j_end = m + slop + δ.
        let m = 50usize;
        let slop = 10usize;
        let band = 20usize;

        // Random-but-deterministic 50 bp "tile" content.
        let tile: Vec<u8> = (0..m).map(|i| b"ACGT"[i % 4]).collect();

        for &delta in &[-6i32, -3, 0, 2, 5] {
            // Build b of length m + 2·slop, with the tile placed at
            // position `slop + delta` inside b. Fill the rest with a
            // distinct base so the kernel can only match the inserted
            // tile.
            let mut b = vec![b'N'; m + 2 * slop];
            let start = (slop as i32 + delta) as usize;
            b[start..start + m].copy_from_slice(&tile);

            let mut s = Scratch::new();
            let r =
                semiglobal_edit_distance_banded(&tile, &b, band, &ScoringConfig::default(), &mut s);
            // The kernel should find the tile inside b (distance 0) at
            // exactly column m + slop + delta.
            assert_eq!(
                r.distance, 0,
                "delta={} expected distance 0, got {}",
                delta, r.distance
            );
            let expected_j_end = (m as i32 + slop as i32 + delta) as u16;
            assert_eq!(
                r.j_end, expected_j_end,
                "delta={} expected j_end={}, got {}",
                delta, expected_j_end, r.j_end
            );
        }
    }

    #[test]
    fn j_end_is_inside_band_window() {
        // Property: for any random pair, j_end is inside [j_lo(m), j_hi(m)]
        // (or J_END_NONE if the band is empty).
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(5);
        let bases = b"ACGT";
        let mut s = Scratch::new();
        for _ in 0..100 {
            let m: usize = rng.random_range(5..30);
            let n: usize = rng.random_range(m..m + 15);
            let a: Vec<u8> = (0..m).map(|_| bases[rng.random_range(0..4)]).collect();
            let b: Vec<u8> = (0..n).map(|_| bases[rng.random_range(0..4)]).collect();
            let band: usize = rng.random_range(0..15);
            let r =
                semiglobal_edit_distance_banded(&a, &b, band, &ScoringConfig::default(), &mut s);
            if r.j_end == J_END_NONE {
                continue;
            }
            let offset_max = n.saturating_sub(m);
            let j_lo = m.saturating_sub(band);
            let j_hi = (m + offset_max + band).min(n);
            assert!(
                (r.j_end as usize) >= j_lo && (r.j_end as usize) <= j_hi,
                "j_end={} outside [{}, {}] (m={}, n={}, band={})",
                r.j_end,
                j_lo,
                j_hi,
                m,
                n,
                band
            );
        }
    }
}

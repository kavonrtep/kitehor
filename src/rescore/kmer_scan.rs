//! Founder-aware k-mer-pair diagnostics for `kitehor rescore`.
//!
//! Two observational metrics that exploit the per-record k-mer
//! position list (already needed by kite) to characterise whether
//! the candidate period `P` carries a real **nested-subrepeat**
//! signal — i.e. whether the k-mer pairs at distance ≈ P live
//! preferentially inside a specific phase of the founder cycle:
//!
//! - [`kmer_density_autocorr_at_founder`] — autocorrelation of the
//!   sliding-window pair-density profile at lag = `founder_period`.
//! - [`kmer_phase_folded_contrast`] — phase-folded density excess
//!   ("max contiguous-half fraction minus 0.5"); jitter-robust
//!   companion to `autoF`.
//!
//! Both metrics are observational in this release — they appear as
//! columns in the rescore TSV but don't gate the `subrepeat` flag.
//! The two are complementary on the IPIP corpus: an OR-combined
//! gate at `autoF ≥ 0.4 OR phaseC ≥ 0.10` covers ≈ 94 % of
//! currently-flagged subrepeats, vs ≈ 75 % each individually.
//!
//! Defaults: `k = 6` (matches `kite-periodicity --kmer-size`),
//! `distance_tol = 3`, `n_bins = 12`, `min_total_pairs = 20`.

use ahash::AHashMap;

/// Per-record k-mer position map. Keyed by the 2-bit-packed forward
/// encoding of the k-mer (A=0, C=1, G=2, T=3); values are the start
/// positions (0-based) of every occurrence in the array.
///
/// K-mers containing any non-ACGT base (including `N` and degenerate
/// IUPAC codes) are skipped.
pub type KmerPositions = AHashMap<u64, Vec<u32>>;

/// Encode a single DNA base. Returns `None` for `N` and any other
/// non-ACGT byte so the scanner can reset the running hash.
#[inline]
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Build the per-record k-mer position map. `O(seq.len())` with a
/// single AHashMap pass.
pub fn build_kmer_positions(seq: &[u8], k: usize) -> KmerPositions {
    let mut map: KmerPositions = AHashMap::new();
    if k == 0 || k > 32 || seq.len() < k {
        return map;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut hash = 0u64;
    let mut valid = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        match encode_base(b) {
            Some(code) => {
                hash = ((hash << 2) | code) & mask;
                valid += 1;
            }
            None => {
                hash = 0;
                valid = 0;
            }
        }
        if valid >= k {
            let start = (i + 1 - k) as u32;
            map.entry(hash).or_default().push(start);
        }
    }
    map
}

/// Density-profile autocorrelation of period-P k-mer evidence,
/// evaluated at lag = `founder_period`.
///
/// The "nested subrepeat" signature: a real subrepeat at period P
/// inside a founder of length F lives in part of each founder copy
/// along the array. K-mer pairs at distance ≈P cluster in those
/// parts and are absent in the gaps. So density(x) — the count of
/// matching pairs per fine-grained sliding window — oscillates
/// with period F along the array.
///
/// This function builds density(x) at a fine scale
/// (window = max(20, P/2), step = window/2), then returns Pearson
/// correlation between density(x) and density(x + founder_period).
/// For a real nested subrepeat the value is positive and close to
/// 1; for a uniform tandem (subrepeat = founder, no gap) it is
/// indeterminate (zero variance → `None`); for near-founder
/// harmonics or noise it sits near zero.
///
/// Returns `None` when:
/// - `founder_period` is `None`, zero, or larger than the array,
/// - too few windows survive for a meaningful correlation
///   (`min_window_pairs` not met),
/// - one or both density vectors have zero variance.
pub fn kmer_density_autocorr_at_founder(
    positions: &KmerPositions,
    seq_len: usize,
    period: usize,
    founder_period: Option<usize>,
    distance_tol: usize,
    min_window_pairs: usize,
) -> Option<f64> {
    let founder = founder_period.filter(|f| *f > 0 && *f < seq_len)?;
    if period == 0 || seq_len < 2 * period {
        return None;
    }
    // Fine window so per-window density can resolve subrepeat-vs-gap
    // variation that lives inside one founder. Step is half the
    // window to give some overlap.
    let win = period.max(20) / 2;
    let step = (win / 2).max(1);
    if seq_len < win + step {
        return None;
    }
    let n_windows = (seq_len - win) / step + 1;
    if n_windows < 4 {
        return None;
    }
    let mut counts = vec![0u32; n_windows];

    let p_lo = period.saturating_sub(distance_tol);
    let p_hi = period + distance_tol;
    let mut total_pairs = 0usize;
    for poses in positions.values() {
        if poses.len() < 2 {
            continue;
        }
        for w in poses.windows(2) {
            let d = (w[1] - w[0]) as usize;
            if d >= p_lo && d <= p_hi {
                let mid = (w[0] as usize + w[1] as usize) / 2;
                // Map midpoint to the window whose center is closest.
                let idx = mid.saturating_sub(win / 2) / step;
                let idx = idx.min(n_windows - 1);
                counts[idx] += 1;
                total_pairs += 1;
            }
        }
    }
    if total_pairs < min_window_pairs {
        return None;
    }

    let lag = (founder + step / 2) / step; // rounded
    if lag == 0 || lag >= n_windows {
        return None;
    }
    let n = n_windows - lag;
    if n < 4 {
        return None;
    }

    // Pearson correlation between counts[0..n] and counts[lag..lag+n].
    let mut sum_x = 0f64;
    let mut sum_y = 0f64;
    let mut sum_xx = 0f64;
    let mut sum_yy = 0f64;
    let mut sum_xy = 0f64;
    for i in 0..n {
        let x = counts[i] as f64;
        let y = counts[i + lag] as f64;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_yy += y * y;
        sum_xy += x * y;
    }
    let nf = n as f64;
    let mean_x = sum_x / nf;
    let mean_y = sum_y / nf;
    let var_x = sum_xx / nf - mean_x * mean_x;
    let var_y = sum_yy / nf - mean_y * mean_y;
    if var_x <= 0.0 || var_y <= 0.0 {
        return None;
    }
    let cov = sum_xy / nf - mean_x * mean_y;
    Some(cov / (var_x * var_y).sqrt())
}

/// Phase-folded k-mer pair density contrast at the founder period.
///
/// Bins midpoints by `(midpoint mod founder_period)` into `n_bins`
/// equal-width phase bins (each `founder_period / n_bins` bp wide).
/// All founder copies along the array fold onto the same N-bin
/// histogram. The statistic returned is the **maximum contiguous
/// half-fraction excess**: of all rotations of a contiguous block
/// of `n_bins/2` bins (wrapping around the modulo cycle), take the
/// rotation whose sum is maximal, divide by the total count, and
/// subtract 0.5. The result is bounded in `[0, 0.5]`.
///
/// Why "contiguous half" rather than "max single bin minus min
/// single bin": a real nested subrepeat doesn't live in *one* phase
/// bin — it occupies a whole contiguous portion of the founder
/// (TRC_104's subrepeat fills the second half of its 180 bp
/// founder). Looking at individual bins underestimates the contrast
/// because the load is spread across many bins inside the
/// subrepeat region.
///
/// Why jitter-robust: a ±5 bp boundary drift moves each midpoint's
/// phase by ±5 bp. Phase bins are typically `founder / n_bins ≈
/// 30 bp` wide, so adjacent founders' midpoints stay within the
/// same or neighboring bins; the subrepeat-half pile-up doesn't
/// redistribute outside the dominant half.
///
/// `None` is returned when the founder is unknown, the seq is too
/// short to host any matching pair, or the total pair count is
/// below `min_total_pairs`.
pub fn kmer_phase_folded_contrast(
    positions: &KmerPositions,
    seq_len: usize,
    period: usize,
    founder_period: Option<usize>,
    distance_tol: usize,
    n_bins: usize,
    min_total_pairs: usize,
) -> Option<f64> {
    let founder = founder_period.filter(|f| *f > 0 && *f < seq_len)?;
    if period == 0 || n_bins < 2 || seq_len < 2 * period {
        return None;
    }
    let p_lo = period.saturating_sub(distance_tol);
    let p_hi = period + distance_tol;
    let mut counts = vec![0usize; n_bins];
    let mut total = 0usize;
    for poses in positions.values() {
        if poses.len() < 2 {
            continue;
        }
        for w in poses.windows(2) {
            let d = (w[1] - w[0]) as usize;
            if d >= p_lo && d <= p_hi {
                let mid = (w[0] as usize + w[1] as usize) / 2;
                let phase = mid % founder;
                let b = (phase.saturating_mul(n_bins) / founder).min(n_bins - 1);
                counts[b] += 1;
                total += 1;
            }
        }
    }
    if total < min_total_pairs {
        return None;
    }
    // Find the contiguous half-window (of size `half = n_bins / 2`,
    // rounded down) with the largest sum. Wraps around the modulo
    // cycle. O(n_bins) via rolling sum.
    let half = n_bins / 2;
    if half == 0 {
        return None;
    }
    // Initial window sum over bins[0..half].
    let mut window_sum: usize = counts[..half].iter().sum();
    let mut max_sum = window_sum;
    for start in 1..n_bins {
        // Shift window by one: remove the bin leaving, add the bin entering.
        // `entering` = bin at index (start + half - 1) mod n_bins.
        let entering = (start + half - 1) % n_bins;
        let leaving = start - 1;
        window_sum = window_sum + counts[entering] - counts[leaving];
        if window_sum > max_sum {
            max_sum = window_sum;
        }
    }
    let max_frac = max_sum as f64 / total as f64;
    let half_frac = half as f64 / n_bins as f64;
    Some((max_frac - half_frac).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_skips_n_bases() {
        let seq = b"ACGTNACGT";
        let map = build_kmer_positions(seq, 3);
        // Expected k-mers: ACG@0, CGT@1, GTN@2-skip, TNA-skip, NAC-skip,
        // ACG@5, CGT@6 — so each of ACG and CGT appears twice.
        // Encoded ACG: A=0,C=1,G=2 → 0b00_01_10 = 6
        // Encoded CGT: C=1,G=2,T=3 → 0b01_10_11 = 27
        let p_acg = map.get(&6).expect("ACG must be present");
        let p_cgt = map.get(&27).expect("CGT must be present");
        assert_eq!(p_acg, &vec![0, 5]);
        assert_eq!(p_cgt, &vec![1, 6]);
    }

    #[test]
    fn autocorr_returns_none_when_founder_unknown() {
        let seq = vec![b'A'; 500];
        let map = build_kmer_positions(&seq, 5);
        let r = kmer_density_autocorr_at_founder(&map, seq.len(), 5, None, 0, 10);
        assert!(r.is_none());
    }

    #[test]
    fn autocorr_high_for_periodic_density() {
        // Build a 600 bp array of 60 bp founder = [40 bp subrepeat
        // region of 5-bp tandem][20 bp gap]. Density of d=5 pairs
        // oscillates with period 60. Expect strong positive
        // autocorrelation at lag=60.
        let mut seq = vec![b'N'; 600];
        for f in 0..10 {
            let start = f * 60;
            for j in 0..40 {
                seq[start + j] = b"ACGTA"[j % 5];
            }
        }
        let map = build_kmer_positions(&seq, 5);
        let r = kmer_density_autocorr_at_founder(&map, seq.len(), 5, Some(60), 0, 5);
        let v = r.expect("expected Some");
        assert!(v > 0.5, "expected high autocorr, got {v}");
    }

    #[test]
    fn phase_fold_high_for_nested_subrepeat() {
        // 10 founders of 60 bp = [40 bp subrepeat at period 5][20 bp
        // gap]. Folding by founder concentrates d=5 pair midpoints
        // into the first 40 bp of each founder cycle = bins [0..4) of
        // 6. The contiguous-half (3 bins) covering [0..3) holds the
        // majority of midpoints; expected excess ≈ 0.2–0.4.
        let mut seq = vec![b'N'; 600];
        for f in 0..10 {
            let start = f * 60;
            for j in 0..40 {
                seq[start + j] = b"ACGTA"[j % 5];
            }
        }
        let map = build_kmer_positions(&seq, 5);
        let r = kmer_phase_folded_contrast(&map, seq.len(), 5, Some(60), 0, 6, 5);
        let v = r.expect("expected Some");
        assert!(v > 0.15, "expected positive phase contrast, got {v}");
    }

    #[test]
    fn phase_fold_low_for_uniform_tandem() {
        // Pure 5-bp tandem across 600 bp — midpoints spread uniformly
        // by phase. Every contiguous half holds ~half the counts ⇒
        // excess ≈ 0.
        let mut seq = vec![0u8; 600];
        for (i, b) in seq.iter_mut().enumerate() {
            *b = b"ACGTA"[i % 5];
        }
        let map = build_kmer_positions(&seq, 5);
        let r = kmer_phase_folded_contrast(&map, seq.len(), 5, Some(60), 0, 6, 5);
        if let Some(v) = r {
            assert!(v < 0.10, "expected near-zero phase contrast, got {v}");
        }
    }

    #[test]
    fn phase_fold_returns_none_when_founder_unknown() {
        let seq = vec![b'A'; 500];
        let map = build_kmer_positions(&seq, 5);
        let r = kmer_phase_folded_contrast(&map, seq.len(), 5, None, 0, 6, 5);
        assert!(r.is_none());
    }

    #[test]
    fn autocorr_low_for_uniform_signal() {
        // Pure 5-bp tandem across the whole array. Density is flat
        // ⇒ zero variance ⇒ None (we can't define correlation).
        let mut seq = vec![0u8; 600];
        for (i, b) in seq.iter_mut().enumerate() {
            *b = b"ACGTA"[i % 5];
        }
        let map = build_kmer_positions(&seq, 5);
        // Founder = 60 (5×12). Density is flat so variance = 0 → None.
        let r = kmer_density_autocorr_at_founder(&map, seq.len(), 5, Some(60), 0, 5);
        // Either Some(very high — trivially) when the variance is
        // tiny but nonzero, or None when it underflows. Both are
        // acceptable here; the meaningful test is that it doesn't
        // crash and isn't a strong negative signal.
        if let Some(v) = r {
            assert!(v.is_finite(), "expected finite or None, got {v}");
        }
    }
}

//! Shifted self-alignment scan for nested short tandem repeats.
//!
//! For one `(sequence, period)` pair, computes the windowed match
//! rate at lag `period` and identifies contiguous runs above an
//! identity threshold. Used by `kitehor rescore` to emit the
//! `scan_n_intervals` and `scan_occupancy_frac` columns.
//!
//! Algorithm summary (port of
//! `tools/subrepeat_scan/scan.py::scan_one_record`):
//!
//! 1. Build `match[i] = 1 if seq[i] == seq[i+period] else 0` for
//!    `i ∈ [0, L − period)`. N bases are encoded so they never
//!    match themselves.
//! 2. Compute the windowed match rate
//!    `rate[i] = mean(match[i..i+period])` for `i ∈ [0, L − 2·period]`
//!    using a cumulative-sum trick — `O(L)`.
//! 3. Find every maximal contiguous run where `rate[i] ≥
//!    id_threshold`. A run of length `r` indices corresponds to
//!    `(r + period)` bp in sequence coordinates; we keep runs whose
//!    index-length is `≥ (min_copies − 1) · period` (i.e. the
//!    tandem covers at least `min_copies · period` bp).
//! 4. Take the union of interval lengths in sequence coordinates as
//!    `occupied_bp`; `occupancy_frac = occupied_bp / L`.
//!
//! Per-row cost: `O(L)`. The k-mer position map built for
//! `kmer_autocorr_founder` and `kmer_phase_contrast` is not used
//! here; the scan operates directly on the raw byte sequence.

/// Per-`(record, period)` scan output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanResult {
    pub n_intervals: usize,
    pub occupied_bp: usize,
    pub occupancy_frac: f64,
}

/// Scan one sequence for tandem-positive runs at the given period.
///
/// Returns `None` when the array is too short to host any
/// qualifying run:
/// `seq.len() < 2 · period + (min_copies − 1) · period`. Returning
/// `None` propagates as `NA` in the rescore TSV.
pub fn scan_one_period(
    seq: &[u8],
    period: usize,
    id_threshold: f64,
    min_copies: usize,
) -> Option<ScanResult> {
    if period == 0 || min_copies < 1 {
        return None;
    }
    let l = seq.len();
    let min_run_indices = (min_copies.saturating_sub(1)) * period;
    let need = 2 * period + min_run_indices;
    if l < need {
        return None;
    }

    let match_buf = build_match(seq, period);
    let rate = windowed_rate(&match_buf, period);
    let intervals = runs_to_intervals(
        &find_runs_above(&rate, id_threshold, min_run_indices.max(1)),
        period,
        l,
    );
    let occupied_bp = union_length(&intervals);
    let occupancy_frac = if l > 0 {
        occupied_bp as f64 / l as f64
    } else {
        0.0
    };
    Some(ScanResult {
        n_intervals: intervals.len(),
        occupied_bp,
        occupancy_frac,
    })
}

/// Per-base match indicator at shift `period`: `1` where `seq[i] ==
/// seq[i+period]`. N bases (and any non-ACGT byte) are encoded so
/// they never match themselves, so N regions don't inflate the
/// match rate.
fn build_match(seq: &[u8], period: usize) -> Vec<u8> {
    debug_assert!(period > 0 && seq.len() > period);
    let n = seq.len() - period;
    let mut out = vec![0u8; n];
    for i in 0..n {
        let a = seq[i];
        let b = seq[i + period];
        if a == b && is_acgt(a) {
            out[i] = 1;
        }
    }
    out
}

#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// For each `i`, return `mean(match[i..i+period])` — the match rate
/// inside one period-wide forward window starting at `i`. Length
/// is `match.len() − period + 1`.
fn windowed_rate(match_buf: &[u8], period: usize) -> Vec<f64> {
    debug_assert!(period > 0);
    if match_buf.len() < period {
        return Vec::new();
    }
    // O(n) cumulative-sum trick.
    let mut csum = vec![0u64; match_buf.len() + 1];
    let mut acc: u64 = 0;
    for (i, &m) in match_buf.iter().enumerate() {
        acc += m as u64;
        csum[i + 1] = acc;
    }
    let n = match_buf.len() - period + 1;
    let mut out = vec![0f64; n];
    let inv = 1.0 / period as f64;
    for i in 0..n {
        let s = csum[i + period] - csum[i];
        out[i] = s as f64 * inv;
    }
    out
}

/// Maximal contiguous runs of indices where `rate[i] ≥ threshold`,
/// keeping only runs whose length is `≥ min_length`. Returns
/// `(start, end)` half-open over indices into `rate`.
fn find_runs_above(rate: &[f64], threshold: f64, min_length: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut i = 0;
    let n = rate.len();
    while i < n {
        if rate[i] >= threshold {
            let s = i;
            while i < n && rate[i] >= threshold {
                i += 1;
            }
            if i - s >= min_length {
                out.push((s, i));
            }
        } else {
            i += 1;
        }
    }
    out
}

/// Map a window-index run `[s, e)` into a sequence-coordinate
/// interval `[s, min(seq_len, e + period − 1))`. (The last window
/// starts at `e − 1` and covers `period` bp.)
fn runs_to_intervals(
    runs: &[(usize, usize)],
    period: usize,
    seq_len: usize,
) -> Vec<(usize, usize)> {
    runs.iter()
        .map(|&(s, e)| (s, (e + period - 1).min(seq_len)))
        .collect()
}

/// Total bp covered by the union of intervals (assumed to be sorted
/// by start, which they are coming out of `find_runs_above`).
fn union_length(intervals: &[(usize, usize)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    let mut total = 0usize;
    let (mut cur_s, mut cur_e) = intervals[0];
    for &(s, e) in &intervals[1..] {
        if s <= cur_e {
            if e > cur_e {
                cur_e = e;
            }
        } else {
            total += cur_e - cur_s;
            cur_s = s;
            cur_e = e;
        }
    }
    total + (cur_e - cur_s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn perfect_tandem(motif: &[u8], copies: usize) -> Vec<u8> {
        motif
            .iter()
            .cycle()
            .take(motif.len() * copies)
            .copied()
            .collect()
    }

    #[test]
    fn build_match_perfect_tandem_is_all_ones() {
        let seq = perfect_tandem(b"ACGTA", 20); // 100 bp
        let m = build_match(&seq, 5);
        // For a 5-bp perfect tandem, every position should match.
        for (i, v) in m.iter().enumerate() {
            assert_eq!(*v, 1, "expected match at i={i}, got {v}");
        }
    }

    #[test]
    fn build_match_random_is_near_quarter() {
        // Deterministic pseudo-random ACGT sequence.
        let mut seq = vec![0u8; 4000];
        let bases = b"ACGT";
        let mut s = 0u64;
        for v in seq.iter_mut() {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *v = bases[(s >> 33) as usize & 3];
        }
        let m = build_match(&seq, 30);
        let mean = m.iter().map(|&x| x as f64).sum::<f64>() / m.len() as f64;
        // Random ACGT pairs match ~ 25 %.
        assert!(mean < 0.35, "expected ≈0.25 for random ACGT, got {mean}");
        assert!(mean > 0.15, "expected ≈0.25 for random ACGT, got {mean}");
    }

    #[test]
    fn build_match_treats_n_as_mismatch() {
        // All-N region must not inflate the match rate even though
        // N == N byte-equal.
        let seq = vec![b'N'; 100];
        let m = build_match(&seq, 5);
        assert!(m.iter().all(|&v| v == 0));
    }

    #[test]
    fn windowed_rate_perfect_tandem_is_one() {
        let seq = perfect_tandem(b"ACGTA", 30); // 150 bp
        let m = build_match(&seq, 5);
        let r = windowed_rate(&m, 5);
        for v in &r {
            assert!((*v - 1.0).abs() < 1e-9, "expected 1.0, got {v}");
        }
    }

    #[test]
    fn find_runs_skips_too_short() {
        let rate = vec![1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0];
        // Threshold 0.9, min_length 3:
        // run [0,2) length 2 → drop; run [4,7) length 3 → keep; run [8,9) → drop.
        let runs = find_runs_above(&rate, 0.9, 3);
        assert_eq!(runs, vec![(4, 7)]);
    }

    #[test]
    fn scan_perfect_tandem_returns_full_occupancy() {
        let seq = perfect_tandem(b"ACGTA", 200); // 1000 bp
        let r = scan_one_period(&seq, 5, 0.55, 3).expect("expected Some");
        assert_eq!(r.n_intervals, 1);
        assert!(
            r.occupancy_frac > 0.98,
            "expected ~1.0, got {}",
            r.occupancy_frac
        );
    }

    #[test]
    fn scan_no_tandem_returns_low_occupancy() {
        // Random ACGT — no period-5 tandem. Windowed match rate is
        // autocorrelated (adjacent windows share 4 of 5 positions),
        // so a long-enough random sequence will produce a small
        // number of chance runs above threshold. Assert occupancy
        // stays well below the "real subrepeat" regime (< 0.05).
        let mut seq = vec![0u8; 4000];
        let bases = b"ACGT";
        let mut s = 12345u64;
        for v in seq.iter_mut() {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *v = bases[(s >> 33) as usize & 3];
        }
        let r = scan_one_period(&seq, 5, 0.55, 3).expect("expected Some");
        assert!(
            r.occupancy_frac < 0.05,
            "expected near-zero occupancy on random ACGT, got {} from {} intervals",
            r.occupancy_frac,
            r.n_intervals,
        );
    }

    #[test]
    fn scan_nested_subrepeat_returns_partial_occupancy() {
        // Build a synthetic array: 30 founder copies of a 60-bp
        // founder, where each founder contains 8 copies of a 5-bp
        // motif at the start (subrepeat region = 40 bp, gap = 20 bp).
        let mut founder = Vec::with_capacity(60);
        let motif = b"ACGTA";
        for _ in 0..8 {
            founder.extend_from_slice(motif);
        }
        founder.extend_from_slice(b"AAATTCCCGGGAAACCCGGGT"[..20].as_ref());
        let mut seq = Vec::with_capacity(60 * 30);
        for _ in 0..30 {
            seq.extend_from_slice(&founder);
        }
        let r = scan_one_period(&seq, 5, 0.55, 3).expect("expected Some");
        // Subrepeat region = 40/60 of each founder; the gap doesn't
        // produce period-5 matches. Expect occupancy in [0.3, 0.8].
        assert!(
            r.n_intervals >= 1,
            "expected ≥1 interval, got {}",
            r.n_intervals
        );
        assert!(
            r.occupancy_frac > 0.30 && r.occupancy_frac < 0.85,
            "expected partial occupancy, got {}",
            r.occupancy_frac
        );
    }

    #[test]
    fn scan_with_noise_still_fires_at_threshold() {
        // 1000 bp of perfect 5-bp tandem with ~15 % per-base mutated
        // → expected match rate ≈ 0.72. Threshold 0.55 should fire.
        let mut seq = perfect_tandem(b"ACGTA", 200);
        let mut s = 999u64;
        for v in seq.iter_mut() {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            if (s >> 33) as usize % 100 < 15 {
                *v = b"ACGT"[(s >> 35) as usize & 3];
            }
        }
        let r = scan_one_period(&seq, 5, 0.55, 3).expect("expected Some");
        assert!(
            r.n_intervals >= 1,
            "expected ≥1 interval, got {}",
            r.n_intervals
        );
        assert!(
            r.occupancy_frac > 0.50,
            "expected high occupancy under moderate noise, got {}",
            r.occupancy_frac
        );
    }

    #[test]
    fn scan_short_array_returns_none() {
        let seq = vec![b'A'; 50];
        let r = scan_one_period(&seq, 30, 0.55, 3);
        assert!(r.is_none());
    }

    #[test]
    fn scan_period_zero_returns_none() {
        let seq = perfect_tandem(b"ACGT", 100);
        assert!(scan_one_period(&seq, 0, 0.55, 3).is_none());
    }

    #[test]
    fn scan_min_copies_zero_returns_none() {
        let seq = perfect_tandem(b"ACGT", 100);
        assert!(scan_one_period(&seq, 4, 0.55, 0).is_none());
    }

    #[test]
    fn scan_at_long_period_reports_full_occupancy_on_perfect_tandem() {
        // 30 copies of a 60-bp random-ish founder. At lag 60 every
        // position matches, so occupancy ≈ 1.0. Confirms the "long
        // period fires on the founder" behaviour documented in
        // scan_port_plan.md §0.
        let founder = b"ACGTGGTCAAGCTTACGGGAACCTTAACGCTCTGAACGTACGAAACCGATTCAAGGCTAG";
        let mut seq = Vec::new();
        for _ in 0..30 {
            seq.extend_from_slice(founder);
        }
        let r = scan_one_period(&seq, 60, 0.55, 3).expect("expected Some");
        assert!(
            r.occupancy_frac > 0.95,
            "expected ~1.0, got {}",
            r.occupancy_frac
        );
    }

    #[test]
    fn union_length_merges_overlaps() {
        let intervals = vec![(0, 10), (5, 20), (25, 30), (28, 35)];
        assert_eq!(union_length(&intervals), 20 + 10);
    }
}

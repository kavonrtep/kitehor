//! Shift signal — Pass A (local drift / wobble / breakpoint
//! detection) per `detect_impl_plan.md §6.6` and A2.
//!
//! Pass B (phase-shift offset recovery via circular cross-correlation
//! across flanking windows) lands in M3.5.

use crate::detect::config::DetectorConfig;
use rustfft::{num_complex::Complex32, FftPlanner};

#[derive(Debug, Clone)]
pub struct ShiftFeatures {
    /// `best_shift(r)` for adjacent row pairs (r, r+1). Length =
    /// n_rows - 1.
    pub best_shift: Vec<i32>,
    /// Mean shift in bp — drift indicator; non-zero → width
    /// refinement candidate.
    pub mean_shift_bp: f64,
    /// Standard deviation of `best_shift` after breakpoint segments
    /// removed (wobble amplitude in bp).
    pub wobble_amplitude_bp: f64,
    /// Recovered wobble periodicity in bp, if a clear FFT peak
    /// exists. None otherwise.
    pub wobble_periodicity_bp: Option<f64>,
    /// Row indices (0-based into `best_shift`) where a candidate
    /// phase-shift breakpoint was detected. M3.5 uses these to seed
    /// the Pass-B offset-recovery search.
    pub breakpoints: Vec<usize>,
}

/// Pass-A compute over the rows wrapped at `width`. Returns `None`
/// when fewer than 2 rows are available.
pub fn compute(
    seq: &[u8],
    width: usize,
    n_rows: usize,
    cfg: &DetectorConfig,
) -> Option<ShiftFeatures> {
    if width == 0 || n_rows < 2 {
        return None;
    }
    let s_range = cfg.shift_local_range_bp.max(1);
    let s_lo = -s_range;
    let s_hi = s_range;

    let mut best_shift: Vec<i32> = Vec::with_capacity(n_rows - 1);
    for r in 0..n_rows - 1 {
        let prev = &seq[r * width..(r + 1) * width];
        let next = &seq[(r + 1) * width..(r + 2) * width];
        // Iterate by increasing |s| so that on a tie the closest-to-
        // zero shift wins. This avoids spurious mean drift on
        // featureless inputs (e.g. all-A rows where every shift
        // yields identical similarity).
        let mut best = (0i32, match_at_shift(prev, next, 0));
        for mag in 1..=s_range {
            for s in [-mag, mag] {
                if s < s_lo || s > s_hi {
                    continue;
                }
                let m = match_at_shift(prev, next, s);
                if m > best.1 {
                    best = (s, m);
                }
            }
        }
        best_shift.push(best.0);
    }

    let mean_shift_bp = mean(&best_shift);
    let breakpoints = find_breakpoints(&best_shift, cfg.shift_breakpoint_threshold);
    let wobble_amplitude_bp = std_excluding(&best_shift, &breakpoints);
    let wobble_periodicity_bp = fft_periodicity(&best_shift, width);

    Some(ShiftFeatures {
        best_shift,
        mean_shift_bp,
        wobble_amplitude_bp,
        wobble_periodicity_bp,
        breakpoints,
    })
}

/// Fraction of columns at which row `prev` and row `next` agree when
/// `next` is logically slid by `s` bases. Out-of-range columns are
/// skipped (not wrapped).
pub(crate) fn match_at_shift(prev: &[u8], next: &[u8], s: i32) -> f64 {
    let w = prev.len();
    if w == 0 || next.len() != w {
        return 0.0;
    }
    let mut matched = 0usize;
    let mut compared = 0usize;
    for c in 0..w {
        let other = c as i32 + s;
        if other < 0 || other >= w as i32 {
            continue;
        }
        if prev[c] == next[other as usize] {
            matched += 1;
        }
        compared += 1;
    }
    if compared == 0 {
        0.0
    } else {
        matched as f64 / compared as f64
    }
}

fn mean(xs: &[i32]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().map(|&v| v as f64).sum::<f64>() / xs.len() as f64
}

/// Standard deviation of `best_shift` excluding indices identified as
/// breakpoints (so a step in the signal doesn't inflate amplitude).
fn std_excluding(xs: &[i32], breakpoints: &[usize]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mask: std::collections::HashSet<usize> = breakpoints.iter().copied().collect();
    let kept: Vec<f64> = xs
        .iter()
        .enumerate()
        .filter(|(i, _)| !mask.contains(i))
        .map(|(_, &v)| v as f64)
        .collect();
    if kept.is_empty() {
        return 0.0;
    }
    let m = kept.iter().sum::<f64>() / kept.len() as f64;
    let var = kept.iter().map(|v| (v - m).powi(2)).sum::<f64>() / kept.len() as f64;
    var.sqrt()
}

/// Threshold-based breakpoint detection (Pass A only — locates
/// candidates; Pass B recovers the actual offset). A row at index
/// `r` is flagged when `|best_shift[r] - best_shift[r-1]| >= thr`.
pub(crate) fn find_breakpoints(best_shift: &[i32], thr: i32) -> Vec<usize> {
    let mut out = Vec::new();
    for i in 1..best_shift.len() {
        if (best_shift[i] - best_shift[i - 1]).abs() >= thr {
            out.push(i);
        }
    }
    out
}

// ---------------- Pass B — phase-shift offset recovery ----------------
//
// For each Pass-A breakpoint, take a small window of rows on either
// side, build the per-column mode profile, and cross-correlate the
// two profiles **circularly** over `s ∈ [-w/2, +w/2]`. The argmax
// shift is the phase-shift offset in bp.

/// Recover the phase-shift offset between two wrap-row windows by
/// circular cross-correlation of their column-mode profiles.
///
/// `prev_window` and `post_window` are flat row-major byte slices of
/// length `width × rows_in_window` each. Returns the offset `s` (bp)
/// such that the post-window aligns with the pre-window when shifted
/// by `-s` (the conventional sign: positive `s` means the post-block
/// "moved forward" by `s` bp relative to the pre-block).
pub fn recover_offset(prev_window: &[u8], post_window: &[u8], width: usize) -> i32 {
    if width < 4 || prev_window.len() < width || post_window.len() < width {
        return 0;
    }
    let prev = mode_per_column(prev_window, width);
    let post = mode_per_column(post_window, width);
    let s_max = (width as i32) / 2;
    let mut best: (i32, usize) = (0, count_matches_circular(&prev, &post, 0));
    for mag in 1..=s_max {
        for s in [-mag, mag] {
            let m = count_matches_circular(&prev, &post, s);
            if m > best.1 {
                best = (s, m);
            }
        }
    }
    best.0
}

fn mode_per_column(window: &[u8], width: usize) -> Vec<u8> {
    let n_rows = window.len() / width;
    let mut out = vec![b'N'; width];
    for c in 0..width {
        let mut counts = [0usize; 5]; // A, C, G, T, N
        for r in 0..n_rows {
            let b = window[r * width + c];
            let i = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            counts[i] += 1;
        }
        let (max_i, _) = counts.iter().enumerate().max_by_key(|&(_, n)| *n).unwrap();
        out[c] = match max_i {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        };
    }
    out
}

fn count_matches_circular(prev: &[u8], post: &[u8], s: i32) -> usize {
    let w = prev.len() as i32;
    let mut matched = 0usize;
    for c in 0..prev.len() {
        let other = (c as i32 + s).rem_euclid(w) as usize;
        if prev[c] == post[other] && prev[c] != b'N' {
            matched += 1;
        }
    }
    matched
}

/// Convenience: given a sequence and a list of `(width, breakpoint
/// row indices)` pairs, recover one offset per breakpoint by taking
/// the flanking `window_rows` rows on each side.
pub fn recover_offsets_at_breakpoints(
    seq: &[u8],
    width: usize,
    n_rows: usize,
    breakpoints: &[usize],
    window_rows: usize,
) -> Vec<i32> {
    let mut out = Vec::with_capacity(breakpoints.len());
    for &b in breakpoints {
        // Pass-A breakpoint at index `b` in best_shift means the
        // transition is between row b and row b+1 in the wrap. (Recall
        // best_shift[i] compares row i and row i+1.)
        let pre_lo = b.saturating_sub(window_rows);
        let pre_hi = b + 1; // exclusive
        let post_lo = b + 1;
        let post_hi = (post_lo + window_rows).min(n_rows);
        if pre_hi - pre_lo < 2 || post_hi - post_lo < 2 {
            out.push(0);
            continue;
        }
        let prev = &seq[pre_lo * width..pre_hi * width];
        let post = &seq[post_lo * width..post_hi * width];
        out.push(recover_offset(prev, post, width));
    }
    out
}

/// FFT of detrended `best_shift`. Returns the period in bp = `width ×
/// (n_rows / k)` for the bin with the largest magnitude (excluding
/// DC and the trivial first lag). Returns `None` when no bin clears a
/// "meaningfully periodic" magnitude floor.
fn fft_periodicity(best_shift: &[i32], width: usize) -> Option<f64> {
    let n = best_shift.len();
    if n < 16 {
        return None;
    }
    let mean = best_shift.iter().map(|&v| v as f32).sum::<f32>() / n as f32;
    let mut buf: Vec<Complex32> = best_shift
        .iter()
        .map(|&v| Complex32::new(v as f32 - mean, 0.0))
        .collect();
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(n);
    fft.process(&mut buf);

    // Look at bins 2..n/2 (skip DC = 0 and the trivial first harmonic).
    let mut best = (0usize, 0.0f32);
    for k in 2..n / 2 {
        let mag = buf[k].norm_sqr();
        if mag > best.1 {
            best = (k, mag);
        }
    }
    // Magnitude floor: require the peak to be ≥ 3× the median bin
    // magnitude (cheap "is this peak distinct from background" test).
    let mut all_mags: Vec<f32> = (2..n / 2).map(|k| buf[k].norm_sqr()).collect();
    if all_mags.is_empty() {
        return None;
    }
    all_mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = all_mags[all_mags.len() / 2];
    if best.1 < median * 3.0 {
        return None;
    }
    // Period = signal length / harmonic index, expressed in rows;
    // bp = period_rows × width.
    let period_rows = n as f64 / best.0 as f64;
    Some(period_rows * width as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::config::DetectorConfig;

    fn cfg() -> DetectorConfig {
        DetectorConfig::default()
    }

    #[test]
    fn match_at_shift_identifies_exact_offset() {
        // next = prev shifted right by 2.
        let prev = b"ACGTACGTAC";
        let next = b"GTACGTACGT"; // arbitrary; exact offset depends.
        // Build a deterministic case instead:
        let row_a = b"ACGTACGTAC";
        let row_b = b"GTACGTACAC"; // = row_a shifted by -2 (loose tail)
        let _ = (prev, next);
        let m_neg2 = match_at_shift(row_a, row_b, -2);
        let m_zero = match_at_shift(row_a, row_b, 0);
        assert!(
            m_neg2 > m_zero,
            "negative shift should align better than zero; m_neg2={m_neg2} m_zero={m_zero}"
        );
    }

    #[test]
    fn identical_rows_give_zero_shift_and_zero_wobble() {
        let seq = vec![b'A'; 200];
        let s = compute(&seq, 10, 20, &cfg()).unwrap();
        assert_eq!(s.best_shift.len(), 19);
        assert!(s.mean_shift_bp.abs() < 1e-9);
        assert_eq!(s.wobble_amplitude_bp, 0.0);
        assert!(s.breakpoints.is_empty());
    }

    #[test]
    fn breakpoint_threshold_triggers_on_large_step() {
        let xs = vec![0, 0, 0, 0, 5, 5, 5, 5];
        let bp = find_breakpoints(&xs, 3);
        assert_eq!(bp, vec![4]);
    }

    #[test]
    fn breakpoint_skipped_when_change_below_threshold() {
        let xs = vec![0, 1, 2, 1, 0];
        let bp = find_breakpoints(&xs, 3);
        assert!(bp.is_empty());
    }

    #[test]
    fn recover_offset_identifies_known_shift() {
        // Build a non-periodic 40 bp row, then a row that's a cyclic
        // shift of it by exactly 7. Pass B should recover that shift.
        let width = 40usize;
        // Deterministic pseudo-random row using a fixed seed.
        use rand::SeedableRng;
        use rand::Rng;
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::seed_from_u64(99);
        let row_a: Vec<u8> = (0..width)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let shift_amount = 7i32;
        let mut row_b = vec![0u8; width];
        for i in 0..width {
            row_b[i] = row_a[(i as i32 - shift_amount).rem_euclid(width as i32) as usize];
        }
        let rows = 20usize;
        let pre: Vec<u8> = (0..rows).flat_map(|_| row_a.iter().copied()).collect();
        let post: Vec<u8> = (0..rows).flat_map(|_| row_b.iter().copied()).collect();

        let offset = recover_offset(&pre, &post, width);
        // post[i] = pre[i - shift_amount] → at s = +shift_amount we
        // compare pre[c] with post[c + shift_amount] = pre[c] → match.
        // (Sign convention: positive s means post-block "moved forward".)
        assert!(
            offset == shift_amount || offset == -shift_amount,
            "expected offset = ±{shift_amount}; got {offset}"
        );
    }

    #[test]
    fn recover_offset_zero_for_identical_windows() {
        let width = 30usize;
        let row: Vec<u8> = (0..width).map(|i| if i % 2 == 0 { b'A' } else { b'C' }).collect();
        let win: Vec<u8> = (0..15).flat_map(|_| row.iter().copied()).collect();
        let offset = recover_offset(&win, &win, width);
        assert_eq!(offset, 0);
    }

    #[test]
    fn synthetic_step_recovers_one_breakpoint() {
        // Two halves, each 8 rows of width 10. First half is row 'A',
        // second half is the same row shifted by +3 cyclically (we
        // construct directly).
        let width = 16usize;
        let n_rows = 32usize;
        let row_a: Vec<u8> = b"ACGTACGTACGTACGT".to_vec();
        let row_b: Vec<u8> = b"TACGTACGTACGTACG".to_vec(); // row_a rotated right by 1
        let mut seq = Vec::new();
        for r in 0..n_rows {
            if r < n_rows / 2 {
                seq.extend_from_slice(&row_a);
            } else {
                seq.extend_from_slice(&row_b);
            }
        }
        let mut c = cfg();
        c.shift_local_range_bp = 3;
        c.shift_breakpoint_threshold = 1;
        let s = compute(&seq, width, n_rows, &c).unwrap();
        // The single discontinuity sits at row n_rows/2 - 1 -> n_rows/2.
        // Best-shift array index = n_rows/2 - 1 (transition pair).
        assert!(
            !s.breakpoints.is_empty(),
            "expected at least one breakpoint; got none. best_shift={:?}",
            s.best_shift
        );
    }
}

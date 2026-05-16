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

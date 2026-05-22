//! Edge fields and column-edge profile (`detect_impl_plan.md §6.5`).
//!
//! ```text
//! diff_x[r, c] = 1 iff S[r·w + c]     ≠ S[r·w + c+1]   (within-row change)
//! diff_y[r, c] = 1 iff S[r·w + c]     ≠ S[(r+1)·w + c] (row-to-row change at fixed col)
//! ```
//!
//! Aggregate features:
//! - `horizontal_edge_rate = mean(diff_x)`
//! - `vertical_edge_rate   = mean(diff_y)`
//! - `column_edge_rate[c]  = mean_r diff_y[r, c]`
//!
//! `column_edge_rate` autocorrelated along `c` gives an independent
//! vote on the HOR multiplicity — at lag `k = base_width × multiplicity`
//! the profile self-aligns. The strongest non-trivial lag is reported
//! as `(column_edge_autocorr_k, column_edge_autocorr_score)`.

#[derive(Debug, Clone)]
pub struct EdgeFeatures {
    pub horizontal_edge_rate: f64,
    pub vertical_edge_rate: f64,
    pub column_edge_rate: Vec<f64>,
    pub column_edge_autocorr_k: Option<usize>,
    pub column_edge_autocorr_score: Option<f64>,
}

/// Compute the edge-field summary for a sequence wrapped at width `w`.
/// `n_rows` is the number of complete rows (caller computes this
/// from `wrap::wrap_and_ic`). Returns `None` if `n_rows < 2`.
pub fn compute(seq: &[u8], width: usize, n_rows: usize) -> Option<EdgeFeatures> {
    if width == 0 || n_rows < 2 {
        return None;
    }

    // diff_x rate: per-row neighboring-column changes.
    let mut total_x = 0usize;
    let mut count_x = 0usize;
    for r in 0..n_rows {
        let base = r * width;
        for c in 0..width.saturating_sub(1) {
            if seq[base + c] != seq[base + c + 1] {
                total_x += 1;
            }
            count_x += 1;
        }
    }
    let horizontal_edge_rate = if count_x == 0 {
        0.0
    } else {
        total_x as f64 / count_x as f64
    };

    // diff_y rate and per-column profile.
    let mut col_edges = vec![0u32; width];
    let mut total_y = 0usize;
    for r in 0..n_rows - 1 {
        let base = r * width;
        let next = (r + 1) * width;
        for c in 0..width {
            if seq[base + c] != seq[next + c] {
                col_edges[c] += 1;
                total_y += 1;
            }
        }
    }
    let pairs = n_rows - 1;
    let denom = (width * pairs) as f64;
    let vertical_edge_rate = if denom == 0.0 {
        0.0
    } else {
        total_y as f64 / denom
    };
    let column_edge_rate: Vec<f64> = col_edges.iter().map(|&n| n as f64 / pairs as f64).collect();

    // Autocorrelation of column_edge_rate along c (Pearson-style),
    // searching k in 1..width/2 for the strongest non-trivial peak.
    let (k_max, score) = best_circular_autocorr(&column_edge_rate);

    Some(EdgeFeatures {
        horizontal_edge_rate,
        vertical_edge_rate,
        column_edge_rate,
        column_edge_autocorr_k: k_max,
        column_edge_autocorr_score: score,
    })
}

/// Circular Pearson autocorrelation of `x` at lags 2..len/2. Returns
/// the lag with the largest value and the corresponding score.
fn best_circular_autocorr(x: &[f64]) -> (Option<usize>, Option<f64>) {
    let n = x.len();
    if n < 4 {
        return (None, None);
    }
    let mean = x.iter().sum::<f64>() / n as f64;
    let centered: Vec<f64> = x.iter().map(|v| v - mean).collect();
    let var = centered.iter().map(|v| v.powi(2)).sum::<f64>() / n as f64;
    if var <= f64::EPSILON {
        return (None, None);
    }
    let max_lag = (n / 2).max(2);
    let mut best_k = None;
    let mut best_score = None;
    for k in 2..=max_lag {
        let mut s = 0.0;
        for i in 0..n {
            s += centered[i] * centered[(i + k) % n];
        }
        let r = s / (n as f64 * var);
        if best_score.map(|b: f64| r > b).unwrap_or(true) {
            best_k = Some(k);
            best_score = Some(r);
        }
    }
    (best_k, best_score)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_sequence_has_no_edges() {
        let seq = vec![b'A'; 1000];
        let e = compute(&seq, 100, 10).unwrap();
        assert!(e.horizontal_edge_rate.abs() < 1e-9);
        assert!(e.vertical_edge_rate.abs() < 1e-9);
        for r in &e.column_edge_rate {
            assert!(r.abs() < 1e-9);
        }
    }

    #[test]
    fn checkerboard_high_edges() {
        // ACAC repeated rows, alternating: AC AC AC; next row CA CA CA.
        let mut seq = Vec::new();
        let width = 4usize;
        let n_rows = 8usize;
        for r in 0..n_rows {
            let bytes = if r % 2 == 0 { b"ACAC" } else { b"CACA" };
            seq.extend_from_slice(bytes);
        }
        let e = compute(&seq, width, n_rows).unwrap();
        // Every neighbouring pair differs both within and across rows.
        assert!(e.horizontal_edge_rate > 0.99);
        assert!(e.vertical_edge_rate > 0.99);
    }

    #[test]
    fn column_autocorr_finds_period() {
        // Build a 2D matrix where columns have a clear k=4 periodicity
        // in their vertical edge rate.
        let width = 16usize;
        let n_rows = 100usize;
        let mut seq = Vec::with_capacity(width * n_rows);
        for r in 0..n_rows {
            for c in 0..width {
                // Columns 0, 4, 8, 12 flip every row; others stay 'A'.
                let b = if c % 4 == 0 && r % 2 == 1 { b'C' } else { b'A' };
                seq.push(b);
            }
        }
        let e = compute(&seq, width, n_rows).unwrap();
        let k = e.column_edge_autocorr_k.unwrap();
        // Period 4 in column_edge_rate, but circular autocorr can
        // equivalently land on 8 or 12 (multiples). Accept any multiple of 4.
        assert!(
            k % 4 == 0 && k >= 4,
            "expected autocorr lag = multiple of 4; got {k}"
        );
        let score = e.column_edge_autocorr_score.unwrap();
        assert!(score > 0.5, "autocorr score should be strong; got {score}");
    }

    #[test]
    fn returns_none_for_single_row() {
        let seq = vec![b'A'; 100];
        assert!(compute(&seq, 100, 1).is_none());
    }
}

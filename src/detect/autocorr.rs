//! Row autocorrelation `R(k)` (`detect_impl_plan.md §6.4`).
//!
//! ```text
//! R(k) = mean over i of dot(emb[i], emb[i + k])    for k = 1..K
//! ```
//!
//! Returns the curve plus the canonical summary: `best_k =
//! argmax_{k ≥ 2} R(k)`, `best_lag_score = R(best_k)`, and `r_lag1 =
//! R(1)`.

use crate::detect::embed::{dot, RowEmbedding};

#[derive(Debug, Clone)]
pub struct AutocorrSummary {
    /// `R(k)` for k = 1..=K (length K, 0-indexed).
    pub r_k: Vec<f64>,
    /// k ≥ 2 with the largest `R(k)`; `None` if K < 2 or all R(k) are
    /// uniformly low (no lag with meaningful similarity).
    pub best_lag: Option<usize>,
    /// `R(best_lag)`.
    pub best_lag_score: Option<f64>,
    /// `R(1)`.
    pub r_lag1: Option<f64>,
}

/// Compute `R(k)` for `k = 1..=k_max`.
///
/// `k_max` is clamped to `embeddings.len() - 1` so we never read out
/// of bounds. Empty input returns a summary with `None` everywhere.
pub fn compute(embeddings: &[RowEmbedding], k_max: usize) -> AutocorrSummary {
    let n = embeddings.len();
    let mut r_k = Vec::new();
    if n < 2 {
        return AutocorrSummary {
            r_k,
            best_lag: None,
            best_lag_score: None,
            r_lag1: None,
        };
    }
    let k_top = k_max.min(n - 1);
    for k in 1..=k_top {
        let pairs = n - k;
        let sum: f64 = (0..pairs)
            .map(|i| dot(&embeddings[i], &embeddings[i + k]))
            .sum();
        r_k.push(sum / pairs as f64);
    }
    let r_lag1 = r_k.first().copied();
    let (best_lag, best_lag_score) = if r_k.len() >= 2 {
        // argmax over k ≥ 2 (i.e. indices ≥ 1 in r_k since r_k[0] = R(1)).
        let mut best_idx = 1usize;
        let mut best_val = r_k[1];
        for (i, &v) in r_k.iter().enumerate().skip(2) {
            if v > best_val {
                best_val = v;
                best_idx = i;
            }
        }
        (Some(best_idx + 1), Some(best_val))
    } else {
        (None, None)
    };
    AutocorrSummary {
        r_k,
        best_lag,
        best_lag_score,
        r_lag1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::config::DetectorConfig;
    use crate::detect::embed::embed_rows;

    #[test]
    fn empty_input_yields_none() {
        let s = compute(&[], 30);
        assert!(s.best_lag.is_none());
        assert!(s.r_lag1.is_none());
        assert!(s.r_k.is_empty());
    }

    #[test]
    fn r_lag1_equals_one_for_identical_rows() {
        // All-A sequence wrapped to width 10 produces identical
        // embeddings → R(1) = 1.0 and R(k) = 1.0 for every k.
        let seq = vec![b'A'; 200];
        let cfg = DetectorConfig::default();
        let embs = embed_rows(&seq, 10, &cfg);
        assert_eq!(embs.len(), 20);
        let s = compute(&embs, 5);
        for v in &s.r_k {
            assert!((v - 1.0).abs() < 1e-6, "R(k) should be 1.0; got {v}");
        }
        assert!((s.r_lag1.unwrap() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn clean_hor_k4_peaks_at_lag_4() {
        // Four distinct slots (all-A, all-C, all-G, all-T) × 50 copies,
        // wrapped to width = monomer length (100). At lag 4, rows are
        // identical (same slot); at lag 1, 2, 3, rows are different.
        let mut seq = Vec::new();
        for _ in 0..50 {
            for c in [b'A', b'C', b'G', b'T'] {
                seq.extend(std::iter::repeat(c).take(100));
            }
        }
        let cfg = DetectorConfig::default();
        let embs = embed_rows(&seq, 100, &cfg);
        assert_eq!(embs.len(), 200);
        let s = compute(&embs, 12);
        assert_eq!(s.best_lag.unwrap(), 4, "best lag should be k=4; r_k={:?}", s.r_k);
        // R(4) should be ~1.0; R(1) < R(4) by a noticeable margin.
        assert!(s.r_k[3] > s.r_k[0] + 0.5);
    }

    #[test]
    fn k_max_clamped_to_n_minus_1() {
        let seq = vec![b'A'; 50];
        let cfg = DetectorConfig::default();
        let embs = embed_rows(&seq, 10, &cfg);
        assert_eq!(embs.len(), 5);
        let s = compute(&embs, 100);
        assert_eq!(s.r_k.len(), 4); // n_rows - 1
    }
}

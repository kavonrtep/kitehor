//! Oriented k-mer row embeddings (`detect_impl_plan.md §6.3`, A6).
//!
//! Default `k = 4` gives 4⁴ = 256 oriented k-mers. **Not**
//! reverse-complement-canonicalised: input is assumed canonically
//! oriented upstream, and strand-aware comparison is a v2 feature.
//! Oriented k-mers keep the embedding sensitive to inversions when
//! that mode arrives.
//!
//! Each row's embedding is the L2-normalised count vector across all
//! 4ᵏ k-mers. k-mers containing N are skipped (don't contribute to
//! any bin). Row similarity is the dot product, which equals 1.0 for
//! identical rows and is ≥ 0 for any pair of count-derived vectors.

use crate::detect::config::DetectorConfig;

/// One row's L2-normalised k-mer count vector. Length = 4ᵏ.
pub type RowEmbedding = Vec<f32>;

/// Build embeddings for every wrap-row at the given width.
/// Returns a `Vec` of `n_rows = seq.len() / width` rows; trailing
/// partial row is dropped exactly as in `wrap::wrap_and_ic`.
pub fn embed_rows(seq: &[u8], width: usize, cfg: &DetectorConfig) -> Vec<RowEmbedding> {
    let k = cfg.embedding_k;
    let dim = 4usize.pow(k as u32);
    if width < k || seq.len() < width {
        return Vec::new();
    }
    let n_rows = seq.len() / width;
    let mut out = Vec::with_capacity(n_rows);
    for r in 0..n_rows {
        let start = r * width;
        let row = &seq[start..start + width];
        out.push(embed_row(row, k, dim));
    }
    out
}

/// Embed a single row.
fn embed_row(row: &[u8], k: usize, dim: usize) -> RowEmbedding {
    let mut counts = vec![0u32; dim];
    let mut idx: u32 = 0;
    let mut filled = 0usize;
    for &b in row {
        let code = base_code(b);
        match code {
            Some(c) => {
                idx = ((idx << 2) | c) & ((dim as u32) - 1);
                filled = filled.saturating_add(1).min(k);
                if filled == k {
                    counts[idx as usize] += 1;
                }
            }
            None => {
                // N or other → reset the rolling window
                idx = 0;
                filled = 0;
            }
        }
    }
    // L2-normalise (drop the squared magnitude into v.iter().map(...).sum::<f64>()).
    let sumsq: f64 = counts.iter().map(|&c| (c as f64).powi(2)).sum();
    let norm = sumsq.sqrt();
    if norm == 0.0 {
        return vec![0.0; dim];
    }
    let inv = 1.0 / norm;
    counts.iter().map(|&c| (c as f64 * inv) as f32).collect()
}

#[inline]
fn base_code(b: u8) -> Option<u32> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Dot product of two row embeddings of equal length.
pub fn dot(a: &RowEmbedding, b: &RowEmbedding) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut s = 0.0;
    for i in 0..a.len() {
        s += a[i] as f64 * b[i] as f64;
    }
    s
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    fn cfg(k: usize) -> DetectorConfig {
        let mut c = DetectorConfig::default();
        c.embedding_k = k;
        c
    }

    #[test]
    fn embed_dim_is_4_to_the_k() {
        let v = embed_row(b"ACGTACGTACGT", 4, 256);
        assert_eq!(v.len(), 256);
        let v3 = embed_row(b"ACGTACGTACGT", 3, 64);
        assert_eq!(v3.len(), 64);
    }

    #[test]
    fn all_a_row_concentrates_at_aaaa_bin() {
        // k=4, all-A row of length 10 → 7 instances of AAAA k-mer.
        // L2-norm makes the bin = 1.0; others 0.
        let v = embed_row(b"AAAAAAAAAA", 4, 256);
        // AAAA = index 0 in base-4 encoding (A=0).
        assert!((v[0] - 1.0).abs() < 1e-6, "AAAA bin should be 1.0");
        for (i, &x) in v.iter().enumerate() {
            if i != 0 {
                assert!(x.abs() < 1e-6, "bin {i} should be 0, got {x}");
            }
        }
    }

    #[test]
    fn identical_rows_dot_to_one() {
        let a = embed_row(b"ACGTACGTACGTACGT", 4, 256);
        let d = dot(&a, &a);
        assert!(
            (d - 1.0).abs() < 1e-6,
            "identical rows must dot to 1.0; got {d}"
        );
    }

    #[test]
    fn n_resets_the_kmer_window() {
        // Row with N in the middle. The k-mers spanning the N are
        // dropped; surrounding ones contribute.
        let v = embed_row(b"AAAANAAAA", 4, 256);
        // Pre-N: positions 0..3 give one AAAA (at offset 3).
        // Post-N: positions 5..8 give one AAAA (at offset 8).
        // Total AAAA count = 2 → after L2 normalisation = 1.0.
        assert!((v[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn empty_or_all_n_row_returns_zero_vec() {
        let v = embed_row(b"NNNN", 4, 256);
        assert!(v.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn embed_rows_drops_trailing_partial() {
        let cfg = cfg(4);
        // 13 bases, width 5 → 2 full rows.
        let seq = b"ACGTAACGTACAC";
        let rows = embed_rows(seq, 5, &cfg);
        assert_eq!(rows.len(), 2);
        for r in &rows {
            assert_eq!(r.len(), 256);
        }
    }

    #[test]
    fn distinct_rows_have_lower_similarity() {
        // Random-ish row A vs random-ish row B should have similarity
        // < 1.0. Using fixed deterministic strings.
        let a = embed_row(b"ACGTACGTACGTACGTACGTACGTACGTACGT", 4, 256);
        let b = embed_row(b"TGCATGCATGCATGCATGCATGCATGCATGCA", 4, 256);
        let s = dot(&a, &b);
        assert!(
            s < 0.95,
            "expected distinct rows to similar < 0.95; got {s}"
        );
    }
}

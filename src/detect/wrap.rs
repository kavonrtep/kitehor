//! Wrap a sequence to width *w* and compute background-corrected
//! column IC (`detect_impl_plan.md §6.2`, A5).
//!
//! Wrap semantics (pinned in the plan):
//! - Bases are already A/C/G/T/N (the existing `crate::io::load_fasta`
//!   normaliser handles uppercasing and non-ACGT → N).
//! - At width `w`, full rows are `n_rows = seq.len() / w`. The
//!   trailing partial row is **dropped** for feature extraction;
//!   `length_bp` in the property table still reflects the original.
//! - Widths producing fewer than `cfg.min_rows_per_width` rows are
//!   skipped (`None` returned).
//! - Per-column IC is computed over A/C/G/T only — Ns are excluded
//!   from the denominator. `IC = Σ p_b · log2(p_b / q_b)`.
//! - Array-wide background `q_b` is computed once per array, over
//!   A/C/G/T only. Any `q_b == 0` falls back to the uniform 0.25 so
//!   `log(0)` never appears.

use crate::detect::config::DetectorConfig;

/// Array-wide A/C/G/T frequencies, computed once per array and
/// reused across every tested width.
#[derive(Debug, Clone)]
pub struct Background {
    pub q: [f64; 4],   // A, C, G, T
    pub n_acgt: usize, // total non-N bases
}

impl Background {
    pub fn compute(seq: &[u8]) -> Self {
        let mut counts = [0usize; 4];
        let mut n_acgt = 0usize;
        for &b in seq {
            let i = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => continue,
            };
            counts[i] += 1;
            n_acgt += 1;
        }
        let mut q = [0.0; 4];
        if n_acgt > 0 {
            for i in 0..4 {
                q[i] = counts[i] as f64 / n_acgt as f64;
                if q[i] == 0.0 {
                    q[i] = 0.25; // fallback to avoid log(0) downstream
                }
            }
        } else {
            q = [0.25; 4];
        }
        Background { q, n_acgt }
    }
}

/// Per-width result. Returns `None` when the width can't produce
/// enough complete rows.
#[derive(Debug, Clone)]
pub struct WrapStats {
    pub n_rows: usize,
    /// IC per column, length = width.
    pub column_ic: Vec<f64>,
    pub mean_column_ic: f64,
    /// Fraction of columns with `column_ic >= cfg.ic_threshold_min`.
    pub fraction_conserved: f64,
}

pub fn wrap_and_ic(
    seq: &[u8],
    width: usize,
    bg: &Background,
    cfg: &DetectorConfig,
) -> Option<WrapStats> {
    if width == 0 {
        return None;
    }
    let n_rows = seq.len() / width;
    if n_rows < cfg.min_rows_per_width {
        return None;
    }

    let mut column_ic = vec![0.0; width];
    for c in 0..width {
        let mut counts = [0usize; 4];
        let mut n_n = 0usize;
        for r in 0..n_rows {
            let b = seq[r * width + c];
            match b {
                b'A' => counts[0] += 1,
                b'C' => counts[1] += 1,
                b'G' => counts[2] += 1,
                b'T' => counts[3] += 1,
                _ => n_n += 1,
            }
        }
        let n_acgt = n_rows - n_n;
        if n_acgt == 0 {
            // All-N column carries no information; IC = 0 by definition.
            column_ic[c] = 0.0;
            continue;
        }
        let denom = n_acgt as f64;
        let mut ic = 0.0;
        for (i, &count) in counts.iter().enumerate() {
            let p = count as f64 / denom;
            if p > 0.0 {
                let q = bg.q[i];
                // q is guaranteed > 0 by Background::compute's fallback.
                ic += p * (p / q).log2();
            }
        }
        column_ic[c] = ic;
    }

    let mean_column_ic = column_ic.iter().sum::<f64>() / width as f64;
    let n_conserved = column_ic
        .iter()
        .filter(|&&x| x >= cfg.ic_threshold_min)
        .count();
    let fraction_conserved = n_conserved as f64 / width as f64;

    Some(WrapStats {
        n_rows,
        column_ic,
        mean_column_ic,
        fraction_conserved,
    })
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;
    use crate::detect::config::DetectorConfig;

    fn cfg_min_rows(n: usize) -> DetectorConfig {
        let mut c = DetectorConfig::default();
        c.min_rows_per_width = n;
        c
    }

    #[test]
    fn background_for_uniform_sequence() {
        let seq = b"ACGTACGTACGT";
        let bg = Background::compute(seq);
        assert_eq!(bg.n_acgt, 12);
        for i in 0..4 {
            assert!((bg.q[i] - 0.25).abs() < 1e-9);
        }
    }

    #[test]
    fn background_handles_all_n() {
        let seq = b"NNNNN";
        let bg = Background::compute(seq);
        assert_eq!(bg.n_acgt, 0);
        // All q fall back to 0.25.
        for i in 0..4 {
            assert!((bg.q[i] - 0.25).abs() < 1e-9);
        }
    }

    #[test]
    fn background_fallback_for_zero_base() {
        let seq = b"AAAAAAAAAA"; // C, G, T all absent
        let bg = Background::compute(seq);
        // q[A] = 1.0, others fall back to 0.25.
        assert!((bg.q[0] - 1.0).abs() < 1e-9);
        for i in 1..4 {
            assert!((bg.q[i] - 0.25).abs() < 1e-9);
        }
    }

    #[test]
    fn wrap_drops_trailing_partial_row() {
        let cfg = cfg_min_rows(2);
        // 7 bases, width 3 → 2 complete rows, 1 trailing leftover.
        let seq = b"ACGACGT";
        let bg = Background::compute(seq);
        let s = wrap_and_ic(seq, 3, &bg, &cfg).unwrap();
        assert_eq!(s.n_rows, 2);
        assert_eq!(s.column_ic.len(), 3);
    }

    #[test]
    fn wrap_returns_none_when_not_enough_rows() {
        let cfg = cfg_min_rows(10);
        let seq = b"ACGTACGT";
        let bg = Background::compute(seq);
        // 8 / 4 = 2 rows < 10 → None.
        assert!(wrap_and_ic(seq, 4, &bg, &cfg).is_none());
    }

    #[test]
    fn all_a_column_against_uniform_bg_is_2_bits() {
        // 100 rows of width 4, every column always 'A'.
        let seq: Vec<u8> = vec![b'A'; 400];
        // Bg from this seq is q[A]=1.0; the other bases fall back to
        // 0.25 each. To get the textbook "2 bits" we want a uniform
        // background, so construct a separate balanced bg.
        let bg = Background {
            q: [0.25, 0.25, 0.25, 0.25],
            n_acgt: 1,
        };
        let cfg = cfg_min_rows(8);
        let s = wrap_and_ic(&seq, 4, &bg, &cfg).unwrap();
        for ic in &s.column_ic {
            assert!((ic - 2.0).abs() < 1e-9, "expected 2.0, got {ic}");
        }
        assert!((s.mean_column_ic - 2.0).abs() < 1e-9);
        assert_eq!(s.fraction_conserved, 1.0); // all conserved
    }

    #[test]
    fn random_sequence_low_ic() {
        // Pseudo-random ACGT sequence — every column is roughly
        // uniform, so per-column IC should be near 0.
        use rand::Rng;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let n = 4000;
        let mut seq = Vec::with_capacity(n);
        for _ in 0..n {
            seq.push(b"ACGT"[rng.random_range(0..4)]);
        }
        let bg = Background::compute(&seq);
        let cfg = cfg_min_rows(10);
        let s = wrap_and_ic(&seq, 40, &bg, &cfg).unwrap();
        assert!(
            s.mean_column_ic < 0.05,
            "random sequence column IC should be ~0, got {}",
            s.mean_column_ic
        );
    }

    #[test]
    fn n_excluded_from_denominator() {
        // 20 rows × width 2; column 0 is all 'N', column 1 is all 'A'.
        let mut seq = Vec::with_capacity(40);
        for _ in 0..20 {
            seq.push(b'N');
            seq.push(b'A');
        }
        let bg = Background {
            q: [0.25, 0.25, 0.25, 0.25],
            n_acgt: 1,
        };
        let cfg = cfg_min_rows(8);
        let s = wrap_and_ic(&seq, 2, &bg, &cfg).unwrap();
        assert_eq!(s.column_ic[0], 0.0, "all-N column IC must be 0");
        assert!(
            (s.column_ic[1] - 2.0).abs() < 1e-9,
            "all-A column IC must be 2.0"
        );
    }

    #[test]
    fn hor_matrix_yields_conserved_columns() {
        // Synthetic clean HOR: 4 slots × 100 bp × 50 copies. At
        // width = 100, every column has high IC (one base each).
        let mut slots = Vec::with_capacity(4);
        for c in [b'A', b'C', b'G', b'T'] {
            slots.push(vec![c; 100]);
        }
        let mut seq = Vec::with_capacity(4 * 100 * 50);
        for _ in 0..50 {
            for s in &slots {
                seq.extend_from_slice(s);
            }
        }
        let bg = Background::compute(&seq);
        let cfg = cfg_min_rows(10);
        // Width = 100 (one slot per row): all-A, all-C, all-G, all-T rows
        // mix uniformly → column entropy looks like background → low IC.
        let s100 = wrap_and_ic(&seq, 100, &bg, &cfg).unwrap();
        // Width = 400 (one HOR unit per row): each column is constant
        // (its single base) → high IC.
        let s400 = wrap_and_ic(&seq, 400, &bg, &cfg).unwrap();
        assert!(
            s400.mean_column_ic > s100.mean_column_ic + 0.5,
            "HOR-unit width should have higher IC than base width on a degenerate-by-slot dataset; got base={} unit={}",
            s100.mean_column_ic, s400.mean_column_ic
        );
        assert!(s400.fraction_conserved > 0.9);
    }
}

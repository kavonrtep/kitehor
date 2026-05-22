//! Local-irregularity score (`detect_impl_plan.md §6.8`,
//! `detect_spec.md §7.9`).
//!
//! Splits the array into blocks of `B` rows where
//! `B = max(block_size_rows_min, n_rows / 50)`. For each block,
//! computes mean column IC (the simplest per-block feature). The
//! irregularity score is the std of per-block IC divided by the
//! global mean IC — a unit-free "how much do block features
//! disagree" measure.
//!
//! Threshold-based use (M4): when irregularity exceeds a small
//! cut-off, an otherwise-HOR call is downgraded to `irregular_HOR`.

use crate::detect::config::DetectorConfig;

/// Compute the irregularity score for a sequence wrapped at `width`.
///
/// Returns `None` when fewer than two blocks fit (too few rows for a
/// variance estimate to be meaningful).
pub fn compute(
    seq: &[u8],
    width: usize,
    bg: &crate::detect::wrap::Background,
    cfg: &DetectorConfig,
) -> Option<f64> {
    if width == 0 {
        return None;
    }
    let n_rows = seq.len() / width;
    if n_rows < 2 * cfg.block_size_rows_min {
        return None;
    }
    let block_rows = cfg.block_size_rows_min.max(n_rows / 50);
    let n_blocks = n_rows / block_rows;
    if n_blocks < 2 {
        return None;
    }

    let mut block_ics: Vec<f64> = Vec::with_capacity(n_blocks);
    let mut cfg_block = cfg.clone();
    cfg_block.min_rows_per_width = 2;
    for b in 0..n_blocks {
        let row_lo = b * block_rows;
        let row_hi = row_lo + block_rows;
        let slice = &seq[row_lo * width..row_hi * width];
        if let Some(stats) = crate::detect::wrap::wrap_and_ic(slice, width, bg, &cfg_block) {
            block_ics.push(stats.mean_column_ic);
        }
    }
    if block_ics.len() < 2 {
        return None;
    }
    let mean = block_ics.iter().sum::<f64>() / block_ics.len() as f64;
    let var = block_ics.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / block_ics.len() as f64;
    let std = var.sqrt();
    if mean.abs() > 1e-6 {
        Some((std / mean).min(2.0))
    } else {
        Some(std)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::wrap::Background;

    #[test]
    fn returns_none_when_too_few_rows() {
        let seq = vec![b'A'; 1000];
        let bg = Background::compute(&seq);
        let cfg = DetectorConfig::default();
        // 10 rows < 200 block-rows × 2 → None.
        assert!(compute(&seq, 100, &bg, &cfg).is_none());
    }

    #[test]
    fn uniform_blocks_yield_low_irregularity() {
        let mut seq = Vec::new();
        let block: Vec<u8> = (0..170)
            .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
            .collect();
        for _ in 0..600 {
            seq.extend_from_slice(&block);
        }
        let bg = Background::compute(&seq);
        let cfg = DetectorConfig::default();
        let irr = compute(&seq, 170, &bg, &cfg).expect("expected a score");
        assert!(irr < 0.05, "uniform array → low irregularity; got {irr}");
    }

    #[test]
    fn varying_blocks_yield_higher_irregularity() {
        use rand::Rng;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::seed_from_u64(7);
        let row_len = 170usize;
        let half_rows = 400usize;
        let mut seq = vec![b'A'; row_len * half_rows];
        for _ in 0..row_len * half_rows {
            seq.push(b"ACGT"[rng.random_range(0..4)]);
        }
        let bg = Background::compute(&seq);
        let cfg = DetectorConfig::default();
        let irr = compute(&seq, row_len, &bg, &cfg).expect("expected a score");
        assert!(
            irr > 0.10,
            "half-conserved / half-random → higher irregularity; got {irr}"
        );
    }
}

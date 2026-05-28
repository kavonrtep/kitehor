//! Per-record anchor sampler for `rescore`.
//!
//! For a candidate period `P`, we sample `K` anchor offsets uniformly from
//! `[0, L − 2P − slop]` and pair the tile at the anchor with the tile that
//! follows it (with ±`slop` boundary slack on the partner). Pairs that
//! straddle too much `N` are dropped and re-drawn up to a configurable
//! retry budget; this gives the kernel clean input without paying the
//! N-handling penalty repeatedly.
//!
//! Seeding: ChaCha20 seeded by FNV-1a of `(top_seed, case_id)`, matching
//! the rest of the project (see `src/synth/rng.rs`). Same `(seed, case_id,
//! period)` triple always yields the same anchor sequence, independent of
//! thread count or per-record ordering.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

const FNV_OFFSET: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

/// FNV-1a hash of `top.to_le_bytes() || ":" || case_id`. Mirrors
/// `synth::rng::derive` so behaviour is consistent across stages.
pub fn derive_seed(top: u64, case_id: &str) -> u64 {
    let mut h = FNV_OFFSET;
    for b in top.to_le_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= b':' as u64;
    h = h.wrapping_mul(FNV_PRIME);
    for &b in case_id.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

/// Sampler configuration shared across all (record, period) calls.
#[derive(Debug, Clone, Copy)]
pub struct SampleConfig {
    /// Target number of pairs.
    pub k: usize,
    /// Boundary slack on B in bp. Must satisfy `slop ≤ period`.
    pub slop: usize,
    /// Reject samples with Ns above this fraction of the total pair size.
    pub max_n_frac: f64,
    /// Additional draws per slot if the initial draw is rejected.
    pub max_retries: usize,
    /// Top-level RNG seed (CLI flag).
    pub seed: u64,
}

impl Default for SampleConfig {
    fn default() -> Self {
        Self {
            k: 200,
            slop: 10,
            max_n_frac: 0.05,
            max_retries: 3,
            seed: 42,
        }
    }
}

/// A sampled pair of windows into the array.
///
/// `a_start..a_end` is the anchor tile (length `period`). `b_start..b_end`
/// is the slop-extended adjacent tile (length `period + 2·slop`). Both
/// ranges are half-open and bounded inside the source array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pair {
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
}

/// Sample up to `cfg.k` adjacent-tile pairs for the given (record, period).
///
/// Returns an empty Vec when the array is too short to host any pair, or
/// when `slop > period` (invalid configuration — caller filtered too late).
/// The returned count may be less than `cfg.k` if the retry budget is
/// exhausted by N-heavy windows; callers should propagate the effective
/// count via `identity_n`.
pub fn sample_pairs(seq: &[u8], period: usize, case_id: &str, cfg: &SampleConfig) -> Vec<Pair> {
    if cfg.slop > period {
        return Vec::new();
    }
    let l = seq.len();
    let span = 2 * period + cfg.slop;
    if l < span {
        return Vec::new();
    }
    let max_anchor = l - span; // inclusive upper bound on anchor

    let mut rng = ChaCha20Rng::seed_from_u64(derive_seed(cfg.seed, case_id));
    let mut pairs = Vec::with_capacity(cfg.k);
    let mut attempts_total = 0usize;
    let attempts_cap = cfg.k.saturating_mul(cfg.max_retries.saturating_add(1));

    while pairs.len() < cfg.k && attempts_total < attempts_cap {
        let s = rng.random_range(0..=max_anchor);
        let pair = Pair {
            a_start: s,
            a_end: s + period,
            b_start: s + period - cfg.slop,
            b_end: s + 2 * period + cfg.slop,
        };
        attempts_total += 1;
        if !is_n_heavy(seq, &pair, cfg.max_n_frac) {
            pairs.push(pair);
        }
    }
    pairs
}

fn is_n_heavy(seq: &[u8], pair: &Pair, max_frac: f64) -> bool {
    let a = &seq[pair.a_start..pair.a_end];
    let b = &seq[pair.b_start..pair.b_end];
    let total = a.len() + b.len();
    if total == 0 {
        return false;
    }
    let n_count =
        a.iter().filter(|&&c| c == b'N').count() + b.iter().filter(|&&c| c == b'N').count();
    (n_count as f64) / (total as f64) > max_frac
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(k: usize) -> SampleConfig {
        SampleConfig {
            k,
            slop: 5,
            max_n_frac: 0.05,
            max_retries: 3,
            seed: 1234,
        }
    }

    #[test]
    fn determinism_same_seed_same_pairs() {
        let seq = vec![b'A'; 1000];
        let a = sample_pairs(&seq, 50, "case-x", &cfg(10));
        let b = sample_pairs(&seq, 50, "case-x", &cfg(10));
        assert_eq!(a, b, "same (seed, case_id) must yield same pairs");
        assert_eq!(a.len(), 10);
    }

    #[test]
    fn different_case_ids_diverge() {
        let seq = vec![b'A'; 1000];
        let a = sample_pairs(&seq, 50, "case-1", &cfg(20));
        let b = sample_pairs(&seq, 50, "case-2", &cfg(20));
        assert_ne!(a, b, "different case_ids should yield different anchors");
    }

    #[test]
    fn pairs_stay_within_bounds() {
        let seq = vec![b'A'; 500];
        let pairs = sample_pairs(&seq, 50, "x", &cfg(100));
        for p in &pairs {
            assert!(p.a_start < p.a_end);
            assert!(p.a_end <= seq.len());
            assert!(p.b_start < p.b_end);
            assert!(p.b_end <= seq.len());
            assert_eq!(p.a_end - p.a_start, 50);
            assert_eq!(p.b_end - p.b_start, 50 + 2 * 5);
            // B starts at A's end minus slop (adjacent tile, with slack)
            assert_eq!(p.b_start, p.a_end - 5);
        }
    }

    #[test]
    fn array_too_short_returns_empty() {
        let seq = vec![b'A'; 100];
        // Need 2*P + slop = 200 + 5 = 205; array is 100.
        let pairs = sample_pairs(&seq, 100, "x", &cfg(50));
        assert!(pairs.is_empty());
    }

    #[test]
    fn array_exactly_span_yields_one_anchor() {
        // span = 2*50 + 5 = 105; max_anchor = 0.
        let seq = vec![b'A'; 105];
        let pairs = sample_pairs(&seq, 50, "x", &cfg(20));
        assert!(!pairs.is_empty());
        for p in &pairs {
            assert_eq!(p.a_start, 0);
        }
    }

    #[test]
    fn slop_exceeds_period_returns_empty() {
        let seq = vec![b'A'; 1000];
        let mut c = cfg(10);
        c.slop = 100; // > period=50
        let pairs = sample_pairs(&seq, 50, "x", &c);
        assert!(pairs.is_empty());
    }

    #[test]
    fn n_heavy_pairs_are_rejected_and_retried() {
        // Sequence is mostly N; expect very few pairs to pass the 5% cap.
        let mut seq = vec![b'N'; 1000];
        // Plant a clean stretch big enough for at least one pair.
        seq[100..400].fill(b'A');
        let pairs = sample_pairs(&seq, 50, "x", &cfg(20));
        // We don't require a specific count, but every returned pair must
        // pass the N-cap.
        for p in &pairs {
            let n_a = seq[p.a_start..p.a_end]
                .iter()
                .filter(|&&c| c == b'N')
                .count();
            let n_b = seq[p.b_start..p.b_end]
                .iter()
                .filter(|&&c| c == b'N')
                .count();
            let total = (p.a_end - p.a_start) + (p.b_end - p.b_start);
            assert!((n_a + n_b) as f64 / total as f64 <= 0.05);
        }
    }

    #[test]
    fn all_n_array_returns_empty_after_retries() {
        let seq = vec![b'N'; 1000];
        let pairs = sample_pairs(&seq, 50, "x", &cfg(20));
        assert!(pairs.is_empty(), "all-N array must yield no pairs");
    }

    #[test]
    fn derive_seed_is_deterministic_and_diverges() {
        assert_eq!(derive_seed(42, "case-1"), derive_seed(42, "case-1"));
        assert_ne!(derive_seed(42, "case-1"), derive_seed(42, "case-2"));
        assert_ne!(derive_seed(42, "case-1"), derive_seed(43, "case-1"));
    }
}

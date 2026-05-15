//! Synthetic array generator.
//!
//! Generates one tandem-repeat or HOR array following `ground_truth/
//! simulate_hor.py`'s parameter convention exactly. Supports per-base
//! substitution + indel mutations, post-construction block-level and
//! monomer-level conversion events, and per-case diagnostic metrics
//! (mean intra-block / homologous / cross-position identity).
//!
//! Per-case pipeline (matches Python's `simulate_array`):
//!  1. Build base monomer M of length `monomer_len`. If `submono_k >= 2`,
//!     M is built by tiling a smaller random motif of length
//!     `monomer_len / submono_k` `submono_k` times.
//!  2. Derive `hor_order` founder monomers from M using
//!     `sub_rate_intra` + `indel_rate_intra` (one per founder).
//!  3. Replicate the HOR block `n_blocks` times. Each monomer in each
//!     block copy is mutated with `sub_rate_inter` + `indel_rate_inter`
//!     independently.
//!  4. Apply `block_conversions` random source→target block copies
//!     (strengthens HOR signal). Then `monomer_conversions` random
//!     monomer-level copies (degrades HOR signal).
//!  5. Compute diagnostic identities and return (array, truth, monomers,
//!     events).
//!
//! NB: param naming follows Python:
//!   intra = within-HOR-unit divergence (used when deriving founders)
//!   inter = between-HOR-unit divergence (per-monomer-copy variation)

use crate::errors::Result;
use crate::sequence::ArrayRecord;
use serde::Serialize;
// Levenshtein computed below (hand-rolled Wagner-Fischer); strsim
// 0.11 doesn't give us a clean &[u8] entry point.

#[derive(Debug, Clone)]
pub struct SimulateParams {
    pub monomer_len: usize,
    /// HOR order (multiplicity). 1 means a plain tandem repeat.
    pub hor_order: usize,
    /// Number of HOR copies (n_blocks in truth.tsv).
    pub n_blocks: usize,
    /// Per-base substitution rate when deriving founders from the base.
    /// (Python `sub_rate_intra`.) Applied once per founder.
    pub sub_rate_intra: f64,
    /// Per-base substitution rate per monomer copy.
    /// (Python `sub_rate_inter`.) Applied independently to each copy.
    pub sub_rate_inter: f64,
    /// Per-position indel rate when deriving founders.
    /// (Python `indel_rate_intra`.)
    pub indel_rate_intra: f64,
    /// Per-position indel rate per monomer copy.
    /// (Python `indel_rate_inter`.)
    pub indel_rate_inter: f64,
    /// Number of block-level conversion events. Each picks two block
    /// indices s, t at random and copies block s' monomers onto block t.
    pub block_conversions: usize,
    /// Number of monomer-level conversion events. Each picks two flat
    /// monomer indices s, t and copies monomer s onto t.
    pub monomer_conversions: usize,
    /// If >= 2, each founder is built by tiling a (monomer_len/submono_k)
    /// bp sub-motif `submono_k` times.
    pub submono_k: usize,
    pub seed: u64,
}

impl Default for SimulateParams {
    fn default() -> Self {
        Self {
            monomer_len: 171,
            hor_order: 12,
            n_blocks: 100,
            sub_rate_intra: 0.05,
            sub_rate_inter: 0.03,
            indel_rate_intra: 0.0,
            indel_rate_inter: 0.0,
            block_conversions: 0,
            monomer_conversions: 0,
            submono_k: 1,
            seed: 0,
        }
    }
}

/// Diagnostic record matching truth.tsv's primary columns.
#[derive(Debug, Clone, Serialize)]
pub struct SimulateTruth {
    pub case_id: String,
    pub monomer_len: usize,
    pub hor_order: usize,
    pub n_blocks: usize,
    pub sub_rate_intra: f64,
    pub sub_rate_inter: f64,
    pub indel_rate_intra: f64,
    pub indel_rate_inter: f64,
    pub block_conversions: usize,
    pub monomer_conversions: usize,
    pub submono_k: usize,
    pub seed: u64,
    pub array_length: usize,
    pub n_monomers: usize,
    /// Mean alignment-based identity between monomers in the same block
    /// (NaN when hor_order=1).
    pub mean_intra_block_id: f64,
    /// Mean identity between monomers in different blocks but at the
    /// same founder_idx (NaN when hor_order=1).
    pub mean_homologous_id: f64,
    /// Mean identity between monomers in different blocks AND different
    /// founder_idx.
    pub mean_cross_position_id: f64,
    /// HOR signal: homologous identity minus cross-position identity.
    pub hor_signal: f64,
}

/// Per-monomer record for monomers.tsv.
#[derive(Debug, Clone)]
pub struct MonomerInfo {
    pub block_idx: usize,
    pub founder_idx: usize,
    pub start: usize,
    pub end: usize,
}

/// Conversion event for events.tsv. `scope = "block"` => indices are
/// block indices; `scope = "monomer"` => flat monomer indices.
#[derive(Debug, Clone)]
pub struct ConversionEvent {
    pub event_order: usize,
    pub scope: &'static str,
    pub source_idx: usize,
    pub target_idx: usize,
}

/// Generate one simulated array with full per-case bookkeeping.
pub fn simulate(
    case_id: &str,
    params: &SimulateParams,
) -> Result<(ArrayRecord, SimulateTruth, Vec<MonomerInfo>, Vec<ConversionEvent>)> {
    if params.monomer_len == 0 || params.hor_order == 0 || params.n_blocks == 0 {
        return Err(crate::errors::HordetectError::InvalidParam(
            "monomer_len, hor_order, n_blocks must all be > 0".into(),
        ));
    }
    let submono_k = params.submono_k.max(1);
    if submono_k > 1 && params.monomer_len % submono_k != 0 {
        return Err(crate::errors::HordetectError::InvalidParam(format!(
            "monomer_len ({}) must be divisible by submono_k ({})",
            params.monomer_len, submono_k
        )));
    }

    let mut rng = Xorshift64::new(params.seed.max(1));

    // 1. Build base monomer M.
    let base: Vec<u8> = if submono_k >= 2 {
        let sub_len = params.monomer_len / submono_k;
        let sub_motif = rng.random_dna(sub_len);
        let mut m = Vec::with_capacity(params.monomer_len);
        for _ in 0..submono_k {
            m.extend_from_slice(&sub_motif);
        }
        m
    } else {
        rng.random_dna(params.monomer_len)
    };

    // 2. Build `hor_order` founders from M using INTRA rates.
    //    For hor_order=1, the single founder equals the base (no intra
    //    mutation — matches Python's `simulate_array`).
    let founders: Vec<Vec<u8>> = if params.hor_order == 1 {
        vec![base.clone()]
    } else {
        (0..params.hor_order)
            .map(|_| rng.mutate(&base, params.sub_rate_intra, params.indel_rate_intra))
            .collect()
    };

    // 3. Build all monomers: n_blocks × hor_order, each mutated
    //    independently with INTER rates from its founder.
    let mut monomer_seqs: Vec<Vec<u8>> = Vec::with_capacity(params.n_blocks * params.hor_order);
    let mut monomer_lattice: Vec<(usize, usize)> =
        Vec::with_capacity(params.n_blocks * params.hor_order);
    for b in 0..params.n_blocks {
        for f in 0..params.hor_order {
            let seq = rng.mutate(&founders[f], params.sub_rate_inter, params.indel_rate_inter);
            monomer_seqs.push(seq);
            monomer_lattice.push((b, f));
        }
    }

    // 4a. Block-level conversions.
    let mut events: Vec<ConversionEvent> = Vec::new();
    for _ in 0..params.block_conversions {
        if params.n_blocks < 2 {
            break;
        }
        let (s, t) = rng.sample_two(params.n_blocks);
        for f in 0..params.hor_order {
            let src_idx = s * params.hor_order + f;
            let tgt_idx = t * params.hor_order + f;
            monomer_seqs[tgt_idx] = monomer_seqs[src_idx].clone();
        }
        events.push(ConversionEvent {
            event_order: events.len(),
            scope: "block",
            source_idx: s,
            target_idx: t,
        });
    }
    // 4b. Monomer-level conversions.
    let total_monomers = monomer_seqs.len();
    for _ in 0..params.monomer_conversions {
        if total_monomers < 2 {
            break;
        }
        let (s, t) = rng.sample_two(total_monomers);
        monomer_seqs[t] = monomer_seqs[s].clone();
        events.push(ConversionEvent {
            event_order: events.len(),
            scope: "monomer",
            source_idx: s,
            target_idx: t,
        });
    }

    // 5. Concatenate; build MonomerInfo with start/end coordinates.
    let mut seq: Vec<u8> = Vec::with_capacity(total_monomers * params.monomer_len);
    let mut monomers: Vec<MonomerInfo> = Vec::with_capacity(total_monomers);
    let mut offset = 0usize;
    for (i, ms) in monomer_seqs.iter().enumerate() {
        let (b, f) = monomer_lattice[i];
        let len = ms.len();
        monomers.push(MonomerInfo {
            block_idx: b,
            founder_idx: f,
            start: offset,
            end: offset + len,
        });
        seq.extend_from_slice(ms);
        offset += len;
    }

    let array = ArrayRecord::from_raw(case_id, &seq);

    // 6. Diagnostic identities (alignment-based, sampled).
    let metrics = diagnostic_metrics(
        &monomer_seqs,
        &monomer_lattice,
        params.hor_order,
        &mut rng,
        /* max_pairs_per_category */ 40,
    );

    let truth = SimulateTruth {
        case_id: case_id.to_string(),
        monomer_len: params.monomer_len,
        hor_order: params.hor_order,
        n_blocks: params.n_blocks,
        sub_rate_intra: params.sub_rate_intra,
        sub_rate_inter: params.sub_rate_inter,
        indel_rate_intra: params.indel_rate_intra,
        indel_rate_inter: params.indel_rate_inter,
        block_conversions: params.block_conversions,
        monomer_conversions: params.monomer_conversions,
        submono_k,
        seed: params.seed,
        array_length: array.length,
        n_monomers: total_monomers,
        mean_intra_block_id: metrics.intra,
        mean_homologous_id: metrics.homol,
        mean_cross_position_id: metrics.cross,
        hor_signal: metrics.signal,
    };
    Ok((array, truth, monomers, events))
}

/// Diagnostic identities computed over sampled monomer pairs.
struct PairwiseMetrics {
    intra: f64,
    homol: f64,
    cross: f64,
    signal: f64,
}

fn diagnostic_metrics(
    seqs: &[Vec<u8>],
    lattice: &[(usize, usize)],
    hor_order: usize,
    rng: &mut Xorshift64,
    max_pairs_per_category: usize,
) -> PairwiseMetrics {
    let n = seqs.len();
    let mut intra_pairs: Vec<(usize, usize)> = Vec::new();
    let mut homol_pairs: Vec<(usize, usize)> = Vec::new();
    let mut cross_pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let (bi, fi) = lattice[i];
            let (bj, fj) = lattice[j];
            if bi == bj {
                if hor_order > 1 {
                    intra_pairs.push((i, j));
                }
            } else if fi == fj {
                homol_pairs.push((i, j));
            } else {
                cross_pairs.push((i, j));
            }
        }
    }
    let mean_id = |pairs: Vec<(usize, usize)>, rng: &mut Xorshift64| -> f64 {
        if pairs.is_empty() {
            return f64::NAN;
        }
        let chosen = if pairs.len() > max_pairs_per_category {
            rng.sample_indices(&pairs, max_pairs_per_category)
        } else {
            pairs
        };
        let mut acc = 0.0f64;
        let mut n_kept = 0usize;
        for (i, j) in chosen {
            let id = aligned_identity(&seqs[i], &seqs[j]);
            if !id.is_nan() {
                acc += id;
                n_kept += 1;
            }
        }
        if n_kept == 0 {
            f64::NAN
        } else {
            acc / n_kept as f64
        }
    };
    let intra = mean_id(intra_pairs, rng);
    let homol = mean_id(homol_pairs, rng);
    let cross = mean_id(cross_pairs, rng);
    let signal = if homol.is_nan() || cross.is_nan() {
        f64::NAN
    } else {
        homol - cross
    };
    PairwiseMetrics {
        intra,
        homol,
        cross,
        signal,
    }
}

/// Identity = 1 - levenshtein(a,b) / max(len(a), len(b)). NaN when both empty.
fn aligned_identity(a: &[u8], b: &[u8]) -> f64 {
    let m = a.len().max(b.len());
    if m == 0 {
        return f64::NAN;
    }
    let d = levenshtein_bytes(a, b) as f64;
    1.0 - d / (m as f64)
}

/// Wagner-Fischer Levenshtein over byte slices. O(|a| * |b|) time,
/// O(min(|a|, |b|)) memory.
fn levenshtein_bytes(a: &[u8], b: &[u8]) -> usize {
    if a.is_empty() {
        return b.len();
    }
    if b.is_empty() {
        return a.len();
    }
    // Use the shorter as the inner loop to minimize memory.
    let (a, b) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let mut prev: Vec<usize> = (0..=a.len()).collect();
    let mut curr: Vec<usize> = vec![0; a.len() + 1];
    for (i, &cb) in b.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &ca) in a.iter().enumerate() {
            let cost = if ca == cb { 0 } else { 1 };
            curr[j + 1] = (curr[j] + 1).min(prev[j + 1] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[a.len()]
}

/// Xorshift64 PRNG. Deterministic given (seed, call order).
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }
    fn next(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }
    fn unit(&mut self) -> f64 {
        // 53-bit mantissa to (0, 1); good enough for sampling.
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn random_base(&mut self) -> u8 {
        match self.next() & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }
    fn random_dna(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.random_base()).collect()
    }
    /// Apply substitution + indel mutation per position, matching
    /// Python `simulate_hor.mutate`. Each position independently:
    ///   - with probability `indel_rate`, an indel event happens:
    ///     50/50 either insert a random base before the position
    ///     (and keep the original), or delete (skip the position);
    ///   - else, with probability `sub_rate`, replace the base with
    ///     one of the other three uniformly.
    fn mutate(&mut self, seq: &[u8], sub_rate: f64, indel_rate: f64) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::with_capacity(seq.len() + (seq.len() / 20));
        for &base in seq {
            if indel_rate > 0.0 && self.unit() < indel_rate {
                if self.unit() < 0.5 {
                    out.push(self.random_base());
                    out.push(base);
                }
                // else: deletion (skip this position)
            } else if sub_rate > 0.0 && self.unit() < sub_rate {
                out.push(alt_base(base, self.next()));
            } else {
                out.push(base);
            }
        }
        out
    }
    /// Sample two distinct indices from [0, n).
    fn sample_two(&mut self, n: usize) -> (usize, usize) {
        debug_assert!(n >= 2);
        let a = (self.next() % n as u64) as usize;
        loop {
            let b = (self.next() % n as u64) as usize;
            if b != a {
                return (a, b);
            }
        }
    }
    /// Random reservoir-style sample of `k` items from `pairs`. Order
    /// of the output is not preserved (we shuffle in place).
    fn sample_indices(
        &mut self,
        pairs: &[(usize, usize)],
        k: usize,
    ) -> Vec<(usize, usize)> {
        let n = pairs.len();
        if k >= n {
            return pairs.to_vec();
        }
        // Fisher-Yates partial.
        let mut idx: Vec<usize> = (0..n).collect();
        for i in 0..k {
            let j = i + (self.next() as usize % (n - i));
            idx.swap(i, j);
        }
        idx[..k].iter().map(|&i| pairs[i]).collect()
    }
}

fn alt_base(b: u8, r: u64) -> u8 {
    let alts: &[u8] = match b {
        b'A' | b'a' => &[b'C', b'G', b'T'],
        b'C' | b'c' => &[b'A', b'G', b'T'],
        b'G' | b'g' => &[b'A', b'C', b'T'],
        b'T' | b't' => &[b'A', b'C', b'G'],
        _ => &[b'A', b'C', b'G'],
    };
    alts[(r as usize) % 3]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simulate_basic_tandem_repeat() {
        let params = SimulateParams {
            monomer_len: 200,
            hor_order: 1,
            n_blocks: 20,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            submono_k: 1,
            seed: 42,
            ..SimulateParams::default()
        };
        let (arr, truth, monomers, events) = simulate("test", &params).unwrap();
        assert_eq!(arr.length, 200 * 20);
        assert_eq!(truth.n_monomers, 20);
        assert_eq!(monomers.len(), 20);
        assert!(events.is_empty());
        for i in 1..20 {
            assert_eq!(&arr.seq[0..200], &arr.seq[i * 200..(i + 1) * 200]);
        }
    }

    #[test]
    fn simulate_hor_intra_makes_distinct_founders() {
        // Python's `sub_intra` makes founders differ. With hor_order=4
        // and intra=0.10, founder 0 vs founder 1 should differ ~10%.
        let params = SimulateParams {
            monomer_len: 200,
            hor_order: 4,
            n_blocks: 10,
            sub_rate_intra: 0.10,
            sub_rate_inter: 0.0,
            submono_k: 1,
            seed: 42,
            ..SimulateParams::default()
        };
        let (arr, _, _, _) = simulate("test", &params).unwrap();
        let f0 = &arr.seq[0..200];
        let f1 = &arr.seq[200..400];
        let diff = f0.iter().zip(f1.iter()).filter(|(a, b)| a != b).count();
        assert!(diff > 10 && diff < 40, "expected ~10% divergence, got {diff}/200");
        // Same founder 0 in block 2 should still equal f0 (no inter).
        let f0_block2 = &arr.seq[800..1000];
        assert_eq!(f0, f0_block2);
    }

    #[test]
    fn simulate_submono_tiles_internally() {
        let params = SimulateParams {
            monomer_len: 180,
            hor_order: 1,
            n_blocks: 10,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            submono_k: 3,
            seed: 42,
            ..SimulateParams::default()
        };
        let (arr, _, _, _) = simulate("sm", &params).unwrap();
        assert_eq!(arr.length, 180 * 10);
        let tile_a = &arr.seq[0..60];
        let tile_b = &arr.seq[60..120];
        let tile_c = &arr.seq[120..180];
        assert_eq!(tile_a, tile_b);
        assert_eq!(tile_b, tile_c);
    }

    #[test]
    fn simulate_inter_divergence_applies_per_copy() {
        // Python's `sub_inter` adds per-copy variation between blocks.
        let params = SimulateParams {
            monomer_len: 200,
            hor_order: 1,
            n_blocks: 20,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.10,
            submono_k: 1,
            seed: 42,
            ..SimulateParams::default()
        };
        let (arr, _, _, _) = simulate("d", &params).unwrap();
        let c0 = &arr.seq[0..200];
        let c1 = &arr.seq[200..400];
        let diff = c0.iter().zip(c1.iter()).filter(|(a, b)| a != b).count();
        // Each copy diverged from the founder by ~10% independently.
        // Expected pairwise diff ~ 2*0.10*(1-1/4) ≈ 15%.
        assert!(diff > 10 && diff < 60, "expected pairwise diff 10-60%, got {diff}/200");
    }

    #[test]
    fn simulate_block_conversion_strengthens_signal() {
        let mut params = SimulateParams {
            monomer_len: 100,
            hor_order: 3,
            n_blocks: 6,
            sub_rate_intra: 0.10,
            sub_rate_inter: 0.05,
            submono_k: 1,
            seed: 1,
            ..SimulateParams::default()
        };
        // Zero conversions.
        params.block_conversions = 0;
        let (_, t0, _, _) = simulate("c0", &params).unwrap();
        // Many block conversions: forces blocks to look identical.
        params.block_conversions = 50;
        let (_, t50, _, e50) = simulate("c50", &params).unwrap();
        assert!(!e50.is_empty(), "expected some block events");
        // Homologous identity should INCREASE with block conversions.
        assert!(
            t50.mean_homologous_id > t0.mean_homologous_id,
            "block conv should raise homol id ({} vs {})",
            t50.mean_homologous_id,
            t0.mean_homologous_id
        );
    }

    #[test]
    fn simulate_monomer_conversion_degrades_hor_signal() {
        let mut params = SimulateParams {
            monomer_len: 100,
            hor_order: 4,
            n_blocks: 8,
            sub_rate_intra: 0.10,
            sub_rate_inter: 0.03,
            submono_k: 1,
            seed: 7,
            ..SimulateParams::default()
        };
        let (_, t0, _, _) = simulate("m0", &params).unwrap();
        params.monomer_conversions = 100;
        let (_, t_hi, _, e) = simulate("m100", &params).unwrap();
        assert!(!e.is_empty());
        // Monomer conversions blur founder distinctions → lower hor_signal.
        assert!(
            t_hi.hor_signal < t0.hor_signal,
            "expected hor_signal to drop; got {} vs {}",
            t_hi.hor_signal,
            t0.hor_signal
        );
    }

    #[test]
    fn simulate_indel_changes_lengths() {
        let params = SimulateParams {
            monomer_len: 200,
            hor_order: 1,
            n_blocks: 10,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            indel_rate_intra: 0.0,
            indel_rate_inter: 0.05,
            submono_k: 1,
            seed: 3,
            ..SimulateParams::default()
        };
        let (arr, truth, monomers, _) = simulate("i", &params).unwrap();
        // With 5% indels each position 50/50 +1/-1, expected length
        // diff vs original = 0. Array length should still be near 2000
        // but not exactly.
        assert!(arr.length > 1500 && arr.length < 2500);
        assert_eq!(truth.n_monomers, 10);
        // Monomer lengths vary.
        let lengths: Vec<usize> = monomers.iter().map(|m| m.end - m.start).collect();
        let unique: std::collections::HashSet<_> = lengths.iter().collect();
        assert!(unique.len() > 1, "expected variable monomer lengths under indels");
    }

    #[test]
    fn simulate_is_deterministic_for_fixed_seed() {
        let params = SimulateParams {
            monomer_len: 200,
            hor_order: 4,
            n_blocks: 8,
            sub_rate_intra: 0.05,
            sub_rate_inter: 0.05,
            indel_rate_intra: 0.005,
            indel_rate_inter: 0.005,
            block_conversions: 2,
            monomer_conversions: 3,
            submono_k: 1,
            seed: 12345,
        };
        let (a, _, _, _) = simulate("x", &params).unwrap();
        let (b, _, _, _) = simulate("x", &params).unwrap();
        assert_eq!(a.seq, b.seq);
    }

    #[test]
    fn simulate_rejects_bad_submono_divisor() {
        let params = SimulateParams {
            monomer_len: 100,
            hor_order: 1,
            n_blocks: 5,
            sub_rate_intra: 0.0,
            sub_rate_inter: 0.0,
            submono_k: 3,
            seed: 1,
            ..SimulateParams::default()
        };
        assert!(simulate("x", &params).is_err());
    }

}

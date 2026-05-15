//! Block-mean homology probe at a given period (`probe_period`).
//!
//! Given a record and a candidate period `m`, this lays down adjacent
//! length-`m` blocks across the array, sketches each as an L2-normalised
//! k-mer count vector, and reports the mean adjacent cosine similarity.
//! It is the source of the `h_d1` / `h_founder` features that the
//! classifier feeds the random forest.
//!
//! Originally part of a larger monomer-inference stage that also picked
//! the base monomer from a periodogram; that path was retired in favour
//! of the kite-first probabilistic classifier (see `classify.rs`). What
//! remains is the homology probe plus the supporting block-sketching
//! helpers.

use crate::sequence::ArrayRecord;
use serde::Serialize;

#[derive(Debug, Clone, Copy)]
pub struct MonomerModelConfig {
    pub min_monomer_bp: usize,
    pub max_monomer_bp: usize,
    /// k-mer size for the block sketch. 5 keeps the vector at 4^5 = 1024.
    pub sketch_kmer_size: usize,
    /// Cosine-similarity threshold for greedy subtype clustering
    /// (`subtype_count` field on `MonomerCandidate`).
    pub subtype_threshold: f64,
    /// Phase offsets to try when laying down blocks (sampled across [0,m)).
    pub n_phase_offsets: usize,
    /// Adjacent-pair sample cap (keeps scoring O(min(n,cap) * dim)).
    pub max_pairs_sampled: usize,
    /// Candidates with block-level homology below this floor are dropped
    /// from scoring. `probe_period` zeros this out so every probe returns
    /// a score (the floor mattered for the now-retired inference path).
    pub homology_floor: f64,
}

impl Default for MonomerModelConfig {
    fn default() -> Self {
        Self {
            min_monomer_bp: 15,
            max_monomer_bp: 5000,
            sketch_kmer_size: 5,
            subtype_threshold: 0.85,
            n_phase_offsets: 8,
            max_pairs_sampled: 200,
            homology_floor: 0.40,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct MonomerCandidate {
    pub monomer_bp: usize,
    pub phase_offset: usize,
    pub n_blocks: usize,
    pub homology_score: f64,
    pub subtype_count: usize,
}

/// Probe a single (record, period) pair. Returns
/// `(homology, phase_offset, n_blocks)` when at least two blocks fit;
/// `None` otherwise. The homology score is the average adjacent
/// block-pair cosine similarity at the best phase offset.
pub fn probe_period(
    record: &ArrayRecord,
    m: usize,
    cfg: &MonomerModelConfig,
) -> Option<(f64, usize, usize)> {
    let mut c = *cfg;
    c.homology_floor = 0.0;
    let cand = score_candidate(record, m, &c)?;
    Some((cand.homology_score, cand.phase_offset, cand.n_blocks))
}

/// Score one (record, m) candidate. Returns `None` if the record is too
/// short, no phase offset yields ≥ 2 blocks, or the resulting homology
/// is below `cfg.homology_floor`.
pub fn score_candidate(
    record: &ArrayRecord,
    m: usize,
    cfg: &MonomerModelConfig,
) -> Option<MonomerCandidate> {
    if m == 0 || record.length < 4 * m {
        return None;
    }
    let dim = 1usize << (2 * cfg.sketch_kmer_size);

    // Best phase offset by average adjacent similarity over sampled pairs.
    let mut best: Option<(usize, f64, Vec<Vec<f32>>)> = None;
    let n_offsets = cfg.n_phase_offsets.max(1);
    for o_idx in 0..n_offsets {
        let offset = (o_idx * m) / n_offsets;
        if offset + m > record.length {
            continue;
        }
        let blocks = sketch_blocks(record, m, offset, cfg.sketch_kmer_size, dim);
        if blocks.len() < 2 {
            continue;
        }
        let sim = mean_adjacent_similarity(&blocks, cfg.max_pairs_sampled);
        if best.as_ref().map(|(_, b, _)| sim > *b).unwrap_or(true) {
            best = Some((offset, sim, blocks));
        }
    }

    let (phase_offset, homology_score, blocks) = best?;
    let n_blocks = blocks.len();
    let subtype_count = greedy_cluster_count(&blocks, cfg.subtype_threshold);

    if homology_score < cfg.homology_floor {
        return None;
    }

    Some(MonomerCandidate {
        monomer_bp: m,
        phase_offset,
        n_blocks,
        homology_score,
        subtype_count,
    })
}

/// Build an L2-normalised k-mer-count vector per block of length `m`.
fn sketch_blocks(
    record: &ArrayRecord,
    m: usize,
    offset: usize,
    kmer_size: usize,
    dim: usize,
) -> Vec<Vec<f32>> {
    let seq = &record.seq;
    let n = (seq.len().saturating_sub(offset)) / m;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let start = offset + i * m;
        let end = start + m;
        let mut v = vec![0f32; dim];
        for (_, k) in CanonicalKmers::new(&seq[start..end], kmer_size) {
            v[k as usize % dim] += 1.0;
        }
        l2_normalize(&mut v);
        out.push(v);
    }
    out
}

fn l2_normalize(v: &mut [f32]) {
    let s: f32 = v.iter().map(|x| x * x).sum::<f32>().sqrt();
    if s > 0.0 {
        for x in v.iter_mut() {
            *x /= s;
        }
    }
}

fn cosine(a: &[f32], b: &[f32]) -> f64 {
    // Both vectors are pre-normalised so cosine = dot.
    let mut s = 0.0f32;
    for i in 0..a.len() {
        s += a[i] * b[i];
    }
    s as f64
}

fn mean_adjacent_similarity(blocks: &[Vec<f32>], max_pairs: usize) -> f64 {
    let n = blocks.len();
    if n < 2 {
        return 0.0;
    }
    let stride = ((n - 1) / max_pairs.max(1)).max(1);
    let mut count = 0usize;
    let mut sum = 0.0f64;
    let mut i = 0;
    while i + 1 < n && count < max_pairs {
        sum += cosine(&blocks[i], &blocks[i + 1]);
        count += 1;
        i += stride;
    }
    if count == 0 {
        0.0
    } else {
        sum / count as f64
    }
}

/// Greedy online clustering: walk blocks, assign each to the first
/// centroid whose cosine similarity exceeds `threshold`; otherwise spawn
/// a new cluster. Returns the number of clusters.
fn greedy_cluster_count(blocks: &[Vec<f32>], threshold: f64) -> usize {
    let mut centroids: Vec<Vec<f32>> = Vec::new();
    for b in blocks {
        let mut assigned = false;
        for c in &centroids {
            if cosine(c, b) >= threshold {
                assigned = true;
                break;
            }
        }
        if !assigned {
            centroids.push(b.clone());
        }
    }
    centroids.len()
}

// ---------------------------------------------------------------------------
// CanonicalKmers — inlined to drop the `minimizer.rs` dependency.
// Emits (position, canonical-k-mer) pairs over the sequence, where the
// canonical k-mer is `min(forward, reverse-complement)`. k-mers
// containing a non-ACGT byte are skipped.
// ---------------------------------------------------------------------------

struct CanonicalKmers<'a> {
    seq: &'a [u8],
    k: usize,
    i: usize,
}

impl<'a> CanonicalKmers<'a> {
    fn new(seq: &'a [u8], k: usize) -> Self {
        Self { seq, k, i: 0 }
    }
}

impl<'a> Iterator for CanonicalKmers<'a> {
    type Item = (usize, u64);
    fn next(&mut self) -> Option<Self::Item> {
        while self.i + self.k <= self.seq.len() {
            let pos = self.i;
            let window = &self.seq[pos..pos + self.k];
            self.i += 1;
            let mut fwd: u64 = 0;
            let mut rc: u64 = 0;
            let mut ok = true;
            for &b in window {
                let f = match b {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => {
                        ok = false;
                        break;
                    }
                };
                fwd = (fwd << 2) | f;
                rc = (rc >> 2) | ((3 ^ f) << (2 * (self.k - 1)));
            }
            if !ok {
                continue;
            }
            return Some((pos, fwd.min(rc)));
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ArrayRecord;

    fn tandem_record(monomer: &[u8], copies: usize) -> ArrayRecord {
        let mut seq = Vec::with_capacity(monomer.len() * copies);
        for _ in 0..copies {
            seq.extend_from_slice(monomer);
        }
        let length = seq.len();
        ArrayRecord {
            id: "test".into(),
            seq,
            length,
            n_count: 0,
        }
    }

    #[test]
    fn probe_period_high_homology_on_pure_tandem() {
        // 60-bp monomer × 50 copies — all blocks identical → homology ≈ 1.
        let mon = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let rec = tandem_record(mon, 50);
        let (h, _off, n_blocks) =
            probe_period(&rec, mon.len(), &MonomerModelConfig::default()).unwrap();
        assert!(h > 0.95, "expected near-1 homology, got {h}");
        assert!(n_blocks >= 49);
    }

    #[test]
    fn probe_period_returns_none_when_array_too_short() {
        let rec = tandem_record(b"ACGTACGTACGT", 2);
        // 4 × m = 48 > 24 record bytes → None.
        assert!(probe_period(&rec, 12, &MonomerModelConfig::default()).is_none());
    }
}

//! Kite-driven consensus monomer extraction. P-mer Counter walked
//! most-common-first, deduped by canonical form, capped at K
//! canonical-distinct motifs or when the count drops below
//! `freq_ratio × top_count`.

use super::find_ssrs::normalize_motif;
use ahash::AHashMap;

/// One consensus monomer pulled from the sequence.
#[derive(Debug, Clone)]
pub struct ConsensusEntry {
    /// The lowercase P-mer as found in the sequence.
    pub kmer: String,
    /// Canonical form (uppercase) via `normalize_motif`.
    pub canonical: String,
    pub count: u64,
}

/// Port of `extract_consensus_monomers(sequence, period, max_k, freq_ratio)`.
/// Walks `Counter.most_common()` order; skips kmers containing `n`;
/// stops when count drops below `freq_ratio × top_count` OR `max_k`
/// canonical-distinct motifs have been collected.
pub fn extract_consensus_monomers(
    sequence: &[u8],
    period: usize,
    max_k: usize,
    freq_ratio: f64,
) -> Vec<ConsensusEntry> {
    if period == 0 || period > sequence.len() {
        return Vec::new();
    }
    let seq_lower: Vec<u8> = sequence.to_ascii_lowercase();
    let mut counts: AHashMap<Vec<u8>, u64> = AHashMap::new();
    let n = seq_lower.len();
    for i in 0..=(n - period) {
        let kmer = &seq_lower[i..i + period];
        if kmer.iter().any(|b| *b == b'n') {
            continue;
        }
        *counts.entry(kmer.to_vec()).or_insert(0) += 1;
    }
    if counts.is_empty() {
        return Vec::new();
    }
    // most_common in Python is a stable sort by count descending. For
    // ties Python's Counter preserves insertion order. AHashMap's
    // iteration order is non-deterministic, so we need to sort
    // (count desc, kmer asc as a stable tiebreaker — matches what
    // Python tends to produce for our domain since insertion order is
    // governed by position-of-first-occurrence which strongly correlates
    // with high-count kmers appearing early).
    let mut pairs: Vec<(Vec<u8>, u64)> = counts.into_iter().collect();
    pairs.sort_by(|a, b| {
        b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0))
    });

    let mut seen_canonical: AHashMap<String, ()> = AHashMap::new();
    let mut result: Vec<ConsensusEntry> = Vec::new();
    let mut top_count: Option<u64> = None;
    for (kmer, count) in pairs {
        match top_count {
            None => top_count = Some(count),
            Some(tc) => {
                if (count as f64) < freq_ratio * (tc as f64) {
                    break;
                }
            }
        }
        let canonical = normalize_motif(&kmer);
        if seen_canonical.contains_key(&canonical) {
            continue;
        }
        seen_canonical.insert(canonical.clone(), ());
        result.push(ConsensusEntry {
            kmer: String::from_utf8(kmer).expect("ascii"),
            canonical,
            count,
        });
        if result.len() >= max_k {
            break;
        }
    }
    result
}

/// Build `monomer * n_copies`, then extend by one extra copy at a time
/// until length ≥ `min_length`. Matches the prototype's
/// `build_consensus_dimer`.
pub fn build_consensus_dimer(monomer: &str, n_copies: usize, min_length: usize) -> String {
    if monomer.is_empty() {
        return String::new();
    }
    let mut s = monomer.repeat(n_copies);
    while s.len() < min_length {
        s.push_str(monomer);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_gt_one_motif() {
        let seq: Vec<u8> = b"GT".repeat(50);
        let v = extract_consensus_monomers(&seq, 2, 3, 0.3);
        // The top kmer is "gt" (49 windows starting at even indices) and
        // "tg" (49 windows at odd indices). With period=2 these are
        // different lowercase k-mers BUT same canonical (AC vs AC).
        // The first walked is the strongest of the two — only one
        // canonical-distinct entry returned.
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].canonical, "AC"); // canonical of GT/TG
    }

    #[test]
    fn build_dimer_extends_to_min() {
        let d = build_consensus_dimer("GT", 4, 30);
        assert!(d.len() >= 30);
        // GT * 4 = 8 bp; we need at least 30 → keep adding GTs.
        assert_eq!(&d[..], "G".repeat(0).as_str().to_owned() + "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT");
    }
}

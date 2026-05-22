//! TideCluster `find_ssrs` port. The Python prototype uses the regex
//! `(([gatc]{L})\2{min-1,})` per motif length. Rust's `regex` crate
//! does not support backreferences, so we hand-roll the equivalent
//! greedy non-overlapping scanner — exact semantics match.

use super::MotifSpec;
use ahash::AHashSet;

/// One raw SSR hit (lowercase motif, 1-based-inclusive start, 0-based-
/// exclusive end — matches TideCluster convention verbatim).
#[derive(Debug, Clone)]
pub struct Hit {
    pub ssr_number: u32,
    pub motif_length: usize,
    /// Lowercase, as found in the input. The canonical form lives in
    /// `normalized_motif`.
    pub motif_sequence: String,
    pub repeats: usize,
    pub start: usize,
    pub end: usize,
    /// Canonical form (uppercase): lex-min over rotations of `motif`
    /// AND `reverse_complement(motif)`.
    pub normalized_motif: String,
}

/// Port of `tools/rule_proto/ssr_scan.py::find_ssrs`.
///
/// **Convention**: `start` is 1-based inclusive; `end` is 0-based
/// exclusive (the TideCluster idiom — preserve verbatim).
///
/// Semantics: per spec `(L, min_reps)`, scan left-to-right; at each
/// position attempt to grow a greedy run of L-base motif repeats; emit
/// a hit when the run reaches `min_reps` complete repeats AND the
/// start hasn't been claimed by a hit at a shorter motif length.
pub fn find_ssrs(sequence: &[u8], specs: &[MotifSpec]) -> Vec<Hit> {
    let mut seq_lower = sequence.to_ascii_lowercase();
    // Anything other than a/c/g/t becomes 'n' so the motif-character
    // filter below can simply reject `n` (matches the prototype's
    // `[gatc]` character class behavior on `seq.lower()`).
    for b in seq_lower.iter_mut() {
        if !matches!(*b, b'a' | b'c' | b'g' | b't') {
            *b = b'n';
        }
    }
    let len = seq_lower.len();

    let mut results: Vec<Hit> = Vec::new();
    let mut homopolymer_buf: Vec<Hit> = Vec::new();
    let mut locations: AHashSet<usize> = AHashSet::new();

    for spec in specs {
        let l = spec.motif_length;
        let min_reps = spec.min_repeats;
        if l == 0 || min_reps == 0 || l > len {
            continue;
        }
        let mut i = 0usize;
        while i + l * min_reps <= len {
            let motif = &seq_lower[i..i + l];
            // motif characters must all be valid ACGT (no 'n').
            if motif.iter().any(|b| !matches!(*b, b'a' | b'c' | b'g' | b't')) {
                i += 1;
                continue;
            }
            // Greedy run length: count contiguous L-base motif copies.
            let mut j = i + l;
            while j + l <= len && &seq_lower[j..j + l] == motif {
                j += l;
            }
            let n_reps = (j - i) / l;
            if n_reps >= min_reps {
                let start = i + 1; // 1-based inclusive
                let end = j; // 0-based exclusive
                let motif_lower = std::str::from_utf8(motif)
                    .expect("ascii motif")
                    .to_string();
                let is_homo = is_homopolymer(&motif_lower);
                if is_homo {
                    if l == 1 {
                        let canonical = normalize_motif(motif);
                        homopolymer_buf.push(Hit {
                            ssr_number: 0, // assigned later
                            motif_length: l,
                            motif_sequence: motif_lower,
                            repeats: n_reps,
                            start,
                            end,
                            normalized_motif: canonical,
                        });
                    }
                } else if locations.insert(start) {
                    let canonical = normalize_motif(motif);
                    results.push(Hit {
                        ssr_number: (results.len() + 1) as u32,
                        motif_length: l,
                        motif_sequence: motif_lower,
                        repeats: n_reps,
                        start,
                        end,
                        normalized_motif: canonical,
                    });
                }
                // Non-overlapping: skip past the entire run.
                i = j;
            } else {
                i += 1;
            }
        }
    }

    if results.is_empty() && !homopolymer_buf.is_empty() {
        for (idx, mut h) in homopolymer_buf.into_iter().enumerate() {
            h.ssr_number = (idx + 1) as u32;
            results.push(h);
        }
    }
    results
}

/// True iff `motif` is a single nucleotide repeated. Matches the
/// prototype's `homopolymer()` (motif_length < 2 → True; else
/// regex `([gatc])\1{L-1}`). For lowercase ASCII motifs only.
pub fn is_homopolymer(motif: &str) -> bool {
    if motif.len() < 2 {
        return true;
    }
    let bytes = motif.as_bytes();
    let first = bytes[0];
    bytes.iter().all(|&b| b == first)
}

/// Reverse complement of an ASCII nucleotide string. Preserves case.
pub fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'N' => b'N',
            b'n' => b'n',
            other => other,
        })
        .collect()
}

/// Canonical form: lex-min over all rotations of `motif.upper()` AND
/// all rotations of `reverse_complement(motif.upper())`. Matches the
/// prototype's `normalize_motif` byte-for-byte.
pub fn normalize_motif(motif: &[u8]) -> String {
    let upper: Vec<u8> = motif.iter().map(|b| b.to_ascii_uppercase()).collect();
    if upper.is_empty() {
        return String::new();
    }
    let l = upper.len();
    let doubled: Vec<u8> = upper.iter().chain(upper.iter()).copied().collect();
    let rc = reverse_complement(&upper);
    let rc_doubled: Vec<u8> = rc.iter().chain(rc.iter()).copied().collect();
    let mut best: Vec<u8> = doubled[0..l].to_vec();
    for i in 1..l {
        let slice = &doubled[i..i + l];
        if slice < best.as_slice() {
            best = slice.to_vec();
        }
    }
    for i in 0..l {
        let slice = &rc_doubled[i..i + l];
        if slice < best.as_slice() {
            best = slice.to_vec();
        }
    }
    String::from_utf8(best).expect("ascii")
}

/// Drop motifs that are an exact integer-repeat of a shorter motif in
/// the set (e.g. ATAT shrinks to AT).
pub fn get_unique_motifs(motifs: &[String]) -> Vec<String> {
    let mut sorted: Vec<String> = motifs.to_vec();
    sorted.sort_by_key(|m| m.len());
    let mut out: Vec<String> = Vec::new();
    for m in sorted {
        let mut is_unique = true;
        for s in &out {
            if !s.is_empty() && m.len() % s.len() == 0 {
                let times = m.len() / s.len();
                let multiplied: String = s.repeat(times);
                if multiplied == m {
                    is_unique = false;
                    break;
                }
            }
        }
        if is_unique {
            out.push(m);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAT"), b"ATTT");
        assert_eq!(reverse_complement(b"GT"), b"AC");
    }

    #[test]
    fn normalize_motif_examples() {
        assert_eq!(normalize_motif(b"ta"), "AT");
        assert_eq!(normalize_motif(b"gt"), "AC");
        assert_eq!(normalize_motif(b"ct"), "AG");
        assert_eq!(normalize_motif(b"AT"), "AT");
        assert_eq!(normalize_motif(b"AAAT"), "AAAT");
    }

    #[test]
    fn homopolymer_detection() {
        assert!(is_homopolymer("a"));
        assert!(is_homopolymer("aa"));
        assert!(is_homopolymer("aaaaa"));
        assert!(!is_homopolymer("at"));
        assert!(!is_homopolymer("aat"));
    }

    #[test]
    fn drop_multi_motifs() {
        let m = ["AT".to_string(), "ATAT".to_string(), "AAT".to_string()];
        let out = get_unique_motifs(&m);
        assert!(out.iter().any(|s| s == "AT"));
        assert!(out.iter().any(|s| s == "AAT"));
        assert!(!out.iter().any(|s| s == "ATAT"));
    }

    #[test]
    fn find_ssrs_basic_at_dinuc() {
        let seq: Vec<u8> = b"AT".repeat(20);
        let cfg = super::super::Config::default();
        let hits = find_ssrs(&seq, &cfg.specs);
        assert!(!hits.is_empty());
        let h = &hits[0];
        assert_eq!(h.motif_length, 2);
        assert_eq!(h.motif_sequence, "at");
        assert_eq!(h.repeats, 20);
        assert_eq!(h.start, 1);
        assert_eq!(h.end, 40);
        assert_eq!(h.normalized_motif, "AT");
    }

    #[test]
    fn find_ssrs_homopolymer_only_emits_when_alone() {
        let seq = b"AAAAAAAAAAAAAAAAAAAA".to_vec(); // 20 A's
        let cfg = super::super::Config::default();
        let hits = find_ssrs(&seq, &cfg.specs);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].motif_length, 1);
    }
}

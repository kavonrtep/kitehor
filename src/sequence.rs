//! Sequence container and base normalization.
//!
//! Bases are normalized to ACGT (uppercase) or N. We store the sequence as a
//! plain `Vec<u8>` of ASCII bytes — every downstream stage treats unrecognised
//! input as N and either skips or down-weights it.

/// Per-record container used throughout the pipeline.
#[derive(Debug, Clone)]
pub struct ArrayRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub length: usize,
    pub n_count: usize,
}

impl ArrayRecord {
    /// Build an `ArrayRecord` from a raw ID and ASCII byte slice. Bases are
    /// normalized in-place: uppercase ACGT pass through, everything else
    /// becomes `b'N'`. The N count is tallied for downstream QC.
    pub fn from_raw(id: impl Into<String>, raw: &[u8]) -> Self {
        let mut seq = Vec::with_capacity(raw.len());
        let mut n_count = 0usize;
        for &b in raw {
            let nb = normalize_base(b);
            if nb == b'N' {
                n_count += 1;
            }
            seq.push(nb);
        }
        let length = seq.len();
        Self {
            id: id.into(),
            seq,
            length,
            n_count,
        }
    }

    pub fn n_fraction(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.n_count as f64 / self.length as f64
        }
    }
}

/// Map any input byte to one of {A, C, G, T, N}.
#[inline]
pub fn normalize_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' | b'U' | b'u' => b'T',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_basics() {
        assert_eq!(normalize_base(b'A'), b'A');
        assert_eq!(normalize_base(b'a'), b'A');
        assert_eq!(normalize_base(b'T'), b'T');
        assert_eq!(normalize_base(b'u'), b'T');
        assert_eq!(normalize_base(b'N'), b'N');
        assert_eq!(normalize_base(b'X'), b'N');
        assert_eq!(normalize_base(b'-'), b'N');
        assert_eq!(normalize_base(b'\n'), b'N');
    }

    #[test]
    fn record_normalizes_and_counts_ns() {
        // Each pair is (upper, lower) of the same base, then two pairs of N-likes.
        let r = ArrayRecord::from_raw("test", b"AaCcGgTtNnXx");
        assert_eq!(r.id, "test");
        assert_eq!(r.length, 12);
        assert_eq!(&r.seq[..], b"AACCGGTTNNNN");
        assert_eq!(r.n_count, 4);
        assert!((r.n_fraction() - 4.0 / 12.0).abs() < 1e-9);
    }

    #[test]
    fn empty_record_is_zero_n_fraction() {
        let r = ArrayRecord::from_raw("empty", b"");
        assert_eq!(r.length, 0);
        assert_eq!(r.n_fraction(), 0.0);
    }

    #[test]
    fn rna_u_becomes_t() {
        let r = ArrayRecord::from_raw("rna", b"AUGCU");
        assert_eq!(&r.seq[..], b"ATGCT");
        assert_eq!(r.n_count, 0);
    }
}

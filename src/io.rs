//! FASTA input via `needletail`, plus per-record QC gates.

use crate::errors::{HordetectError, Result};
use crate::sequence::ArrayRecord;
use needletail::parse_fastx_file;
use std::path::Path;

/// QC thresholds applied at load time. Records that fail are reported but
/// do NOT abort the run; callers decide whether to skip or pass through.
#[derive(Debug, Clone, Copy)]
pub struct LoadQc {
    pub min_array_bp: usize,
    pub max_n_fraction: f64,
}

impl Default for LoadQc {
    fn default() -> Self {
        Self {
            min_array_bp: 5_000,
            max_n_fraction: 0.20,
        }
    }
}

/// Status attached to every loaded record so the caller can route it.
#[derive(Debug, Clone)]
pub enum LoadStatus {
    Ok,
    TooShort { length: usize, min: usize },
    TooManyNs { n_fraction: f64, limit: f64 },
}

#[derive(Debug, Clone)]
pub struct LoadedRecord {
    pub record: ArrayRecord,
    pub status: LoadStatus,
}

/// Stream a FASTA file, normalizing bases and applying QC. Returns all
/// records (even failing ones) so the summary writer can report failures.
pub fn load_fasta<P: AsRef<Path>>(path: P, qc: LoadQc) -> Result<Vec<LoadedRecord>> {
    let mut reader = parse_fastx_file(path.as_ref())
        .map_err(|e| HordetectError::Fasta(format!("opening {:?}: {}", path.as_ref(), e)))?;

    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| HordetectError::Fasta(format!("record: {e}")))?;
        let id = std::str::from_utf8(rec.id())
            .map_err(|e| HordetectError::Fasta(format!("non-utf8 id: {e}")))?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = rec.seq();
        let arr = ArrayRecord::from_raw(id, &seq);

        let status = if arr.length < qc.min_array_bp {
            LoadStatus::TooShort {
                length: arr.length,
                min: qc.min_array_bp,
            }
        } else if arr.n_fraction() > qc.max_n_fraction {
            LoadStatus::TooManyNs {
                n_fraction: arr.n_fraction(),
                limit: qc.max_n_fraction,
            }
        } else {
            LoadStatus::Ok
        };

        out.push(LoadedRecord {
            record: arr,
            status,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn fasta(records: &[(&str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        for (id, seq) in records {
            writeln!(f, ">{}", id).unwrap();
            // Wrap at 60 bp to exercise multi-line parsing.
            for chunk in seq.as_bytes().chunks(60) {
                f.write_all(chunk).unwrap();
                f.write_all(b"\n").unwrap();
            }
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn loads_single_record() {
        let f = fasta(&[("chr1", &"ACGT".repeat(2000))]); // 8 kb
        let qc = LoadQc::default();
        let recs = load_fasta(f.path(), qc).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].record.id, "chr1");
        assert_eq!(recs[0].record.length, 8000);
        assert!(matches!(recs[0].status, LoadStatus::Ok));
    }

    #[test]
    fn flags_short_record() {
        let f = fasta(&[("tiny", "ACGTACGTAC")]);
        let qc = LoadQc::default();
        let recs = load_fasta(f.path(), qc).unwrap();
        assert!(matches!(recs[0].status, LoadStatus::TooShort { .. }));
    }

    #[test]
    fn flags_too_many_ns() {
        let f = fasta(&[("ny", &"N".repeat(6000))]);
        let qc = LoadQc::default();
        let recs = load_fasta(f.path(), qc).unwrap();
        assert!(matches!(recs[0].status, LoadStatus::TooManyNs { .. }));
    }

    #[test]
    fn id_stripped_to_first_whitespace_token() {
        let f = fasta(&[("chr1 extra description", &"ACGT".repeat(2000))]);
        let recs = load_fasta(f.path(), LoadQc::default()).unwrap();
        assert_eq!(recs[0].record.id, "chr1");
    }

    #[test]
    fn normalizes_lowercase_and_n() {
        let f = fasta(&[("mix", &"acgtnACGTN".repeat(800))]);
        let recs = load_fasta(f.path(), LoadQc::default()).unwrap();
        // Half of the bases are Ns → fraction 0.2, right at the limit.
        let n_frac = recs[0].record.n_fraction();
        assert!((n_frac - 0.2).abs() < 1e-9);
    }

    #[test]
    fn handles_multiple_records() {
        let f = fasta(&[
            ("a", &"ACGT".repeat(2000)),
            ("b", &"GGCC".repeat(2000)),
            ("c", "TINY"),
        ]);
        let recs = load_fasta(f.path(), LoadQc::default()).unwrap();
        assert_eq!(recs.len(), 3);
        assert!(matches!(recs[0].status, LoadStatus::Ok));
        assert!(matches!(recs[1].status, LoadStatus::Ok));
        assert!(matches!(recs[2].status, LoadStatus::TooShort { .. }));
    }
}

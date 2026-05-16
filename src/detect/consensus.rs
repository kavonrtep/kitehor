//! Column-vote consensus monomer and HOR-unit
//! (`detect_impl_plan.md §6.12`).
//!
//! At the chosen width, take the plurality base per column (A/C/G/T,
//! breaking ties alphabetically; N is excluded from the vote and used
//! as the fallback when a column has no A/C/G/T). The HOR-unit
//! consensus is built **from the HOR-unit-width wrap directly** —
//! NOT by repeating the base-width consensus (A9 in the plan).

use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Compute the consensus byte string for a sequence wrapped at
/// `width`. Drops trailing partial rows. Returns `None` if fewer
/// than 2 rows fit at this width.
pub fn consensus(seq: &[u8], width: usize) -> Option<Vec<u8>> {
    if width == 0 || seq.len() < 2 * width {
        return None;
    }
    let n_rows = seq.len() / width;
    let mut out = Vec::with_capacity(width);
    for c in 0..width {
        let mut counts = [0usize; 4]; // A, C, G, T
        for r in 0..n_rows {
            match seq[r * width + c] {
                b'A' => counts[0] += 1,
                b'C' => counts[1] += 1,
                b'G' => counts[2] += 1,
                b'T' => counts[3] += 1,
                _ => {}
            }
        }
        // argmax with alphabetic tie-breaking (A < C < G < T).
        let (max_i, _) = counts
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(&a.0)))
            .unwrap();
        let base = match max_i {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        // If all four are zero (all-N column), output N.
        if counts.iter().sum::<usize>() == 0 {
            out.push(b'N');
        } else {
            out.push(base);
        }
    }
    Some(out)
}

/// One per-array consensus block — used for batched writes.
#[derive(Debug, Clone)]
pub struct ConsensusRecord {
    pub array_id: String,
    pub monomer: Vec<u8>,
    pub hor_unit: Option<Vec<u8>>,
    pub hor_k: Option<usize>,
}

/// Write the consensus FASTA. Writes one or two records per array
/// (`.monomer` always, `.hor_unit` when present).
pub fn write_fasta(out_prefix: &Path, records: &[ConsensusRecord]) -> Result<()> {
    let path = consensus_path(out_prefix);
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {:?}", path))?;
    for r in records {
        writeln!(f, ">{}.monomer  length={}", r.array_id, r.monomer.len())?;
        write_wrapped(&mut f, &r.monomer)?;
        if let Some(unit) = &r.hor_unit {
            let k = r.hor_k.unwrap_or(0);
            writeln!(
                f,
                ">{}.hor_unit  length={}  k={}",
                r.array_id,
                unit.len(),
                k
            )?;
            write_wrapped(&mut f, unit)?;
        }
    }
    Ok(())
}

pub fn consensus_path(prefix: &Path) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(".consensus.fa");
    PathBuf::from(s)
}

fn write_wrapped(f: &mut std::fs::File, seq: &[u8]) -> std::io::Result<()> {
    for chunk in seq.chunks(80) {
        f.write_all(chunk)?;
        f.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consensus_returns_none_for_one_row() {
        // 5 rows × width 100 would need at least 200; 50 < 200 → None.
        let seq = vec![b'A'; 50];
        assert!(consensus(&seq, 100).is_none());
    }

    #[test]
    fn consensus_picks_plurality_base() {
        // 10 rows × width 4: column 0 is mostly A (7 A, 3 C), column 1
        // is all G, column 2 is all T, column 3 is N+C tie → C wins.
        let mut seq = Vec::with_capacity(40);
        for r in 0..10 {
            seq.push(if r < 7 { b'A' } else { b'C' });
            seq.push(b'G');
            seq.push(b'T');
            seq.push(if r % 2 == 0 { b'N' } else { b'C' });
        }
        let cs = consensus(&seq, 4).unwrap();
        assert_eq!(cs, b"AGTC".to_vec());
    }

    #[test]
    fn consensus_emits_n_for_all_n_column() {
        let seq: Vec<u8> = (0..40).map(|_| b'N').collect();
        let cs = consensus(&seq, 4).unwrap();
        assert_eq!(cs, b"NNNN".to_vec());
    }

    #[test]
    fn consensus_ties_break_alphabetically() {
        // 10 rows × width 1, alternating A and C → tie. A wins.
        let seq: Vec<u8> = (0..10)
            .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
            .collect();
        let cs = consensus(&seq, 1).unwrap();
        assert_eq!(cs, vec![b'A']);
    }

    #[test]
    fn write_fasta_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        let rec = ConsensusRecord {
            array_id: "arr1".into(),
            monomer: b"ACGTACGTAC".to_vec(),
            hor_unit: Some(b"ACGTACGTACTGCATGCATGC".to_vec()),
            hor_k: Some(4),
        };
        write_fasta(&prefix, std::slice::from_ref(&rec)).unwrap();
        let text = std::fs::read_to_string(consensus_path(&prefix)).unwrap();
        assert!(text.contains(">arr1.monomer  length=10"));
        assert!(text.contains("ACGTACGTAC"));
        assert!(text.contains(">arr1.hor_unit  length=21  k=4"));
        assert!(text.contains("ACGTACGTACTGCATGCATGC"));
    }

    #[test]
    fn write_fasta_multi_record() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        let recs = vec![
            ConsensusRecord {
                array_id: "a".into(),
                monomer: b"AAAA".to_vec(),
                hor_unit: None,
                hor_k: None,
            },
            ConsensusRecord {
                array_id: "b".into(),
                monomer: b"CCCC".to_vec(),
                hor_unit: Some(b"CCCCGGGG".to_vec()),
                hor_k: Some(2),
            },
        ];
        write_fasta(&prefix, &recs).unwrap();
        let s = std::fs::read_to_string(consensus_path(&prefix)).unwrap();
        assert!(s.contains(">a.monomer"));
        assert!(s.contains(">b.monomer"));
        assert!(s.contains(">b.hor_unit"));
    }
}

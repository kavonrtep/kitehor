//! FASTA-like periodogram bundle writer.
//!
//! Mirrors the data shape of TideCluster's in-memory `profile_list`
//! (`tarean/kite.R`): per input record, the dense neighbour-distance
//! histogram `H[d]` for `d ∈ [1, N]` (raw integer counts) and the
//! composition-matched random background envelope `bg[d]`.
//!
//! TideCluster never persists this data as text — only as `peaks_list.RDS`.
//! kitehor invents a single bundle file with FASTA-style records so the
//! output is easy to load from any language and easy to concatenate.
//!
//! Format (one file per run, passed via `--periodogram <path>`):
//!
//! ```text
//! >case_id|H length=<N> kmer=<K>
//! <H[1]> <H[2]> ... <H[N]>
//! >case_id|bg length=<N> kmer=<K>
//! <bg[1]> <bg[2]> ... <bg[N]>
//! ```
//!
//! - Two records per input sequence (`|H` raw counts, `|bg` smoothed
//!   background).
//! - Header fields after the record id are space-separated `key=value`
//!   pairs. A parser that splits on whitespace gets `id|channel` first
//!   then any number of metadata tokens.
//! - Values are whitespace-separated on a single line per record. `H[d]`
//!   formatted as integer (`{:.0}`; the upstream histogram is integral),
//!   `bg[d]` formatted with 6 fractional digits.
//! - The vector covers `d = 1..=N` where `N` is the record's
//!   `length_bp`. `H[0]` is unused upstream and not emitted.
//! - Records whose array failed kite analysis (empty `profile`) are
//!   skipped silently.

use anyhow::{Context, Result};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::kite::{KiteConfig, KiteResult};

/// Write a periodogram bundle covering every `KiteResult` in `results`.
///
/// `kite_cfg` is read only to record `kmer=<K>` in each header.
pub fn write_periodogram_bundle(
    path: &Path,
    results: &[KiteResult],
    kite_cfg: &KiteConfig,
) -> Result<usize> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating parent of {:?}", path))?;
        }
    }
    let f = std::fs::File::create(path).with_context(|| format!("creating {:?}", path))?;
    let mut w = BufWriter::new(f);

    let mut n_written = 0usize;
    for r in results {
        // Skip records with empty profile (kite refused to analyze, e.g.
        // L < k + 2). Emitting an empty record would force consumers to
        // special-case it; skipping keeps the format dense.
        if r.profile.is_empty() || r.background.is_empty() {
            continue;
        }
        let n = r.length_bp;
        write_channel(&mut w, &r.array_id, "H", n, kite_cfg.k, &r.profile, false)?;
        write_channel(
            &mut w,
            &r.array_id,
            "bg",
            n,
            kite_cfg.k,
            &r.background,
            true,
        )?;
        n_written += 1;
    }
    w.flush().with_context(|| format!("flushing {:?}", path))?;
    Ok(n_written)
}

fn write_channel<W: Write>(
    w: &mut W,
    array_id: &str,
    channel: &str,
    n: usize,
    k: usize,
    values: &[f64],
    fractional: bool,
) -> Result<()> {
    writeln!(w, ">{}|{} length={} kmer={}", array_id, channel, n, k)?;
    // Emit d=1..=n. Sources are length L+1 with index 0 unused; clamp
    // defensively in case a future change shortens them.
    let mut first = true;
    for d in 1..=n {
        if !first {
            w.write_all(b" ")?;
        }
        first = false;
        let v = values.get(d).copied().unwrap_or(0.0);
        if fractional {
            write!(w, "{:.6}", v)?;
        } else {
            // H is the result of integer counts — round to nearest and
            // emit without trailing zeros.
            write!(w, "{}", v.round() as i64)?;
        }
    }
    w.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kite::{KitePeak, KiteResult};

    fn synth_result(id: &str, n: usize) -> KiteResult {
        // d=0..n; we only emit d=1..=n.
        let mut profile = vec![0.0f64; n + 1];
        let mut background = vec![0.0f64; n + 1];
        for d in 1..=n {
            profile[d] = (d as f64).rem_euclid(7.0); // deterministic small ints
            background[d] = (d as f64) * 0.125; // exact in f64
        }
        KiteResult {
            array_id: id.into(),
            length_bp: n,
            peaks: Vec::<KitePeak>::new(),
            profile,
            background,
        }
    }

    #[test]
    fn round_trip_two_records() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path();
        let results = vec![synth_result("rec_a", 12), synth_result("rec_b", 5)];
        let cfg = KiteConfig::default();
        let n = write_periodogram_bundle(path, &results, &cfg).unwrap();
        assert_eq!(n, 2);

        let body = std::fs::read_to_string(path).unwrap();
        let mut lines = body.lines();
        // rec_a|H
        let h_a = lines.next().unwrap();
        assert_eq!(h_a, ">rec_a|H length=12 kmer=6");
        let h_a_vals: Vec<i64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse::<i64>().unwrap())
            .collect();
        assert_eq!(h_a_vals.len(), 12);
        for (d, v) in h_a_vals.iter().enumerate() {
            let expected = ((d + 1) as f64).rem_euclid(7.0).round() as i64;
            assert_eq!(*v, expected, "H mismatch at d={}", d + 1);
        }
        // rec_a|bg
        let bg_a = lines.next().unwrap();
        assert_eq!(bg_a, ">rec_a|bg length=12 kmer=6");
        let bg_a_vals: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        assert_eq!(bg_a_vals.len(), 12);
        for (d, v) in bg_a_vals.iter().enumerate() {
            let expected = ((d + 1) as f64) * 0.125;
            assert!((v - expected).abs() < 1e-9);
        }
        // rec_b|H header next.
        let h_b = lines.next().unwrap();
        assert_eq!(h_b, ">rec_b|H length=5 kmer=6");
    }

    #[test]
    fn empty_profile_record_skipped() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let kr = KiteResult {
            array_id: "too_short".into(),
            length_bp: 0,
            peaks: Vec::new(),
            profile: Vec::new(),
            background: Vec::new(),
        };
        let cfg = KiteConfig::default();
        let n = write_periodogram_bundle(tmp.path(), &[kr], &cfg).unwrap();
        assert_eq!(n, 0);
        let body = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(body.is_empty(), "expected empty file, got {:?}", body);
    }

    #[test]
    fn header_records_kmer_size() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let kr = synth_result("x", 3);
        let cfg = KiteConfig {
            k: 7,
            ..KiteConfig::default()
        };
        write_periodogram_bundle(tmp.path(), &[kr], &cfg).unwrap();
        let body = std::fs::read_to_string(tmp.path()).unwrap();
        assert!(body.starts_with(">x|H length=3 kmer=7"));
    }
}

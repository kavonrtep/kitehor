//! Irregularity scan — Rust port of
//! `tools/rule_proto/irregularity_v2.py` (distance-residual approach
//! with phase-bin k-mer clustering and explicit dropout/indel split).
//!
//! Quantifies indel-like phase-shift events in a tandem repeat array.
//! For each register-locked k-mer, consecutive distances are compared
//! to the nearest multiple of the k-mer's modal spacing; non-zero
//! residuals signal coordinated displacement (indel events). The
//! algorithm runs at the structural period reported by kite.
//!
//! Caveat: only non-P-multiple length disruptions are detected.
//! Whole-monomer copy gain/loss is invisible.
//!
//! See `docs/irregularity_and_subrepeat_v0_12.md` §3 for the full
//! algorithm description and `ground_truth3` calibration.

pub mod io;
pub mod scan;

use anyhow::Result;
use std::path::Path;

pub use scan::{analyse_record, Config, RecordResult};

/// Standalone subcommand entry point —
/// `kitehor irregularity --fasta <fa> --kite <kite.tsv> -o <prefix>`.
///
/// Reads the FASTA + kite-summary TSV (per-record top period), runs
/// the scan in parallel over records, and emits
/// `<prefix>.irregularity.tsv` (14 columns matching the Python
/// prototype).
pub fn run_subcommand(
    fasta: &Path,
    kite_summary: &Path,
    out_prefix: &Path,
    cfg: &Config,
) -> Result<usize> {
    use rayon::prelude::*;

    let records = crate::ssr::io::read_fasta_ordered(fasta)?;
    let periods = io::read_kite_top_periods(kite_summary)?;

    let results: Vec<RecordResult> = records
        .par_iter()
        .map(|(rid, seq)| {
            let p = periods.get(rid).copied();
            scan::analyse_record(rid, seq, p, cfg)
        })
        .collect();

    let path = io::irregularity_path(out_prefix);
    io::write_irregularity(&path, &results)?;
    Ok(results.len())
}

//! Rule-based HOR/simple_tr/unresolved classifier — Rust port of
//! `tools/rule_proto/rule_proto.py`.
//!
//! Operates on the long-format kite peaks TSV (one row per kept peak per
//! record). For each record we cluster peaks by relative period gap,
//! score-filter the clusters, and apply a first-match-wins decision tree
//! (Case A: top is `k × shorter cluster`; Case B: top is the founder
//! with a non-monotonic bump at `k × top`; fallbacks: monotonic-multiples
//! simple_tr, lone-significant-cluster simple_tr, else unresolved).
//!
//! Two integration paths:
//!
//! 1. Standalone `kitehor rule-classify <peaks.tsv> -o <prefix>` →
//!    emits `<prefix>.verdicts.tsv` (10 columns, `%.6g` float format —
//!    byte-equivalent with the Python prototype).
//! 2. Library call from `kite-periodicity --classify` and from
//!    `detect::auto_periods`: `classify(&KiteResult, &Config) -> Verdict`
//!    + the `LegacyVerdict` adapter consumed by `emit_periods::build_rows`.

pub mod cluster;
pub mod decide;
pub mod io;

pub use cluster::{cluster_peaks, Cluster};
pub use decide::{decide, harmonic_confirms_hor, Config, Verdict, VerdictKind};

use crate::kite::KiteResult;

/// Apply the rule classifier to a single record's kite peaks.
pub fn classify(kite: &KiteResult, cfg: &Config) -> Verdict {
    let rows: Vec<decide::PeakRow> = kite
        .peaks
        .iter()
        .enumerate()
        .map(|(i, p)| decide::PeakRow {
            rank: (i + 1) as u32,
            period: p.period,
            score2_norm: p.score2_norm,
        })
        .collect();
    decide::decide(&kite.array_id, &rows, cfg)
}

/// Backward-compatible enum consumed by `emit_periods::build_rows` and
/// the `kite-periodicity --classify` per-row output formatter. Adapts
/// the richer [`Verdict`] to the small shape the bridge module needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LegacyVerdict {
    /// No kite peaks survived clustering (`no_clusters`).
    NoSignal,
    /// Classifier could not decide between simple_tr and HOR
    /// (`no_multiples` with no lone significant cluster).
    Unresolved,
    /// `simple_tr` — either `monotonic_multiples` or
    /// `lone_significant_cluster`. The monomer is the founder period.
    Tandem { monomer_bp: usize },
    /// `hor` — founder × k = tile.
    Hor {
        founder: usize,
        tile: usize,
        k: usize,
        /// Algorithmic confidence (0..=1). Sourced from
        /// [`Verdict::confidence`].
        share: f64,
    },
}

impl LegacyVerdict {
    pub fn from_verdict(v: &Verdict) -> Self {
        match v.kind {
            VerdictKind::Hor => LegacyVerdict::Hor {
                founder: v.founder.map(|f| f.round() as usize).unwrap_or(0),
                tile: v.tile.map(|t| t.round() as usize).unwrap_or(0),
                k: v.multiplicity.unwrap_or(0) as usize,
                share: v.confidence.unwrap_or(0.0),
            },
            VerdictKind::SimpleTr => LegacyVerdict::Tandem {
                monomer_bp: v.founder.map(|f| f.round() as usize).unwrap_or(0),
            },
            VerdictKind::Unresolved => {
                if v.founder.is_none() {
                    LegacyVerdict::NoSignal
                } else {
                    LegacyVerdict::Unresolved
                }
            }
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            LegacyVerdict::Hor { .. } => "hor",
            LegacyVerdict::Tandem { .. } => "simple_tr",
            LegacyVerdict::Unresolved => "unresolved",
            LegacyVerdict::NoSignal => "no_signal",
        }
    }

    pub fn founder(&self) -> Option<usize> {
        match self {
            LegacyVerdict::Hor { founder, .. } => Some(*founder),
            _ => None,
        }
    }

    pub fn tile(&self) -> Option<usize> {
        match self {
            LegacyVerdict::Hor { tile, .. } => Some(*tile),
            LegacyVerdict::Tandem { monomer_bp } => Some(*monomer_bp),
            _ => None,
        }
    }

    pub fn multiplicity(&self) -> Option<usize> {
        match self {
            LegacyVerdict::Hor { k, .. } => Some(*k),
            LegacyVerdict::Tandem { .. } => Some(1),
            _ => None,
        }
    }

    pub fn share(&self) -> Option<f64> {
        match self {
            LegacyVerdict::Hor { share, .. } => Some(*share),
            _ => None,
        }
    }
}

/// Entry point for `kitehor rule-classify <peaks.tsv> -o <prefix>`.
pub fn run_subcommand(
    peaks: &std::path::Path,
    out_prefix: &std::path::Path,
    cfg: &Config,
    dump_clusters: Option<&std::path::Path>,
) -> anyhow::Result<usize> {
    let grouped = io::read_peaks_grouped(peaks)?;
    let mut verdicts: Vec<Verdict> = Vec::with_capacity(grouped.len());
    if let Some(dir) = dump_clusters {
        std::fs::create_dir_all(dir)?;
    }
    for (case_id, rows) in &grouped {
        let (verdict, clusters) = decide::decide_with_clusters(case_id, rows, cfg);
        verdicts.push(verdict);
        if let Some(dir) = dump_clusters {
            if !clusters.is_empty() {
                io::write_clusters_dump(dir, case_id, &clusters)?;
            }
        }
    }
    let verdicts_path = io::verdicts_path(out_prefix);
    io::write_verdicts(&verdicts_path, &verdicts)?;
    Ok(verdicts.len())
}

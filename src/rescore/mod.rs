//! kite rescore — sampled pairwise tile-identity rescoring of kite peaks.
//!
//! Reads a kite-periodicity peaks TSV (long format) plus the source FASTA(s),
//! and for each (record, candidate period) computes the median pairwise
//! identity between adjacent tiles. Appends `identity_med`, `identity_iqr`,
//! `identity_p25`, `identity_n` columns to the peaks TSV.
//!
//! Wired as the `kitehor rescore` subcommand. The metric is additive only —
//! downstream stages (rule-classify, analyze) see the new columns but
//! continue to drive decisions from kite's `score2_norm`.
//!
//! Sampling is adjacent-tile only (`d=1`). Cross-distance probing
//! (`d=2,3,…`) is reserved for a future flag; the current goal is HOR vs.
//! monomer rescoring where adjacent comparison is sufficient.

pub mod aligner;
pub mod io;
pub mod sample;

use crate::io::{load_fasta, LoadQc, LoadStatus};
use crate::sequence::ArrayRecord;
use ahash::AHashMap;
use anyhow::{anyhow, Context, Result};
use log::info;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

pub use sample::SampleConfig;

/// Rescore stage configuration.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub samples: usize,
    pub slop: usize,
    /// Indel-deviation tolerance for the banded kernel. `0` = auto, which
    /// resolves to `max(20, 2 · slop)` at run time.
    pub band: usize,
    pub max_n_frac: f64,
    pub max_retries: usize,
    pub min_period: usize,
    /// Skip candidates with `period > max_period`. `0` = unlimited.
    pub max_period: usize,
    pub seed: u64,
    /// 0 = all peaks; otherwise only rows with `rank ≤ top_n` are rescored.
    pub top_n: usize,
    /// Per-cell costs used by the alignment kernel. Defaults to unit
    /// edit distance (mismatch=1, gap=1). Match cost is always 0.
    pub scoring: aligner::ScoringConfig,
    pub load_qc: LoadQc,
    pub force: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            samples: 200,
            slop: 10,
            band: 0,
            max_n_frac: 0.05,
            max_retries: 3,
            min_period: 20,
            max_period: 5000,
            seed: 42,
            top_n: 10,
            scoring: aligner::ScoringConfig::default(),
            load_qc: LoadQc::default(),
            force: false,
        }
    }
}

impl Config {
    /// Resolve the effective band given the `band == 0` auto convention.
    pub fn resolved_band(&self) -> usize {
        if self.band == 0 {
            (2 * self.slop).max(20)
        } else {
            self.band
        }
    }
}

impl Config {
    fn sample_cfg(&self) -> SampleConfig {
        SampleConfig {
            k: self.samples,
            slop: self.slop,
            max_n_frac: self.max_n_frac,
            max_retries: self.max_retries,
            seed: self.seed,
        }
    }
}

/// Per-row identity statistics. `None` fields render as `NA`.
#[derive(Debug, Clone, Copy)]
pub struct RowStats {
    pub identity_med: Option<f64>,
    pub identity_iqr: Option<f64>,
    pub identity_p25: Option<f64>,
    pub identity_n: usize,
}

impl RowStats {
    pub fn na() -> Self {
        Self {
            identity_med: None,
            identity_iqr: None,
            identity_p25: None,
            identity_n: 0,
        }
    }

    fn format_row(&self) -> String {
        let f = |o: Option<f64>| {
            o.map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".into())
        };
        format!(
            "{}\t{}\t{}\t{}",
            f(self.identity_med),
            f(self.identity_iqr),
            f(self.identity_p25),
            self.identity_n
        )
    }
}

/// Stage entry point used by the CLI.
///
/// Reads `<peaks_in>`, loads sequences from `fastas`, computes the four
/// identity columns per row, and writes the augmented file to
/// `<peaks_out>`. Returns the number of rows processed.
pub fn run_subcommand(
    fastas: &[PathBuf],
    peaks_in: &Path,
    peaks_out: &Path,
    cfg: &Config,
) -> Result<usize> {
    if peaks_in == peaks_out && !cfg.force {
        return Err(anyhow!(
            "refusing to overwrite input peaks file {:?}; pass --force or pick a different -o prefix",
            peaks_in
        ));
    }
    if peaks_out.exists() && !cfg.force {
        return Err(anyhow!(
            "output {:?} already exists; pass --force to overwrite",
            peaks_out
        ));
    }

    let mut records: AHashMap<String, ArrayRecord> = AHashMap::new();
    for path in fastas {
        let loaded = load_fasta(path, cfg.load_qc)?;
        for lr in loaded {
            // Only Ok records get an entry; non-Ok records (TooShort,
            // TooManyNs) yield NA rows downstream via the missing-lookup
            // path.
            if matches!(lr.status, LoadStatus::Ok) {
                records.insert(lr.record.id.clone(), lr.record);
            }
        }
    }

    let loaded = io::load_peaks(peaks_in)?;
    let sample_cfg = cfg.sample_cfg();
    let band = cfg.resolved_band();
    let top_n = if cfg.top_n == 0 {
        usize::MAX
    } else {
        cfg.top_n
    };

    let max_period_eff = if cfg.max_period == 0 {
        usize::MAX
    } else {
        cfg.max_period
    };

    let n_records = records.len();
    let n_rows = loaded.rows.len();
    let n_to_rescore = loaded
        .rows
        .iter()
        .filter(|r| {
            r.rank <= top_n
                && r.period >= cfg.min_period
                && r.period <= max_period_eff
                && records.contains_key(&r.case_id)
        })
        .count();
    let threads = rayon::current_num_threads();

    info!(
        "rescore: loaded {} record(s), {} peak row(s); {} to rescore (filters: min_period={}, max_period={}, top_n={})",
        n_records,
        n_rows,
        n_to_rescore,
        cfg.min_period,
        if cfg.max_period == 0 { "all".to_string() } else { cfg.max_period.to_string() },
        if cfg.top_n == 0 { "all".to_string() } else { cfg.top_n.to_string() },
    );
    info!(
        "rescore: K={} slop={} band={} mismatch_cost={} gap_cost={} max_retries={} seed={} threads={}",
        cfg.samples,
        cfg.slop,
        band,
        cfg.scoring.mismatch_cost,
        cfg.scoring.gap_cost,
        cfg.max_retries,
        cfg.seed,
        threads,
    );

    let start = Instant::now();
    let processed = AtomicUsize::new(0);
    let last_log_sec = AtomicU64::new(0);
    const LOG_INTERVAL_SEC: u64 = 10;

    let stats: Vec<RowStats> = loaded
        .rows
        .par_iter()
        .map_init(aligner::Scratch::new, |scratch, row| {
            let (result, did_rescore) =
                if row.rank > top_n || row.period < cfg.min_period || row.period > max_period_eff {
                    (RowStats::na(), false)
                } else if let Some(record) = records.get(&row.case_id) {
                    let r = rescore_one(
                        &record.seq,
                        row.period,
                        &row.case_id,
                        &sample_cfg,
                        band,
                        &cfg.scoring,
                        scratch,
                    );
                    (r, true)
                } else {
                    (RowStats::na(), false)
                };

            if did_rescore {
                let done = processed.fetch_add(1, Ordering::Relaxed) + 1;
                let elapsed = start.elapsed();
                let elapsed_sec = elapsed.as_secs();
                let prev = last_log_sec.load(Ordering::Relaxed);
                if elapsed_sec >= prev + LOG_INTERVAL_SEC
                    && last_log_sec
                        .compare_exchange(prev, elapsed_sec, Ordering::Relaxed, Ordering::Relaxed)
                        .is_ok()
                {
                    let secs = elapsed.as_secs_f64();
                    let pct = if n_to_rescore == 0 {
                        100.0
                    } else {
                        100.0 * done as f64 / n_to_rescore as f64
                    };
                    let rate = if secs > 0.0 { done as f64 / secs } else { 0.0 };
                    let remaining = n_to_rescore.saturating_sub(done);
                    let eta = if rate > 0.0 {
                        remaining as f64 / rate
                    } else {
                        0.0
                    };
                    info!(
                        "rescore: {}/{} ({:.1}%) elapsed={:.0}s rate={:.0}/s eta={:.0}s",
                        done, n_to_rescore, pct, secs, rate, eta,
                    );
                }
            }

            result
        })
        .collect();

    let total_elapsed = start.elapsed();
    let total_na = stats.iter().filter(|s| s.identity_med.is_none()).count();
    let filtered = n_rows.saturating_sub(n_to_rescore);
    let kernel_na = total_na.saturating_sub(filtered);
    let mut ns: Vec<usize> = stats
        .iter()
        .filter_map(|s| (s.identity_n > 0).then_some(s.identity_n))
        .collect();
    ns.sort_unstable();
    let med_n = if ns.is_empty() { 0 } else { ns[ns.len() / 2] };
    info!(
        "rescore: done in {:.1}s — rescored {}, filtered {}, kernel-NA {}, identity_n median={}",
        total_elapsed.as_secs_f64(),
        n_to_rescore.saturating_sub(kernel_na),
        filtered,
        kernel_na,
        med_n,
    );

    let file =
        std::fs::File::create(peaks_out).with_context(|| format!("creating {:?}", peaks_out))?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "{}\tidentity_med\tidentity_iqr\tidentity_p25\tidentity_n",
        loaded.header
    )?;
    for (row, s) in loaded.rows.iter().zip(stats.iter()) {
        writeln!(w, "{}\t{}", row.line, s.format_row())?;
    }
    w.flush()?;

    Ok(loaded.rows.len())
}

/// Compute median/IQR/p25/n for one (record, period) over `cfg.k` pairs.
pub fn rescore_one(
    seq: &[u8],
    period: usize,
    case_id: &str,
    cfg: &SampleConfig,
    band: usize,
    scoring: &aligner::ScoringConfig,
    scratch: &mut aligner::Scratch,
) -> RowStats {
    let pairs = sample::sample_pairs(seq, period, case_id, cfg);
    if pairs.is_empty() {
        return RowStats::na();
    }
    let mut ids: Vec<f64> = pairs
        .iter()
        .map(|p| {
            let a = &seq[p.a_start..p.a_end];
            let b = &seq[p.b_start..p.b_end];
            let d = aligner::semiglobal_edit_distance_banded(a, b, band, scoring, scratch);
            aligner::identity_from_distance(d, period)
        })
        .collect();
    ids.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = quantile_sorted(&ids, 0.5);
    let p25 = quantile_sorted(&ids, 0.25);
    let p75 = quantile_sorted(&ids, 0.75);
    RowStats {
        identity_med: Some(med),
        identity_iqr: Some((p75 - p25).max(0.0)),
        identity_p25: Some(p25),
        identity_n: ids.len(),
    }
}

/// Linear-interpolated quantile of a pre-sorted slice. `q ∈ [0, 1]`.
fn quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted[0];
    }
    let pos = q * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_basic() {
        let s = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        assert!((quantile_sorted(&s, 0.5) - 0.5).abs() < 1e-12);
        assert!((quantile_sorted(&s, 0.25) - 0.25).abs() < 1e-12);
        assert!((quantile_sorted(&s, 0.75) - 0.75).abs() < 1e-12);
        assert!((quantile_sorted(&s, 0.0) - 0.0).abs() < 1e-12);
        assert!((quantile_sorted(&s, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn na_row_formats_as_na() {
        let s = RowStats::na();
        assert_eq!(s.format_row(), "NA\tNA\tNA\t0");
    }

    #[test]
    fn rescore_one_short_array_is_na() {
        let seq = vec![b'A'; 100];
        let mut scratch = aligner::Scratch::new();
        let s = rescore_one(
            &seq,
            200,
            "x",
            &SampleConfig::default(),
            20,
            &aligner::ScoringConfig::default(),
            &mut scratch,
        );
        assert!(s.identity_med.is_none());
        assert_eq!(s.identity_n, 0);
    }

    #[test]
    fn rescore_one_perfect_repeat_is_high_identity() {
        // 100 copies of a 100 bp monomer = 10kb of perfect tandem repeat.
        let monomer: Vec<u8> = b"ACGT".iter().cycle().take(100).copied().collect();
        let mut seq = Vec::new();
        for _ in 0..100 {
            seq.extend_from_slice(&monomer);
        }
        let mut scratch = aligner::Scratch::new();
        let s = rescore_one(
            &seq,
            100,
            "x",
            &SampleConfig {
                k: 50,
                slop: 10,
                max_n_frac: 0.05,
                max_retries: 3,
                seed: 1,
            },
            20,
            &aligner::ScoringConfig::default(),
            &mut scratch,
        );
        // Perfect tandem ⇒ identity = 1.0 across all pairs.
        assert!(
            s.identity_med.unwrap() > 0.99,
            "expected ~1.0, got {:?}",
            s.identity_med
        );
        assert_eq!(s.identity_n, 50);
    }

    #[test]
    fn resolved_band_auto_uses_max_20_2slop() {
        let cfg = Config {
            slop: 5,
            band: 0,
            ..Config::default()
        };
        // 2*5 = 10, max(20, 10) = 20
        assert_eq!(cfg.resolved_band(), 20);
        let cfg = Config {
            slop: 50,
            band: 0,
            ..Config::default()
        };
        // 2*50 = 100, max(20, 100) = 100
        assert_eq!(cfg.resolved_band(), 100);
        let cfg = Config {
            slop: 5,
            band: 7,
            ..Config::default()
        };
        // Explicit override
        assert_eq!(cfg.resolved_band(), 7);
    }
}

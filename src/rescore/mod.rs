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

/// Sub-period phantom-flag thresholds. The flag fires when the per-pair
/// alignment shift relative to the natural mapping is systematically
/// non-zero — evidence that the claimed period is a sub-tile of the true
/// repeat. See `docs/rescore.md` for the worked example.
#[derive(Debug, Clone, Copy)]
pub struct PhantomConfig {
    /// Only include pairs with identity ≥ this threshold in the shift
    /// aggregate. Below this, the per-pair "best offset" is noise.
    pub identity_min: f64,
    /// Minimum number of high-identity pairs required for shift_med to
    /// be non-NA.
    pub min_pairs: usize,
    /// Phantom fires when `|shift_med| / period > tol_frac`.
    pub tol_frac: f64,
    /// Phantom fires only when the fraction of pairs whose shift is
    /// within ±1 bp of `shift_med` exceeds this.
    pub consistency_min: f64,
}

impl Default for PhantomConfig {
    fn default() -> Self {
        Self {
            identity_min: 0.5,
            min_pairs: 5,
            tol_frac: 0.05,
            consistency_min: 0.5,
        }
    }
}

/// Subrepeat heuristic thresholds. Flags candidate periods whose identity
/// distribution is bimodal — high in part of the array (where the short
/// motif tandemly repeats inside the founder monomer) and near-random
/// elsewhere — without being phantoms or real periods.
#[derive(Debug, Clone, Copy)]
pub struct SubrepeatConfig {
    /// Subrepeat fires only when `identity_p75 ≥ p75_min`. Captures
    /// "some pairs hit hard".
    pub p75_min: f64,
    /// Subrepeat fires only when `identity_iqr ≥ iqr_min`. Captures
    /// "bimodal spread".
    pub iqr_min: f64,
    /// Subrepeat fires only when `identity_med < med_max`. Distinguishes
    /// the bimodal case from a real period that happens to have a wide
    /// IQR.
    pub med_max: f64,
    /// Per-pair identity threshold for the `coverage_frac` column —
    /// pairs at or above this count as "hits". Default 0.70.
    pub coverage_threshold: f64,
    /// Subrepeat fires only when `coverage_frac ≥ cov_min`. Excludes
    /// noise periods where no pairs actually hit. Default 0.10.
    pub cov_min: f64,
    /// Subrepeat fires only when `coverage_frac ≤ cov_max`. Excludes
    /// real periods where most pairs hit. Default 0.50.
    pub cov_max: f64,
    /// Minimum `identity_med` for a row to qualify as the per-record
    /// founder when applying the founder gate (subrepeat must have
    /// period < founder period). Default 0.70.
    pub founder_id_min: f64,
}

impl Default for SubrepeatConfig {
    fn default() -> Self {
        Self {
            p75_min: 0.70,
            iqr_min: 0.15,
            med_max: 0.70,
            coverage_threshold: 0.70,
            cov_min: 0.10,
            cov_max: 0.50,
            founder_id_min: 0.70,
        }
    }
}

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
    pub phantom: PhantomConfig,
    pub subrepeat: SubrepeatConfig,
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
            phantom: PhantomConfig::default(),
            subrepeat: SubrepeatConfig::default(),
            load_qc: LoadQc::default(),
            force: false,
        }
    }
}

impl Config {
    /// Resolve the effective band for a candidate period given the
    /// `band == 0` auto convention. The auto formula is
    /// `max(20, 2·slop, ⌈0.02·period⌉)`: short periods keep the old
    /// `max(20, 2·slop)` default unchanged, while long monomers
    /// (P > ~1000) get a band scaled to ~2 % of the period — enough to
    /// cover realistic internal indel drift in centromeric satellites
    /// without the band-saturation artifact that otherwise drops
    /// identity_med to ~0.5 on long-monomer arrays.
    ///
    /// A user-set `--band N` (non-zero) bypasses the formula.
    pub fn resolved_band(&self, period: usize) -> usize {
        if self.band == 0 {
            let slop_floor = (2 * self.slop).max(20);
            let period_relative = ((period as f64) * 0.02).ceil() as usize;
            slop_floor.max(period_relative)
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
    /// Median shift (in bp) between the optimal alignment offset and the
    /// natural mapping, aggregated over pairs with identity above
    /// `PhantomConfig::identity_min`. `None` when fewer than
    /// `min_pairs` high-identity pairs were available.
    pub shift_med: Option<i32>,
    /// Fraction of high-identity pairs whose shift is within ±1 bp of
    /// `shift_med`. `None` when `shift_med` is `None`.
    pub shift_consistency: Option<f64>,
    /// Derived phantom flag — `true` when the candidate period is likely
    /// a sub-tile of a longer real period (see `PhantomConfig`).
    pub phantom: Option<bool>,
    /// Derived subrepeat flag — `true` when the candidate period
    /// appears to be a short tandem motif localized in part of the
    /// array (typically within the founder monomer). Always `false`
    /// when `phantom == Some(true)`.
    pub subrepeat: Option<bool>,
    /// Fraction of pairs with identity at or above
    /// `SubrepeatConfig::coverage_threshold`. Independent diagnostic of
    /// what fraction of the array the candidate period actually tiles.
    pub coverage_frac: Option<f64>,
}

impl RowStats {
    pub fn na() -> Self {
        Self {
            identity_med: None,
            identity_iqr: None,
            identity_p25: None,
            identity_n: 0,
            shift_med: None,
            shift_consistency: None,
            phantom: None,
            subrepeat: None,
            coverage_frac: None,
        }
    }

    fn format_row(&self) -> String {
        let f = |o: Option<f64>| {
            o.map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".into())
        };
        let fi = |o: Option<i32>| o.map(|v| v.to_string()).unwrap_or_else(|| "NA".into());
        let fb = |o: Option<bool>| o.map(|v| v.to_string()).unwrap_or_else(|| "NA".into());
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f(self.identity_med),
            f(self.identity_iqr),
            f(self.identity_p25),
            self.identity_n,
            fi(self.shift_med),
            f(self.shift_consistency),
            fb(self.phantom),
            fb(self.subrepeat),
            f(self.coverage_frac),
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
    let band_display = if cfg.band == 0 {
        "auto (max(20, 2·slop, ⌈0.02·P⌉))".to_string()
    } else {
        cfg.band.to_string()
    };
    info!(
        "rescore: K={} slop={} band={} mismatch_cost={} gap_cost={} max_retries={} seed={} threads={}",
        cfg.samples,
        cfg.slop,
        band_display,
        cfg.scoring.mismatch_cost,
        cfg.scoring.gap_cost,
        cfg.max_retries,
        cfg.seed,
        threads,
    );
    info!(
        "rescore: phantom_flag identity_min={} min_pairs={} tol_frac={} consistency_min={}",
        cfg.phantom.identity_min,
        cfg.phantom.min_pairs,
        cfg.phantom.tol_frac,
        cfg.phantom.consistency_min,
    );
    info!(
        "rescore: subrepeat_flag p75_min={} iqr_min={} med_max={} coverage_threshold={} cov_min={} cov_max={} founder_id_min={}",
        cfg.subrepeat.p75_min,
        cfg.subrepeat.iqr_min,
        cfg.subrepeat.med_max,
        cfg.subrepeat.coverage_threshold,
        cfg.subrepeat.cov_min,
        cfg.subrepeat.cov_max,
        cfg.subrepeat.founder_id_min,
    );

    let start = Instant::now();
    let processed = AtomicUsize::new(0);
    let last_log_sec = AtomicU64::new(0);
    const LOG_INTERVAL_SEC: u64 = 10;

    let mut stats: Vec<RowStats> = loaded
        .rows
        .par_iter()
        .map_init(aligner::Scratch::new, |scratch, row| {
            let (result, did_rescore) =
                if row.rank > top_n || row.period < cfg.min_period || row.period > max_period_eff {
                    (RowStats::na(), false)
                } else if let Some(record) = records.get(&row.case_id) {
                    let band = cfg.resolved_band(row.period);
                    let r = rescore_one(
                        &record.seq,
                        row.period,
                        &row.case_id,
                        &sample_cfg,
                        band,
                        &cfg.scoring,
                        &cfg.phantom,
                        &cfg.subrepeat,
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

    // Founder gate: subrepeat must have period < founder period within
    // the same record. Founder = lowest-rank row with identity_med ≥
    // subrepeat.founder_id_min AND phantom != true.
    let founder_overrides =
        enforce_subrepeat_founder_gate(&loaded.rows, &mut stats, cfg.subrepeat.founder_id_min);

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
        "rescore: done in {:.1}s — rescored {}, filtered {}, kernel-NA {}, identity_n median={}, founder-gated {}",
        total_elapsed.as_secs_f64(),
        n_to_rescore.saturating_sub(kernel_na),
        filtered,
        kernel_na,
        med_n,
        founder_overrides,
    );

    let file =
        std::fs::File::create(peaks_out).with_context(|| format!("creating {:?}", peaks_out))?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "{}\tidentity_med\tidentity_iqr\tidentity_p25\tidentity_n\tshift_med\tshift_consistency\tphantom\tsubrepeat\tcoverage_frac",
        loaded.header
    )?;
    for (row, s) in loaded.rows.iter().zip(stats.iter()) {
        writeln!(w, "{}\t{}", row.line, s.format_row())?;
    }
    w.flush()?;

    Ok(loaded.rows.len())
}

/// Compute median/IQR/p25/n + shift_med/consistency/phantom for one
/// (record, period) over `cfg.k` pairs.
#[allow(clippy::too_many_arguments)]
pub fn rescore_one(
    seq: &[u8],
    period: usize,
    case_id: &str,
    cfg: &SampleConfig,
    band: usize,
    scoring: &aligner::ScoringConfig,
    phantom_cfg: &PhantomConfig,
    subrepeat_cfg: &SubrepeatConfig,
    scratch: &mut aligner::Scratch,
) -> RowStats {
    let pairs = sample::sample_pairs(seq, period, case_id, cfg);
    if pairs.is_empty() {
        return RowStats::na();
    }

    // Per-pair: (identity, optional shift). Shift = j_end - period - slop;
    // the natural mapping has shift = 0.
    let per_pair: Vec<(f64, Option<i32>)> = pairs
        .iter()
        .map(|p| {
            let a = &seq[p.a_start..p.a_end];
            let b = &seq[p.b_start..p.b_end];
            let r = aligner::semiglobal_edit_distance_banded(a, b, band, scoring, scratch);
            let identity = aligner::identity_from_distance(r.distance, period);
            let shift = if r.j_end == aligner::J_END_NONE {
                None
            } else {
                Some(r.j_end as i32 - period as i32 - cfg.slop as i32)
            };
            (identity, shift)
        })
        .collect();

    // Identity aggregate over all K pairs (unchanged behaviour).
    let mut ids: Vec<f64> = per_pair.iter().map(|(id, _)| *id).collect();
    ids.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = quantile_sorted(&ids, 0.5);
    let p25 = quantile_sorted(&ids, 0.25);
    let p75 = quantile_sorted(&ids, 0.75);

    // Shift aggregate over high-identity pairs only.
    let (shift_med, shift_consistency, phantom) = aggregate_shifts(&per_pair, period, phantom_cfg);

    let iqr = (p75 - p25).max(0.0);
    let coverage_frac = if ids.is_empty() {
        0.0
    } else {
        ids.iter()
            .filter(|x| **x >= subrepeat_cfg.coverage_threshold)
            .count() as f64
            / ids.len() as f64
    };
    let subrepeat = subrepeat_flag(med, iqr, p75, coverage_frac, phantom, subrepeat_cfg);

    RowStats {
        identity_med: Some(med),
        identity_iqr: Some(iqr),
        identity_p25: Some(p25),
        identity_n: ids.len(),
        shift_med,
        shift_consistency,
        phantom,
        subrepeat,
        coverage_frac: Some(coverage_frac),
    }
}

/// Apply the founder gate: for each record, identify the lowest-rank
/// "founder" row (identity_med ≥ `founder_id_min` AND `phantom != true`)
/// and override `subrepeat` to `Some(false)` on any row whose period is
/// not strictly below that founder's period. By construction a subrepeat
/// is shorter than the founder monomer; this catches the false-positives
/// we'd otherwise emit on long-period bimodal rows (kite harmonics that
/// happen to look bimodal across the array).
///
/// Returns the number of subrepeat rows that were overridden.
fn enforce_subrepeat_founder_gate(
    rows: &[io::RawRow],
    stats: &mut [RowStats],
    founder_id_min: f64,
) -> usize {
    // First pass: build case_id → (lowest_rank, period) of qualifying founder.
    let mut founders: AHashMap<&str, (usize, usize)> = AHashMap::new();
    for (row, stat) in rows.iter().zip(stats.iter()) {
        let Some(id_med) = stat.identity_med else {
            continue;
        };
        if id_med < founder_id_min || matches!(stat.phantom, Some(true)) {
            continue;
        }
        let entry = founders
            .entry(row.case_id.as_str())
            .or_insert((usize::MAX, 0));
        if row.rank < entry.0 {
            *entry = (row.rank, row.period);
        }
    }

    // Second pass: override subrepeat where the candidate period is not
    // strictly below the founder.
    let mut overridden = 0usize;
    for (row, stat) in rows.iter().zip(stats.iter_mut()) {
        if !matches!(stat.subrepeat, Some(true)) {
            continue;
        }
        if let Some(&(_, founder_p)) = founders.get(row.case_id.as_str()) {
            if row.period >= founder_p {
                stat.subrepeat = Some(false);
                overridden += 1;
            }
        }
    }
    overridden
}

/// Subrepeat heuristic: bimodal identity distribution + intermediate
/// coverage + not phantom + not a real period.
fn subrepeat_flag(
    med: f64,
    iqr: f64,
    p75: f64,
    coverage_frac: f64,
    phantom: Option<bool>,
    cfg: &SubrepeatConfig,
) -> Option<bool> {
    if matches!(phantom, Some(true)) {
        return Some(false);
    }
    Some(
        p75 >= cfg.p75_min
            && iqr >= cfg.iqr_min
            && med < cfg.med_max
            && coverage_frac >= cfg.cov_min
            && coverage_frac <= cfg.cov_max,
    )
}

/// Aggregate per-pair shifts into (shift_med, shift_consistency, phantom).
fn aggregate_shifts(
    per_pair: &[(f64, Option<i32>)],
    period: usize,
    cfg: &PhantomConfig,
) -> (Option<i32>, Option<f64>, Option<bool>) {
    let mut shifts: Vec<i32> = per_pair
        .iter()
        .filter_map(|(id, s)| if *id >= cfg.identity_min { *s } else { None })
        .collect();
    if shifts.len() < cfg.min_pairs || period == 0 {
        return (None, None, None);
    }
    shifts.sort_unstable();
    let med = shifts[shifts.len() / 2];
    let within_1 = shifts.iter().filter(|s| (**s - med).abs() <= 1).count();
    let consistency = within_1 as f64 / shifts.len() as f64;
    let frac = (med.unsigned_abs() as f64) / (period as f64);
    let phantom = frac > cfg.tol_frac && consistency >= cfg.consistency_min;
    (Some(med), Some(consistency), Some(phantom))
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
        assert_eq!(s.format_row(), "NA\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA");
    }

    // --- subrepeat heuristic -----------------------------------------------

    #[test]
    fn subrepeat_flag_classic_bimodal_fires() {
        // Wide IQR + moderate median + high p75 + intermediate coverage,
        // not phantom.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.55, 0.40, 0.95, 0.30, Some(false), &cfg);
        assert_eq!(r, Some(true));
    }

    #[test]
    fn subrepeat_flag_phantom_blocks_subrepeat() {
        // Same bimodal stats but the row is already a phantom: never flag.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.55, 0.40, 0.95, 0.30, Some(true), &cfg);
        assert_eq!(r, Some(false));
    }

    #[test]
    fn subrepeat_flag_real_period_does_not_fire() {
        // High median ⇒ real period; med_max gate blocks subrepeat.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.95, 0.05, 0.97, 0.95, Some(false), &cfg);
        assert_eq!(r, Some(false));
    }

    #[test]
    fn subrepeat_flag_noise_does_not_fire() {
        // Uniformly low identity, narrow IQR ⇒ noise; p75/iqr gates block.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.42, 0.07, 0.46, 0.0, Some(false), &cfg);
        assert_eq!(r, Some(false));
    }

    #[test]
    fn subrepeat_flag_high_iqr_low_p75_does_not_fire() {
        // Wide IQR but the top isn't high enough — p75_min gate blocks.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.40, 0.25, 0.62, 0.10, Some(false), &cfg);
        assert_eq!(r, Some(false));
    }

    #[test]
    fn subrepeat_flag_phantom_unknown_does_not_block() {
        // phantom = None ⇒ no phantom info, but the row passes the
        // bimodality + coverage gates: still fire.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.55, 0.40, 0.95, 0.30, None, &cfg);
        assert_eq!(r, Some(true));
    }

    #[test]
    fn subrepeat_flag_too_high_coverage_does_not_fire() {
        // Bimodality gates pass but coverage > cov_max ⇒ real period
        // (coverage_frac too high), not a subrepeat.
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.55, 0.40, 0.95, 0.80, Some(false), &cfg);
        assert_eq!(r, Some(false));
    }

    #[test]
    fn subrepeat_flag_too_low_coverage_does_not_fire() {
        // Bimodality gates pass but coverage < cov_min ⇒ too few hits
        // for a real subrepeat (noise-leaning).
        let cfg = SubrepeatConfig::default();
        let r = subrepeat_flag(0.55, 0.40, 0.95, 0.05, Some(false), &cfg);
        assert_eq!(r, Some(false));
    }

    // --- founder gate ------------------------------------------------------

    fn mk_row(case: &str, rank: usize, period: usize) -> io::RawRow {
        io::RawRow {
            line: String::new(),
            case_id: case.to_string(),
            rank,
            period,
        }
    }

    fn mk_stats(
        identity_med: Option<f64>,
        phantom: Option<bool>,
        subrepeat: Option<bool>,
    ) -> RowStats {
        RowStats {
            identity_med,
            identity_iqr: identity_med.map(|_| 0.1),
            identity_p25: identity_med.map(|v| v - 0.05),
            identity_n: if identity_med.is_some() { 200 } else { 0 },
            shift_med: None,
            shift_consistency: None,
            phantom,
            subrepeat,
            coverage_frac: identity_med.map(|_| 0.3),
        }
    }

    #[test]
    fn founder_gate_overrides_when_period_exceeds_founder() {
        // TRC_41-style: rank-1 P=470 is the real founder; a later
        // long-period row was flagged subrepeat=true by the bimodality
        // heuristic; the founder gate should knock it back to false.
        let rows = vec![mk_row("TRC_41", 1, 470), mk_row("TRC_41", 4, 3289)];
        let mut stats = vec![
            mk_stats(Some(0.94), Some(false), Some(false)),
            mk_stats(Some(0.66), Some(false), Some(true)),
        ];
        let n = enforce_subrepeat_founder_gate(&rows, &mut stats, 0.70);
        assert_eq!(n, 1);
        assert_eq!(stats[0].subrepeat, Some(false));
        assert_eq!(stats[1].subrepeat, Some(false)); // overridden
    }

    #[test]
    fn founder_gate_keeps_short_subrepeat() {
        // TRC_104-style: rank-1 P=180 is the founder; rank-2 P=36 is a
        // genuine subrepeat (period < founder). Flag stays true.
        let rows = vec![mk_row("TRC_104", 1, 180), mk_row("TRC_104", 2, 36)];
        let mut stats = vec![
            mk_stats(Some(0.86), Some(false), Some(false)),
            mk_stats(Some(0.60), Some(false), Some(true)),
        ];
        let n = enforce_subrepeat_founder_gate(&rows, &mut stats, 0.70);
        assert_eq!(n, 0);
        assert_eq!(stats[1].subrepeat, Some(true));
    }

    #[test]
    fn founder_gate_skips_when_no_qualifying_founder() {
        // No row has identity_med ≥ founder_id_min ⇒ no founder
        // identified ⇒ subrepeat untouched.
        let rows = vec![mk_row("CASE_X", 1, 500), mk_row("CASE_X", 2, 40)];
        let mut stats = vec![
            mk_stats(Some(0.55), Some(false), Some(false)),
            mk_stats(Some(0.55), Some(false), Some(true)),
        ];
        let n = enforce_subrepeat_founder_gate(&rows, &mut stats, 0.70);
        assert_eq!(n, 0);
        assert_eq!(stats[1].subrepeat, Some(true));
    }

    #[test]
    fn founder_gate_ignores_phantom_when_picking_founder() {
        // Rank 1 is high-id but phantom-flagged (a sub-period harmonic
        // that fooled kite); rank 2 is the true founder. The candidate
        // at rank 3 has period between the phantom's and the founder's
        // — should be gated against rank 2, not rank 1.
        let rows = vec![
            mk_row("CASE_Y", 1, 100), // phantom — ignored as founder
            mk_row("CASE_Y", 2, 470), // true founder
            mk_row("CASE_Y", 3, 600), // candidate; would-be subrepeat
        ];
        let mut stats = vec![
            mk_stats(Some(0.90), Some(true), Some(false)),
            mk_stats(Some(0.90), Some(false), Some(false)),
            mk_stats(Some(0.55), Some(false), Some(true)),
        ];
        let n = enforce_subrepeat_founder_gate(&rows, &mut stats, 0.70);
        // 600 ≥ 470 (founder), so override fires.
        assert_eq!(n, 1);
        assert_eq!(stats[2].subrepeat, Some(false));
    }

    #[test]
    fn founder_gate_uses_lowest_qualifying_rank_not_earliest_iter() {
        // Builder order is rank 5 then rank 2; the founder must be the
        // rank-2 row regardless of iteration order.
        let rows = vec![
            mk_row("CASE_Z", 5, 800),
            mk_row("CASE_Z", 2, 200),
            mk_row("CASE_Z", 7, 500),
        ];
        let mut stats = vec![
            mk_stats(Some(0.90), Some(false), Some(false)),
            mk_stats(Some(0.90), Some(false), Some(false)),
            mk_stats(Some(0.55), Some(false), Some(true)),
        ];
        let n = enforce_subrepeat_founder_gate(&rows, &mut stats, 0.70);
        // Founder is rank-2 with period 200. Candidate is rank-7 at 500.
        // 500 ≥ 200 ⇒ override.
        assert_eq!(n, 1);
        assert_eq!(stats[2].subrepeat, Some(false));
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
            &PhantomConfig::default(),
            &SubrepeatConfig::default(),
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
            &PhantomConfig::default(),
            &SubrepeatConfig::default(),
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
    fn resolved_band_short_period_uses_slop_floor() {
        let cfg = Config {
            slop: 5,
            band: 0,
            ..Config::default()
        };
        // 2·5 = 10, max(20, 10) = 20; ⌈0.02·100⌉ = 2 → still 20.
        assert_eq!(cfg.resolved_band(100), 20);
        assert_eq!(cfg.resolved_band(500), 20); // ⌈10⌉ < 20
    }

    #[test]
    fn resolved_band_wider_slop_dominates() {
        let cfg = Config {
            slop: 50,
            band: 0,
            ..Config::default()
        };
        // 2·50 = 100; ⌈0.02·100⌉ = 2 → 100 wins.
        assert_eq!(cfg.resolved_band(100), 100);
    }

    #[test]
    fn resolved_band_long_period_scales_with_period() {
        let cfg = Config {
            slop: 10,
            band: 0,
            ..Config::default()
        };
        // slop_floor = 20; ⌈0.02·2870⌉ = 58 → 58.
        assert_eq!(cfg.resolved_band(2870), 58);
        // ⌈0.02·5000⌉ = 100.
        assert_eq!(cfg.resolved_band(5000), 100);
        // ⌈0.02·1000⌉ = 20 → unchanged.
        assert_eq!(cfg.resolved_band(1000), 20);
        // ⌈0.02·1050⌉ = 21 → just over.
        assert_eq!(cfg.resolved_band(1050), 21);
    }

    #[test]
    fn resolved_band_user_override_bypasses_formula() {
        let cfg = Config {
            slop: 10,
            band: 7,
            ..Config::default()
        };
        assert_eq!(cfg.resolved_band(100), 7);
        assert_eq!(cfg.resolved_band(2870), 7);
        assert_eq!(cfg.resolved_band(5000), 7);
    }

    #[test]
    fn long_monomer_with_internal_indels_recovers_with_auto_band() {
        // Regression test for the band=20 saturation artifact on long
        // monomers (the TRC_463 case). Construct 8 copies of a ~2000 bp
        // template; each copy gets a 40 bp insertion at a tile-specific
        // random position, so adjacent tiles share ~98 % of their
        // content but the optimal DP path between them drifts ~40 cells
        // off the diagonal at the location of either insertion.
        //
        // Period = 2040 (template + insertion), so:
        //   - band = 20 ⇒ DP saturates at the 40-cell drift, identity
        //     collapses to roughly the band-cap floor (≈ 0.5)
        //   - band = auto = ⌈0.02 · 2040⌉ = 41 ⇒ DP follows the path,
        //     identity recovers above 0.95
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let mut rng = StdRng::seed_from_u64(11);
        let bases = b"ACGT";
        let template_len = 2000usize;
        let insertion_len = 40usize;
        let period = template_len + insertion_len; // 2040

        let template: Vec<u8> = (0..template_len)
            .map(|_| bases[rng.random_range(0..4)])
            .collect();

        // 8 tiles, each = template + insertion at a tile-specific position.
        let mut seq: Vec<u8> = Vec::new();
        for _ in 0..8 {
            let insert_pos = rng.random_range(200..(template_len - 200));
            let insertion: Vec<u8> = (0..insertion_len)
                .map(|_| bases[rng.random_range(0..4)])
                .collect();
            let mut tile = template.clone();
            tile.splice(insert_pos..insert_pos, insertion);
            seq.extend_from_slice(&tile);
        }

        let case_id = "synthetic_long_indel";
        let cfg_auto = Config::default();
        let band_auto = cfg_auto.resolved_band(period);
        // auto band = max(20, 20, ⌈40.8⌉) = 41
        assert!(
            band_auto >= 41,
            "expected auto band ≥ 41, got {}",
            band_auto
        );

        let cfg_tight = Config {
            band: 20,
            ..Config::default()
        };

        let sample_cfg = cfg_auto.sample_cfg();
        let scoring = aligner::ScoringConfig::default();
        let phantom_cfg = PhantomConfig::default();
        let subrepeat_cfg = SubrepeatConfig::default();
        let mut scratch = aligner::Scratch::new();

        let stats_auto = rescore_one(
            &seq,
            period,
            case_id,
            &sample_cfg,
            band_auto,
            &scoring,
            &phantom_cfg,
            &subrepeat_cfg,
            &mut scratch,
        );
        let stats_tight = rescore_one(
            &seq,
            period,
            case_id,
            &sample_cfg,
            cfg_tight.resolved_band(period),
            &scoring,
            &phantom_cfg,
            &subrepeat_cfg,
            &mut scratch,
        );

        let id_auto = stats_auto.identity_med.unwrap();
        let id_tight = stats_tight.identity_med.unwrap();
        assert!(
            id_auto > 0.90,
            "auto-band identity should recover, got {}",
            id_auto
        );
        // The gap on synthetic data (8 pp here) is smaller than what
        // shows up on real long-monomer satellites like TRC_463
        // (~40 pp), because the tight band still finds a suboptimal
        // path through forced gaps on simple inputs. The direction is
        // what matters: auto > tight by a measurable margin.
        assert!(
            id_auto - id_tight >= 0.05,
            "expected ≥ 5 pp gap between auto={} and tight={}",
            id_auto,
            id_tight
        );
    }

    // --- shift aggregation -------------------------------------------------

    fn pp(id: f64, shift: Option<i32>) -> (f64, Option<i32>) {
        (id, shift)
    }

    #[test]
    fn aggregate_shifts_concentrated_nonzero_flags_phantom() {
        // Mirror the TRC_755 P=56 case: most pairs at +6 shift, some at +5,
        // a few elsewhere; all have high identity.
        let mut per_pair = vec![];
        for _ in 0..150 {
            per_pair.push(pp(0.85, Some(6)));
        }
        for _ in 0..40 {
            per_pair.push(pp(0.85, Some(5)));
        }
        for _ in 0..10 {
            per_pair.push(pp(0.85, Some(0)));
        }
        let (shift, cons, phantom) = aggregate_shifts(&per_pair, 56, &PhantomConfig::default());
        assert_eq!(shift, Some(6));
        let cons = cons.unwrap();
        assert!(cons > 0.9, "expected high consistency, got {}", cons);
        assert_eq!(phantom, Some(true));
    }

    #[test]
    fn aggregate_shifts_concentrated_zero_does_not_flag() {
        // Real period: shift is sharply at 0.
        let mut per_pair = vec![];
        for _ in 0..180 {
            per_pair.push(pp(0.95, Some(0)));
        }
        for _ in 0..20 {
            per_pair.push(pp(0.95, Some(2)));
        }
        let (shift, cons, phantom) = aggregate_shifts(&per_pair, 62, &PhantomConfig::default());
        assert_eq!(shift, Some(0));
        assert!(cons.unwrap() > 0.5);
        assert_eq!(phantom, Some(false));
    }

    #[test]
    fn aggregate_shifts_scattered_does_not_flag() {
        // Non-concentrated shifts ⇒ phantom off even if median is nonzero.
        let mut per_pair = vec![];
        for s in -8..=8i32 {
            for _ in 0..12 {
                per_pair.push(pp(0.7, Some(s)));
            }
        }
        let (shift, cons, phantom) = aggregate_shifts(&per_pair, 50, &PhantomConfig::default());
        assert!(shift.is_some());
        let cons = cons.unwrap();
        assert!(cons < 0.5, "expected low consistency, got {}", cons);
        assert_eq!(phantom, Some(false));
    }

    #[test]
    fn aggregate_shifts_excludes_low_identity_pairs() {
        // Most pairs are below threshold; only a handful contribute to
        // the shift. Below the min_pairs floor → NA.
        let mut per_pair = vec![];
        for _ in 0..195 {
            per_pair.push(pp(0.2, Some(6))); // would have flagged but below threshold
        }
        for _ in 0..3 {
            per_pair.push(pp(0.8, Some(0))); // only 3 pass identity_min
        }
        let (shift, cons, phantom) = aggregate_shifts(&per_pair, 56, &PhantomConfig::default());
        assert_eq!(shift, None);
        assert_eq!(cons, None);
        assert_eq!(phantom, None);
    }

    #[test]
    fn aggregate_shifts_relative_tol_does_not_flag_large_period() {
        // |shift|/period below threshold ⇒ no phantom even when concentrated.
        let mut per_pair = vec![];
        for _ in 0..200 {
            per_pair.push(pp(0.9, Some(6))); // |6|/1000 = 0.006 < 0.05
        }
        let (shift, _, phantom) = aggregate_shifts(&per_pair, 1000, &PhantomConfig::default());
        assert_eq!(shift, Some(6));
        assert_eq!(phantom, Some(false));
    }
}

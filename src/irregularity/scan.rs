//! Core irregularity scan — distance-residual + phase-bin clustering
//! algorithm. Mirrors `tools/rule_proto/irregularity_v2.py` v2.

use ahash::AHashMap;

// ---------- defaults (match the v2 Python prototype exactly) ----------

const K: usize = 6;
const TOP_KMERS: usize = 100;
const TOL_DK: f64 = 0.05;
const NEAR_DK_TOL: f64 = 0.05;
const MIN_NEAR_DK_FRAC: f64 = 0.6;
const MAX_N_SKIPPED: usize = 5;
const MIN_KMER_GROUPS: usize = 3;
const MIN_COPIES_FOR_SCAN: usize = 10;
const NOISE_FLOOR_MIN_BP: f64 = 1.0;
const NOISE_FLOOR_K: f64 = 3.0;
const SAME_SIGN_FRAC_MIN: f64 = 0.7;
const MIN_SUPPORT_FRAC: f64 = 0.5;
const EVENT_MERGE_GAP: i64 = 1;
const MAX_COPY_INDICES: usize = 50_000;

const DK_MULTIPLIERS: &[i64] = &[1, 2, 3];
const DK_DIVISORS: &[i64] = &[1, 2, 3, 4, 5];

/// Configuration for the scan. Defaults match the Python prototype.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// K-mer length (must match kite for register-marker consistency).
    pub k: usize,
    /// Top-N most frequent k-mers per record to consider.
    pub top_kmers: usize,
    /// Minimum array_length / period ratio. Below this → `too_short`.
    pub min_copies_for_scan: usize,
    /// Fraction of structural period below which an indel call is not
    /// "meaningful" (step magnitude floor). Default 0.05.
    pub step_min_frac_of_p: f64,
    /// Minimum independent k-mer groups required to enable analysis.
    pub min_kmer_groups: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            k: K,
            top_kmers: TOP_KMERS,
            min_copies_for_scan: MIN_COPIES_FOR_SCAN,
            step_min_frac_of_p: 0.05,
            min_kmer_groups: MIN_KMER_GROUPS,
        }
    }
}

/// Per-record output row (14 columns).
#[derive(Debug, Clone)]
pub struct RecordResult {
    pub record_id: String,
    pub array_length: usize,
    pub period_p: Option<f64>,
    pub n_kmer_groups: usize,
    pub n_pairs_total: usize,
    pub baseline_jitter_bp: Option<f64>,
    pub indel_event_count: Option<usize>,
    pub indel_burden_pct: Option<f64>,
    pub indel_max_shift_bp: Option<f64>,
    pub indel_drift_bp_per_kb: Option<f64>,
    pub dropout_event_count: Option<usize>,
    pub dropout_rate_per_pair: Option<f64>,
    pub flag: &'static str,
    pub notes: String,
}

impl RecordResult {
    fn flagged(
        record_id: &str,
        array_length: usize,
        period_p: Option<f64>,
        n_kmer_groups: usize,
        flag: &'static str,
        notes: String,
    ) -> Self {
        Self {
            record_id: record_id.to_string(),
            array_length,
            period_p,
            n_kmer_groups,
            n_pairs_total: 0,
            baseline_jitter_bp: None,
            indel_event_count: None,
            indel_burden_pct: None,
            indel_max_shift_bp: None,
            indel_drift_bp_per_kb: None,
            dropout_event_count: None,
            dropout_rate_per_pair: None,
            flag,
            notes,
        }
    }
}

// ---------- k-mer indexing ----------

/// Extract top-N frequent k-mers and their (sorted) positions. ACGT
/// only; non-ACGT k-mers skipped. Matches the Python prototype.
pub fn kmer_positions(seq: &[u8], k: usize, top: usize) -> Vec<(Vec<u8>, Vec<usize>)> {
    let l = seq.len();
    if l < k {
        return Vec::new();
    }
    let mut counts: AHashMap<Vec<u8>, usize> = AHashMap::new();
    let mut positions: AHashMap<Vec<u8>, Vec<usize>> = AHashMap::new();
    for i in 0..=(l - k) {
        let kmer = &seq[i..i + k];
        if !kmer.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let key = kmer.to_vec();
        *counts.entry(key.clone()).or_insert(0) += 1;
        positions.entry(key).or_default().push(i);
    }
    let mut by_count: Vec<(Vec<u8>, usize)> = counts.into_iter().collect();
    by_count.sort_by_key(|p| std::cmp::Reverse(p.1));
    by_count.truncate(top);
    by_count
        .into_iter()
        .map(|(km, _)| {
            let pos = positions.remove(&km).unwrap_or_default();
            (km, pos)
        })
        .collect()
}

// ---------- per-k-mer modal-distance + filter ----------

/// Mode of consecutive intervals (1-bp bins).
fn modal_distance(positions: &[usize]) -> Option<i64> {
    if positions.len() < 3 {
        return None;
    }
    let mut hist: AHashMap<i64, usize> = AHashMap::new();
    for i in 0..positions.len() - 1 {
        let d = (positions[i + 1] as i64) - (positions[i] as i64);
        *hist.entry(d).or_insert(0) += 1;
    }
    hist.into_iter().max_by_key(|kv| kv.1).map(|kv| kv.0)
}

/// True iff d_k is approximately m·P (m ∈ DK_MULTIPLIERS) or P/m
/// (m ∈ DK_DIVISORS) within `tol` relative.
fn compatible_with_p(d_k: f64, period: f64, tol: f64) -> bool {
    if d_k <= 0.0 || period <= 0.0 {
        return false;
    }
    for &m in DK_MULTIPLIERS {
        let target = (m as f64) * period;
        if (d_k - target).abs() / target <= tol {
            return true;
        }
    }
    for &m in DK_DIVISORS {
        let target = period / (m as f64);
        if (d_k - target).abs() / target <= tol {
            return true;
        }
    }
    false
}

/// Filter result: keep iff modal d_k compatible with P AND high
/// fraction of consecutive intervals are near n·d_k for small n.
/// Returns the d_k regardless of acceptance (caller may inspect).
pub fn select_register_locked(positions: &[usize], period: f64) -> (f64, bool) {
    let md = match modal_distance(positions) {
        Some(d) => d,
        None => return (0.0, false),
    };
    let d_k = md as f64;
    if !compatible_with_p(d_k, period, TOL_DK) {
        return (d_k, false);
    }
    let n = positions.len();
    if n < 3 {
        return (d_k, false);
    }
    let mut near = 0usize;
    let total = n - 1;
    for i in 0..total {
        let d = (positions[i + 1] as i64 - positions[i] as i64) as f64;
        let m_real = if d_k > 0.0 { d / d_k } else { 0.0 };
        let m_round = m_real.round() as i64;
        if m_round >= 1 && (m_round as usize) <= MAX_N_SKIPPED + 1 {
            let target = (m_round as f64) * d_k;
            if (d - target).abs() / target <= NEAR_DK_TOL {
                near += 1;
            }
        }
    }
    let frac = near as f64 / total as f64;
    (d_k, frac >= MIN_NEAR_DK_FRAC)
}

// ---------- phase-bin clustering ----------

/// Group k-mers into fixed-width phase bins on `phase mod P`.
/// Returns the list of bins (each bin = list of k-mer indices into
/// `selected`). K-mers in the same bin share monomer-internal position
/// and are NOT independent observations.
pub fn cluster_by_phase(
    selected: &[(Vec<u8>, Vec<usize>, f64)],
    period: f64,
    bin_width: usize,
) -> Vec<Vec<usize>> {
    if period <= 0.0 || selected.is_empty() {
        return Vec::new();
    }
    let p_int = period.round() as usize;
    if p_int == 0 {
        return Vec::new();
    }
    let bw = bin_width.max(1);
    let mut bins: AHashMap<usize, Vec<usize>> = AHashMap::new();
    for (idx, (_, pos, _)) in selected.iter().enumerate() {
        if pos.is_empty() {
            continue;
        }
        // Modal phase mod P (1-bp bins).
        let mut ph_counts: AHashMap<usize, usize> = AHashMap::new();
        for &x in pos {
            *ph_counts.entry(x % p_int).or_insert(0) += 1;
        }
        let modal_phase = ph_counts
            .into_iter()
            .max_by_key(|kv| kv.1)
            .map(|kv| kv.0)
            .unwrap_or(0);
        bins.entry(modal_phase / bw).or_default().push(idx);
    }
    bins.into_values().collect()
}

// ---------- per-pair residuals ----------

/// A single consecutive-pair observation: midpoint, residual,
/// n_skipped (= round(d/d_k) − 1, ≥ 0).
#[derive(Debug, Clone, Copy)]
pub struct PairObs {
    pub x_mid: usize,
    pub residual: f64,
    pub n_skipped: usize,
}

/// Per-pair residuals for one k-mer, using its own d_k.
pub fn per_pair_residuals(positions: &[usize], d_k: f64) -> Vec<PairObs> {
    if d_k <= 0.0 || positions.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(positions.len() - 1);
    for i in 0..positions.len() - 1 {
        let x_a = positions[i];
        let x_b = positions[i + 1];
        let d = (x_b as i64 - x_a as i64) as f64;
        let n = (d / d_k).round() as i64;
        if n < 1 || n as usize > MAX_N_SKIPPED + 1 {
            continue;
        }
        let residual = d - (n as f64) * d_k;
        out.push(PairObs {
            x_mid: (x_a + x_b) / 2,
            residual,
            n_skipped: (n - 1) as usize,
        });
    }
    out
}

// ---------- per-bin aggregation ----------

/// Aggregated stats for one genomic bin (one monomer-wide).
#[derive(Debug, Clone)]
pub struct BinAggregate {
    pub bin_idx: i64,
    pub n_groups_present: usize,
    pub consensus_residual: f64,
    pub same_sign_count: usize,
    pub n_dropouts: usize,
}

/// Aggregate per group → per bin. Each group contributes at most one
/// residual per bin (its median, in case multiple pairs in same bin).
pub fn aggregate_by_bin(group_obs: &[Vec<PairObs>], bin_size: f64) -> Vec<BinAggregate> {
    if bin_size <= 0.0 {
        return Vec::new();
    }
    // bin → (group_idx → list of residuals)
    let mut by_bin_group: AHashMap<i64, AHashMap<usize, Vec<f64>>> = AHashMap::new();
    let mut by_bin_dropout: AHashMap<i64, usize> = AHashMap::new();
    for (g_idx, obs_list) in group_obs.iter().enumerate() {
        for o in obs_list {
            let b = (o.x_mid as f64 / bin_size).floor() as i64;
            by_bin_group
                .entry(b)
                .or_default()
                .entry(g_idx)
                .or_default()
                .push(o.residual);
            *by_bin_dropout.entry(b).or_insert(0) += o.n_skipped;
        }
    }
    if by_bin_group.is_empty() {
        return Vec::new();
    }
    let mut bins: Vec<i64> = by_bin_group.keys().copied().collect();
    bins.sort();

    let mut out = Vec::with_capacity(bins.len());
    for b in bins {
        let groups_at_bin = &by_bin_group[&b];
        let mut per_group_medians: Vec<f64> = groups_at_bin
            .values()
            .map(|vs| {
                let mut sorted = vs.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                median_sorted(&sorted)
            })
            .collect();
        per_group_medians.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let cons = median_sorted(&per_group_medians);
        let sign_ref = if cons > 0.0 {
            1
        } else if cons < 0.0 {
            -1
        } else {
            0
        };
        let mut same_sign = 0usize;
        if sign_ref != 0 {
            for &v in &per_group_medians {
                if (v > 0.0 && sign_ref > 0) || (v < 0.0 && sign_ref < 0) {
                    same_sign += 1;
                }
            }
        }
        out.push(BinAggregate {
            bin_idx: b,
            n_groups_present: per_group_medians.len(),
            consensus_residual: cons,
            same_sign_count: same_sign,
            n_dropouts: *by_bin_dropout.get(&b).unwrap_or(&0),
        });
    }
    out
}

fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

// ---------- event detection ----------

#[derive(Debug, Clone)]
pub struct Event {
    pub bin_start: i64,
    pub bin_end: i64,
    pub consensus_shift: f64,
    pub max_support: usize,
    pub same_sign_max_frac: f64,
}

/// Bin-merge event detection: walk the aggregates, accept bins that
/// clear `step_floor` + support + same-sign-fraction gates, merging
/// adjacent same-sign passing bins (with at most `EVENT_MERGE_GAP`
/// non-passing bins between).
pub fn detect_events(aggs: &[BinAggregate], n_groups_total: usize, step_floor: f64) -> Vec<Event> {
    let min_support = std::cmp::max(
        3,
        (MIN_SUPPORT_FRAC * n_groups_total as f64).ceil() as usize,
    );
    let n = aggs.len();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < n {
        let a = &aggs[i];
        let passes = a.consensus_residual.abs() >= step_floor
            && a.n_groups_present >= min_support
            && (a.same_sign_count as f64 / a.n_groups_present.max(1) as f64) >= SAME_SIGN_FRAC_MIN;
        if !passes {
            i += 1;
            continue;
        }
        let sign: i64 = if a.consensus_residual > 0.0 { 1 } else { -1 };
        let mut last_passing = i;
        let mut j = i;
        while j + 1 < n && aggs[j + 1].bin_idx - aggs[last_passing].bin_idx <= EVENT_MERGE_GAP + 1 {
            let a2 = &aggs[j + 1];
            let s2: i64 = if a2.consensus_residual > 0.0 {
                1
            } else if a2.consensus_residual < 0.0 {
                -1
            } else {
                0
            };
            if s2 != sign {
                break;
            }
            let a2_passes = a2.consensus_residual.abs() >= step_floor
                && a2.n_groups_present >= min_support
                && (a2.same_sign_count as f64 / a2.n_groups_present.max(1) as f64)
                    >= SAME_SIGN_FRAC_MIN;
            if a2_passes {
                last_passing = j + 1;
            }
            j += 1;
        }
        let mut cons_vals: Vec<f64> = (i..=last_passing)
            .map(|k| aggs[k].consensus_residual)
            .filter(|v| v.abs() >= step_floor)
            .collect();
        cons_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let cons = if cons_vals.is_empty() {
            0.0
        } else {
            median_sorted(&cons_vals)
        };
        let max_support: usize = (i..=last_passing)
            .map(|k| aggs[k].n_groups_present)
            .max()
            .unwrap_or(0);
        let max_ss_frac: f64 = (i..=last_passing)
            .map(|k| aggs[k].same_sign_count as f64 / aggs[k].n_groups_present.max(1) as f64)
            .fold(0.0_f64, f64::max);
        out.push(Event {
            bin_start: aggs[i].bin_idx,
            bin_end: aggs[last_passing].bin_idx,
            consensus_shift: cons,
            max_support,
            same_sign_max_frac: max_ss_frac,
        });
        i = last_passing + 1;
    }
    out
}

// ---------- baseline noise floor ----------

/// Adaptive baseline jitter from MAD of all per-pair residuals.
/// `max(NOISE_FLOOR_MIN_BP, 1.4826 · MAD · NOISE_FLOOR_K)`.
pub fn adaptive_baseline(residuals: &[f64]) -> f64 {
    if residuals.len() < 2 {
        return NOISE_FLOOR_MIN_BP;
    }
    let mut sorted = residuals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = median_sorted(&sorted);
    let mut deviations: Vec<f64> = sorted.iter().map(|v| (*v - med).abs()).collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = median_sorted(&deviations);
    let sigma = 1.4826 * mad * NOISE_FLOOR_K;
    sigma.max(NOISE_FLOOR_MIN_BP)
}

// ---------- top-level per-record driver ----------

/// Analyse one record. `period` is the kite-supplied top period
/// (`monomer_size` from kite.tsv). `None` → flag = `no_period`.
pub fn analyse_record(
    record_id: &str,
    seq: &[u8],
    period: Option<f64>,
    cfg: &Config,
) -> RecordResult {
    let l = seq.len();
    let p = match period {
        Some(p) if p > 0.0 => p,
        _ => {
            return RecordResult::flagged(
                record_id,
                l,
                None,
                0,
                "no_period",
                "no kite period".into(),
            );
        }
    };
    if (l as f64) < cfg.min_copies_for_scan as f64 * p {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            0,
            "too_short",
            format!(
                "L={} < {}*P={:.0}",
                l,
                cfg.min_copies_for_scan,
                cfg.min_copies_for_scan as f64 * p
            ),
        );
    }

    // 1. Top-N k-mer positions.
    let pos_by_kmer = kmer_positions(seq, cfg.k, cfg.top_kmers);
    // 2. Per-k-mer register-lock filter.
    let mut selected: Vec<(Vec<u8>, Vec<usize>, f64)> = Vec::new();
    for (km, pos) in &pos_by_kmer {
        let (d_k, ok) = select_register_locked(pos, p);
        if ok {
            selected.push((km.clone(), pos.clone(), d_k));
        }
    }
    if selected.len() < cfg.min_kmer_groups {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            selected.len(),
            "no_register_lock",
            format!(
                "only {} register-locked k-mers at P={:.0}",
                selected.len(),
                p
            ),
        );
    }
    // 3. Phase-bin clustering.
    let bins = cluster_by_phase(&selected, p, cfg.k);
    if bins.len() < cfg.min_kmer_groups {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            bins.len(),
            "no_register_lock",
            format!("only {} independent phase groups", bins.len()),
        );
    }
    // 4. Per-group: each k-mer keeps its own d_k; pool the per-pair
    //    observations across all k-mers in the group.
    let mut group_obs: Vec<Vec<PairObs>> = Vec::new();
    for grp in &bins {
        let mut pooled: Vec<PairObs> = Vec::new();
        for &km_idx in grp {
            let (_, pos, d_k) = &selected[km_idx];
            if pos.len() < 3 || *d_k <= 0.0 {
                continue;
            }
            pooled.extend(per_pair_residuals(pos, *d_k));
        }
        if !pooled.is_empty() {
            group_obs.push(pooled);
        }
    }
    if group_obs.len() < cfg.min_kmer_groups {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            group_obs.len(),
            "no_register_lock",
            format!("after grouping, only {} usable groups", group_obs.len()),
        );
    }
    let n_pairs_total: usize = group_obs.iter().map(|g| g.len()).sum();

    // 5. Aggregate by P-wide bin.
    let aggs = aggregate_by_bin(&group_obs, p);
    if aggs.is_empty() {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            group_obs.len(),
            "no_register_lock",
            "no per-bin aggregates".into(),
        );
    }
    if aggs.len() > MAX_COPY_INDICES {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            group_obs.len(),
            "too_long",
            format!(
                "{} bins > {} cap (P̂={:.1})",
                aggs.len(),
                MAX_COPY_INDICES,
                p
            ),
        );
    }

    // 6. Baseline jitter from all per-pair residuals across groups.
    let all_resid: Vec<f64> = group_obs.iter().flatten().map(|o| o.residual).collect();
    if all_resid.is_empty() {
        return RecordResult::flagged(
            record_id,
            l,
            Some(p),
            group_obs.len(),
            "no_register_lock",
            "no residuals".into(),
        );
    }
    let baseline = adaptive_baseline(&all_resid);
    let step_min_abs = cfg.step_min_frac_of_p * p;
    let step_floor = baseline.max(step_min_abs);

    // 7. Event detection.
    let events = detect_events(&aggs, group_obs.len(), step_floor);
    let n_events = events.len();
    let abs_shifts: Vec<f64> = events.iter().map(|e| e.consensus_shift.abs()).collect();
    let signed_shifts: Vec<f64> = events.iter().map(|e| e.consensus_shift).collect();
    let sum_abs: f64 = abs_shifts.iter().sum();
    let sum_signed: f64 = signed_shifts.iter().sum();
    let burden_pct = if l > 0 {
        sum_abs / (l as f64) * 100.0
    } else {
        0.0
    };
    let drift_per_kb = if l > 0 {
        sum_signed.abs() / (l as f64) * 1000.0
    } else {
        0.0
    };
    let max_shift = abs_shifts.iter().fold(0.0_f64, |a, &b| a.max(b));
    let total_drops: usize = aggs.iter().map(|a| a.n_dropouts).sum();
    let dropout_rate = total_drops as f64 / n_pairs_total.max(1) as f64;

    RecordResult {
        record_id: record_id.to_string(),
        array_length: l,
        period_p: Some(p),
        n_kmer_groups: group_obs.len(),
        n_pairs_total,
        baseline_jitter_bp: Some(baseline),
        indel_event_count: Some(n_events),
        indel_burden_pct: Some(burden_pct),
        indel_max_shift_bp: Some(max_shift),
        indel_drift_bp_per_kb: Some(drift_per_kb),
        dropout_event_count: Some(total_drops),
        dropout_rate_per_pair: Some(dropout_rate),
        flag: "ok",
        notes: String::new(),
    }
}

// ---------- unit tests ----------

#[cfg(test)]
mod tests {
    use super::*;

    fn pos_vec(spec: &[(usize, usize, usize)]) -> Vec<usize> {
        // (start, step, count)
        let mut out = Vec::new();
        for &(start, step, count) in spec {
            for i in 0..count {
                out.push(start + step * i);
            }
        }
        out.sort();
        out
    }

    #[test]
    fn modal_distance_basic() {
        let p = pos_vec(&[(0, 100, 5)]);
        assert_eq!(modal_distance(&p), Some(100));
    }

    #[test]
    fn modal_distance_with_skip() {
        // four steps of 100, one skip → 200
        let p = vec![0, 100, 200, 400, 500];
        assert_eq!(modal_distance(&p), Some(100));
    }

    #[test]
    fn compatible_with_p_basics() {
        assert!(compatible_with_p(100.0, 100.0, 0.05));
        assert!(compatible_with_p(200.0, 100.0, 0.05)); // 2*P
        assert!(compatible_with_p(50.0, 100.0, 0.05)); // P/2
        assert!(!compatible_with_p(73.0, 100.0, 0.05));
    }

    #[test]
    fn per_pair_residual_self_corrects_dropout() {
        // Distance = 2*P (one missing) → residual = 0, n_skipped = 1
        let positions = vec![0, 200, 300];
        let obs = per_pair_residuals(&positions, 100.0);
        assert_eq!(obs.len(), 2);
        assert!(obs[0].residual.abs() < 1e-9);
        assert_eq!(obs[0].n_skipped, 1); // 200 = 2*P, one skip
        assert!(obs[1].residual.abs() < 1e-9);
        assert_eq!(obs[1].n_skipped, 0); // 100 = 1*P
    }

    #[test]
    fn per_pair_residual_picks_up_indel() {
        // 100, 100, 113 (+13 indel), 100
        let positions = vec![0, 100, 200, 313, 413];
        let obs = per_pair_residuals(&positions, 100.0);
        assert_eq!(obs.len(), 4);
        assert!(obs[0].residual.abs() < 1e-9);
        assert!(obs[1].residual.abs() < 1e-9);
        assert!((obs[2].residual - 13.0).abs() < 1e-9);
        assert_eq!(obs[2].n_skipped, 0);
        assert!(obs[3].residual.abs() < 1e-9);
    }

    #[test]
    fn adaptive_baseline_floor() {
        // All zeros → MAD = 0 → floor at NOISE_FLOOR_MIN_BP.
        let r = vec![0.0; 50];
        assert!((adaptive_baseline(&r) - NOISE_FLOOR_MIN_BP).abs() < 1e-9);
    }

    #[test]
    fn adaptive_baseline_picks_up_spread() {
        // Symmetric ±5 → MAD = 5, σ = 1.4826*5*3 ≈ 22.2
        let r: Vec<f64> = (0..50)
            .map(|i| if i % 2 == 0 { 5.0 } else { -5.0 })
            .collect();
        let b = adaptive_baseline(&r);
        assert!(b > 15.0);
        assert!(b < 30.0);
    }

    #[test]
    fn detect_events_finds_single_spike() {
        // One bin with a +20 consensus, support 10, same-sign 9 → call it.
        let agg = BinAggregate {
            bin_idx: 5,
            n_groups_present: 10,
            consensus_residual: 20.0,
            same_sign_count: 9,
            n_dropouts: 0,
        };
        let events = detect_events(&[agg], 10, 5.0);
        assert_eq!(events.len(), 1);
        assert!((events[0].consensus_shift - 20.0).abs() < 1e-9);
    }

    #[test]
    fn detect_events_rejects_below_floor() {
        let agg = BinAggregate {
            bin_idx: 5,
            n_groups_present: 10,
            consensus_residual: 3.0,
            same_sign_count: 9,
            n_dropouts: 0,
        };
        let events = detect_events(&[agg], 10, 5.0);
        assert!(events.is_empty());
    }

    #[test]
    fn detect_events_rejects_insufficient_support() {
        let agg = BinAggregate {
            bin_idx: 5,
            n_groups_present: 2,
            consensus_residual: 20.0,
            same_sign_count: 2,
            n_dropouts: 0,
        };
        let events = detect_events(&[agg], 10, 5.0);
        assert!(events.is_empty());
    }

    #[test]
    fn cluster_by_phase_short_period() {
        // P=60, k=6 → 10 phase bins. K-mers at modal phases 0, 7, 14, 21 →
        // map to bins 0, 1, 2, 3 (one k-mer per bin).
        let selected: Vec<(Vec<u8>, Vec<usize>, f64)> = vec![
            (b"AAAAAA".to_vec(), vec![0, 60, 120], 60.0),
            (b"BBBBBB".to_vec(), vec![7, 67, 127], 60.0),
            (b"CCCCCC".to_vec(), vec![14, 74, 134], 60.0),
            (b"DDDDDD".to_vec(), vec![21, 81, 141], 60.0),
        ];
        let bins = cluster_by_phase(&selected, 60.0, K);
        assert_eq!(bins.len(), 4);
        for b in &bins {
            assert_eq!(b.len(), 1);
        }
    }

    #[test]
    fn cluster_by_phase_same_bin_for_close_phases() {
        // P=60, bin_width=6: phases 0 and 5 → same bin (0/6=0, 5/6=0).
        let selected: Vec<(Vec<u8>, Vec<usize>, f64)> = vec![
            (b"AAAAAA".to_vec(), vec![0, 60], 60.0),
            (b"BBBBBB".to_vec(), vec![5, 65], 60.0),
        ];
        let bins = cluster_by_phase(&selected, 60.0, 6);
        assert_eq!(bins.len(), 1);
        assert_eq!(bins[0].len(), 2);
    }

    #[test]
    fn no_period_flag() {
        let r = analyse_record("x", b"AAACGT", None, &Config::default());
        assert_eq!(r.flag, "no_period");
    }

    #[test]
    fn too_short_flag() {
        // P=100, L=500 < 10*P=1000 → too_short
        let r = analyse_record("x", &[b'A'; 500], Some(100.0), &Config::default());
        assert_eq!(r.flag, "too_short");
    }
}

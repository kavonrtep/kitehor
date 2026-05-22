//! Decision tree: Case A → Case B → fallbacks → unresolved.
//! Port of `tools/rule_proto/rule_proto.py::classify_case`.

use super::cluster::{cluster_peaks, Cluster};

#[derive(Debug, Clone, Copy)]
pub struct PeakRow {
    pub rank: u32,
    pub period: usize,
    pub score2_norm: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub tol: f64,
    pub min_period: usize,
    pub min_cluster_frac: f64,
    pub k_max: usize,
    pub non_mono_ratio: f64,
    pub founder_floor: f64,
    pub high_k_tile_floor: f64,
    pub lone_significant_frac: f64,
}

impl Default for Config {
    fn default() -> Self {
        // Prototype defaults from rule_proto.py.
        Self {
            tol: 0.015,
            min_period: 20,
            min_cluster_frac: 0.01,
            k_max: 30,
            non_mono_ratio: 0.5,
            founder_floor: 0.1,
            high_k_tile_floor: 0.05,
            lone_significant_frac: 0.1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerdictKind {
    Hor,
    SimpleTr,
    Unresolved,
}

impl VerdictKind {
    pub fn as_str(self) -> &'static str {
        match self {
            VerdictKind::Hor => "hor",
            VerdictKind::SimpleTr => "simple_tr",
            VerdictKind::Unresolved => "unresolved",
        }
    }
}

/// Verdict row matching the prototype's 10-column TSV schema verbatim.
/// `None` for `founder` / `tile` / `multiplicity` / etc. round-trips to
/// an empty cell when serialised (Python `None` convention).
#[derive(Debug, Clone)]
pub struct Verdict {
    pub case_id: String,
    pub kind: VerdictKind,
    pub founder: Option<f64>,
    pub multiplicity: Option<u32>,
    pub tile: Option<f64>,
    pub founder_score: Option<f64>,
    pub tile_score: Option<f64>,
    pub confidence: Option<f64>,
    pub n_clusters: u32,
    pub reason: String,
}

/// Convenience wrapper used by `classify(&KiteResult, &Config)`.
pub fn decide(case_id: &str, rows: &[PeakRow], cfg: &Config) -> Verdict {
    decide_with_clusters(case_id, rows, cfg).0
}

/// Apply the classifier; return the verdict and the post-filter cluster
/// list (used by `--dump-clusters`).
pub fn decide_with_clusters(
    case_id: &str,
    rows: &[PeakRow],
    cfg: &Config,
) -> (Verdict, Vec<Cluster>) {
    // 1. Drop peaks below the kmer-floor.
    let filtered: Vec<PeakRow> = rows
        .iter()
        .filter(|r| r.period >= cfg.min_period)
        .copied()
        .collect();

    // 2. Cluster.
    let mut clusters = cluster_peaks(&filtered, cfg.tol);

    // 3. Drop low-score clusters relative to max.
    if !clusters.is_empty() {
        let max_s = clusters
            .iter()
            .map(|c| c.total_score)
            .fold(f64::NEG_INFINITY, f64::max);
        clusters.retain(|c| c.total_score >= cfg.min_cluster_frac * max_s);
    }

    if clusters.is_empty() {
        return (
            Verdict {
                case_id: case_id.to_string(),
                kind: VerdictKind::Unresolved,
                founder: None,
                multiplicity: None,
                tile: None,
                founder_score: None,
                tile_score: None,
                confidence: None,
                n_clusters: 0,
                reason: "no_clusters".into(),
            },
            clusters,
        );
    }

    // 4. Sort by total_score desc; top is the strongest cluster.
    // `sort_by` with `partial_cmp` gives stable ordering on ties.
    let mut by_score: Vec<usize> = (0..clusters.len()).collect();
    by_score.sort_by(|&a, &b| {
        clusters[b]
            .total_score
            .partial_cmp(&clusters[a].total_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let top_idx = by_score[0];
    let top = clusters[top_idx].clone();

    // 5a. Case A — top = k × shorter cluster.
    let mut case_a: Vec<(usize, u32)> = Vec::new();
    for (i, c) in clusters.iter().enumerate() {
        if i == top_idx {
            continue;
        }
        if c.rep_period >= top.rep_period {
            continue;
        }
        let Some(k) = multiplicity_match(c.rep_period, top.rep_period, cfg.k_max, cfg.tol) else {
            continue;
        };
        if c.total_score < cfg.founder_floor * top.total_score {
            continue;
        }
        if !harmonic_confirms_hor(&clusters, c.rep_period, k, cfg.tol) {
            continue;
        }
        case_a.push((i, k));
    }
    if !case_a.is_empty() {
        // Pick smallest rep_period among qualifiers.
        let (founder_idx, k) = case_a
            .into_iter()
            .min_by(|a, b| {
                clusters[a.0]
                    .rep_period
                    .partial_cmp(&clusters[b.0].rep_period)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        let founder = clusters[founder_idx].clone();
        let v = hor_verdict(case_id, &founder, k, &top, "top_is_multiple_of_founder", &clusters);
        return (v, clusters);
    }

    // 5b. Case B — top is founder; walk k=2..k_max for a non-monotonic bump.
    let mut prev_seen_score = top.total_score;
    for k in 2..=(cfg.k_max as u32) {
        let expected = (k as f64) * top.rep_period;
        // matching = clusters near k * top.period, excluding top itself
        let mut matching: Vec<&Cluster> = Vec::new();
        for (i, c) in clusters.iter().enumerate() {
            if i == top_idx {
                continue;
            }
            let rel = (c.rep_period - expected).abs() / expected;
            if rel <= cfg.tol {
                matching.push(c);
            }
        }
        if matching.is_empty() {
            continue;
        }
        let best = matching
            .into_iter()
            .max_by(|a, b| {
                a.total_score
                    .partial_cmp(&b.total_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        let s_k = best.total_score;

        let (score_test, extra_tile_doubling, reason) = if k == 2 {
            (s_k > cfg.non_mono_ratio * top.total_score, true, "k2_ratio")
        } else {
            let st = s_k > prev_seen_score && s_k >= cfg.high_k_tile_floor * top.total_score;
            let etd =
                cluster_score_near(&clusters, 2.0 * (k as f64) * top.rep_period, cfg.tol) > 0.0;
            (st, etd, "non_monotonic_bump")
        };

        if score_test
            && extra_tile_doubling
            && harmonic_confirms_hor(&clusters, top.rep_period, k, cfg.tol)
        {
            let v = hor_verdict(case_id, &top, k, best, reason, &clusters);
            return (v, clusters);
        }
        prev_seen_score = s_k;
    }

    // 6. Has any larger multiple at all?
    let has_any = (2..=(cfg.k_max as u32)).any(|k| {
        let expected = (k as f64) * top.rep_period;
        clusters.iter().enumerate().any(|(i, c)| {
            if i == top_idx {
                return false;
            }
            let rel = (c.rep_period - expected).abs() / expected;
            rel <= cfg.tol
        })
    });
    let sum_s: f64 = clusters.iter().map(|c| c.total_score).sum();
    let total_s: f64 = if sum_s > 0.0 { sum_s } else { 1.0 };
    if has_any {
        let conf = (top.total_score / total_s).clamp(0.0, 1.0);
        return (
            simple_tr_verdict(case_id, &top, conf, clusters.len() as u32, "monotonic_multiples"),
            clusters,
        );
    }
    // 7. Lone-significant-cluster fallback.
    let n_significant = clusters
        .iter()
        .filter(|c| c.total_score >= cfg.lone_significant_frac * top.total_score)
        .count();
    if n_significant == 1 {
        let conf = (top.total_score / total_s).clamp(0.0, 1.0);
        return (
            simple_tr_verdict(case_id, &top, conf, clusters.len() as u32, "lone_significant_cluster"),
            clusters,
        );
    }

    // 8. Unresolved with founder.
    let conf = (top.total_score / total_s).clamp(0.0, 1.0);
    (
        Verdict {
            case_id: case_id.to_string(),
            kind: VerdictKind::Unresolved,
            founder: Some(round2(top.rep_period)),
            multiplicity: None,
            tile: None,
            founder_score: Some(top.total_score),
            tile_score: None,
            confidence: Some(conf),
            n_clusters: clusters.len() as u32,
            reason: "no_multiples".into(),
        },
        clusters,
    )
}

fn hor_verdict(
    case_id: &str,
    founder: &Cluster,
    k: u32,
    tile: &Cluster,
    reason: &str,
    all_clusters: &[Cluster],
) -> Verdict {
    let sum: f64 = all_clusters.iter().map(|c| c.total_score).sum();
    let total: f64 = if sum > 0.0 { sum } else { 1.0 };
    // Founder is `tile` only in the (degenerate) case Case-B-k=2 where
    // the loop matched top itself — that's excluded by the `i == top_idx`
    // skip. So in practice founder != tile here.
    let conf = ((founder.total_score + tile.total_score) / total).clamp(0.0, 1.0);
    Verdict {
        case_id: case_id.to_string(),
        kind: VerdictKind::Hor,
        founder: Some(round2(founder.rep_period)),
        multiplicity: Some(k),
        tile: Some(round2(tile.rep_period)),
        founder_score: Some(founder.total_score),
        tile_score: Some(tile.total_score),
        confidence: Some(conf),
        n_clusters: all_clusters.len() as u32,
        reason: reason.into(),
    }
}

fn simple_tr_verdict(
    case_id: &str,
    top: &Cluster,
    confidence: f64,
    n_clusters: u32,
    reason: &str,
) -> Verdict {
    let p = round2(top.rep_period);
    Verdict {
        case_id: case_id.to_string(),
        kind: VerdictKind::SimpleTr,
        founder: Some(p),
        multiplicity: Some(1),
        tile: Some(p),
        founder_score: Some(top.total_score),
        tile_score: Some(top.total_score),
        confidence: Some(confidence),
        n_clusters,
        reason: reason.into(),
    }
}

/// `k = round(p_long / p_short)`; reject if `k < 2 || k > k_max`; then
/// check `|p_long − k·p_short| / p_long <= tol`.
pub fn multiplicity_match(p_short: f64, p_long: f64, k_max: usize, tol: f64) -> Option<u32> {
    if p_short <= 0.0 || p_long <= 0.0 {
        return None;
    }
    let k_f = (p_long / p_short).round();
    if !(2.0..=(k_max as f64)).contains(&k_f) {
        return None;
    }
    let k = k_f as u32;
    let expected = (k as f64) * p_short;
    if (p_long - expected).abs() / p_long <= tol {
        Some(k)
    } else {
        None
    }
}

/// Return the highest `total_score` of any cluster within `tol` of
/// `target`, or 0.0 if none. `target <= 0` returns 0.0.
pub fn cluster_score_near(clusters: &[Cluster], target: f64, tol: f64) -> f64 {
    if target <= 0.0 {
        return 0.0;
    }
    let mut best = 0.0f64;
    for c in clusters {
        let rel = (c.rep_period - target).abs() / target;
        if rel <= tol && c.total_score > best {
            best = c.total_score;
        }
    }
    best
}

/// Harmonic-series check: a real HOR(founder=p, k) shows
/// `score(2k·p) ≥ score((k+1)·p)`. A simple TR (monotonic decay) shows
/// the opposite. Zero scores treated as 0 (not skipped); both missing
/// passes trivially.
pub fn harmonic_confirms_hor(
    clusters: &[Cluster],
    founder_period: f64,
    k: u32,
    tol: f64,
) -> bool {
    let kf = k as f64;
    let s_double_tile = cluster_score_near(clusters, 2.0 * kf * founder_period, tol);
    let s_off_mult = cluster_score_near(clusters, (kf + 1.0) * founder_period, tol);
    s_double_tile >= s_off_mult
}

/// Round to 2 decimals (matches Python `round(x, 2)`).
fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(rank: u32, period: usize, s: f64) -> PeakRow {
        PeakRow {
            rank,
            period,
            score2_norm: s,
        }
    }

    #[test]
    fn no_clusters_is_no_signal_verdict() {
        let cfg = Config::default();
        let v = decide("empty", &[], &cfg);
        assert_eq!(v.kind, VerdictKind::Unresolved);
        assert_eq!(v.reason, "no_clusters");
        assert_eq!(v.n_clusters, 0);
        assert!(v.founder.is_none());
    }

    #[test]
    fn clean_k3_hor_tile_dominant() {
        // Tile dominates (top by total_score), founder visible below.
        // Case A: top=900 = 3 × founder=300, founder.score (0.20+0.10=0.30)
        // ≥ 0.1 × top.score (0.50). Harmonic check passes since no
        // (k+1)×founder = 1200 cluster exists.
        let rows = vec![
            row(1, 900, 0.50),
            row(2, 298, 0.20),
            row(3, 300, 0.10),
        ];
        let cfg = Config::default();
        let v = decide("hor_k3", &rows, &cfg);
        assert_eq!(v.kind, VerdictKind::Hor);
        assert_eq!(v.multiplicity, Some(3));
        let f = v.founder.unwrap();
        assert!((f - 299.0).abs() <= 2.0, "founder = {f}");
        let t = v.tile.unwrap();
        assert!((t - 900.0).abs() <= 2.0, "tile = {t}");
    }

    #[test]
    fn monotonic_decay_is_simple_tr() {
        // top at 300, harmonics at 600/900 with monotonically decreasing scores.
        let rows = vec![
            row(1, 300, 0.50),
            row(2, 600, 0.30),
            row(3, 900, 0.15),
        ];
        let cfg = Config::default();
        let v = decide("tr_mono", &rows, &cfg);
        assert_eq!(v.kind, VerdictKind::SimpleTr);
        assert_eq!(v.reason, "monotonic_multiples");
    }

    #[test]
    fn lone_significant_falls_through_to_simple_tr() {
        // Lone strong peak, no multiples at all → lone_significant_cluster.
        let rows = vec![row(1, 400, 0.50)];
        let cfg = Config::default();
        let v = decide("tr_lone", &rows, &cfg);
        assert_eq!(v.kind, VerdictKind::SimpleTr);
        assert_eq!(v.reason, "lone_significant_cluster");
    }

    #[test]
    fn multiplicity_match_round_and_tolerance() {
        // 2 × 300 = 600 exact
        assert_eq!(multiplicity_match(300.0, 600.0, 30, 0.015), Some(2));
        // 3 × 300 = 900, drift to 906 = 0.66% off — within tol=0.015
        assert_eq!(multiplicity_match(300.0, 906.0, 30, 0.015), Some(3));
        // 3 × 300 = 900, drift to 920 → 2.2% off, rejected at tol=0.015
        assert_eq!(multiplicity_match(300.0, 920.0, 30, 0.015), None);
        // k=1 rejected
        assert_eq!(multiplicity_match(300.0, 300.0, 30, 0.015), None);
        // k>k_max rejected
        assert_eq!(multiplicity_match(10.0, 350.0, 30, 0.015), None);
    }
}

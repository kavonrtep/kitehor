//! HOR-on-top classification layer.
//!
//! Consumes a `KiteResult` (the kite.R port output) and decides
//! one of:
//!
//! - `Hor`        — peaks form an integer-multiple family at a founder.
//!   Outputs (founder, tile, multiplicity).
//! - `Tandem`     — single dominant periodicity, no HOR structure.
//! - `Unresolved` — peaks present but no clean classification
//!   (variable-length tandem, scattered peaks, or
//!   near-tied peaks with no family).
//! - `NoSignal`   — kite produced no peaks above the noise envelope.
//!
//! Pure classification — does not re-detect periodicities. See
//! `HOR_LAYER_PLAN.md` for the full design.

use crate::kite::KiteResult;

#[derive(Debug, Clone, Copy)]
pub struct HorCallConfig {
    /// Max HOR multiplicity considered.
    pub qmax: usize,
    /// Minimum peaks in the family for an HOR call (incl. the founder).
    pub min_family_size: usize,
    /// family_score / total_score must be ≥ this for HOR.
    pub min_family_share: f64,
    /// Top-N peaks examined for the jitter detector.
    pub n_jitter: usize,
    /// Relative tolerance for the jitter band (±jitter_tol × top1).
    pub jitter_tol: f64,
    /// Min peaks in band to trigger `Unresolved` via jitter.
    pub jitter_thr: usize,
    /// Absolute tolerance for `|p − k·m_f| ≤ tol_bp`.
    pub tol_bp: usize,
    /// OR relative tolerance for the same check.
    pub tol_rel: f64,
    /// Top1 / top2 score ratio above which we call `Tandem`.
    pub dominance: f64,
    /// Lowest founder period considered when enumerating divisors.
    pub lo_period: usize,
    /// Min founder score as a fraction of top1's score. Rejects spurious
    /// small-period kite peaks that happen to divide many larger peaks
    /// (the TRC_1__15071287 m_f=60 case).
    pub min_founder_top1_share: f64,
    /// Require the top-K peaks (by score) to all belong to the proposed
    /// family. Catches cases where a kite peak family exists but the
    /// dominant peaks don't fit — i.e., there's a competing periodicity.
    pub require_top_k_in_family: usize,
    /// Required `tile_score / founder_score` for HOR. Real HOR has
    /// tile (HOR-unit clones) at least comparable to founder; pure
    /// tandems have a much weaker peak at 2·monomer / 3·monomer than
    /// at the monomer itself. Synthetic: HOR median 0.55, null median
    /// 0.22 → threshold ~0.30-0.40 separates them.
    pub min_tile_founder_ratio: f64,
}

impl Default for HorCallConfig {
    fn default() -> Self {
        Self {
            qmax: 30,
            min_family_size: 3,
            min_family_share: 0.50,
            n_jitter: 8,
            jitter_tol: 0.15,
            jitter_thr: 4,
            tol_bp: 5,
            tol_rel: 0.02,
            dominance: 3.0,
            lo_period: 15,
            min_founder_top1_share: 0.50,
            require_top_k_in_family: 3,
            min_tile_founder_ratio: 0.15,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HorVerdict {
    Hor,
    Tandem,
    Unresolved,
    NoSignal,
}

impl HorVerdict {
    pub fn as_str(self) -> &'static str {
        match self {
            HorVerdict::Hor => "hor",
            HorVerdict::Tandem => "tandem",
            HorVerdict::Unresolved => "unresolved",
            HorVerdict::NoSignal => "no_signal",
        }
    }
}

#[derive(Debug, Clone)]
pub struct HorCallResult {
    pub verdict: HorVerdict,
    pub founder_bp: Option<usize>,
    pub multiplicity: Option<usize>,
    pub tile_bp: Option<usize>,
    pub family_score: f64,
    pub family_size: usize,
    pub jitter: usize,
    pub reason: String,
}

/// Top-level entry point.
pub fn classify(kite: &KiteResult, cfg: &HorCallConfig) -> HorCallResult {
    let peaks = &kite.peaks;
    if peaks.is_empty() {
        return HorCallResult {
            verdict: HorVerdict::NoSignal,
            founder_bp: None,
            multiplicity: None,
            tile_bp: None,
            family_score: 0.0,
            family_size: 0,
            jitter: 0,
            reason: "no_signal".into(),
        };
    }
    let top1 = &peaks[0];
    if peaks.len() == 1 {
        return HorCallResult {
            verdict: HorVerdict::Tandem,
            founder_bp: Some(top1.period),
            multiplicity: Some(1),
            tile_bp: Some(top1.period),
            family_score: top1.score,
            family_size: 1,
            jitter: 1,
            reason: format!(
                "tandem:monomer={},score={:.6};single_peak",
                top1.period, top1.score
            ),
        };
    }
    let total_score: f64 = peaks.iter().map(|p| p.score).sum();

    // §2c — fit an HOR family.
    let (
        best,
        family_score_best,
        family_size_best,
        best_tile_k,
        best_tile_period,
        best_tile_score,
        best_founder_score,
    ) = find_best_family(peaks, cfg);

    // §2b — jitter detector.
    let jitter = compute_jitter(peaks, cfg);

    // Decision:
    if let Some(m_f) = best {
        let tile_founder_ratio = if best_founder_score > 0.0 {
            best_tile_score / best_founder_score
        } else {
            0.0
        };
        if family_size_best >= cfg.min_family_size
            && (total_score == 0.0 || family_score_best / total_score >= cfg.min_family_share)
            && best_tile_score > 0.0
            && tile_founder_ratio >= cfg.min_tile_founder_ratio
        {
            let mult = best_tile_k.max(2);
            let tile = best_tile_period;
            // Tile-centered jitter: count peaks within ±jitter_tol of
            // the proposed tile. If many peaks cluster around the tile
            // (variable-length tandem signature), call Unresolved
            // instead of HOR.
            let tile_jitter = peaks_within_band(peaks, tile, cfg.jitter_tol);
            let fam_str = describe_family(peaks, m_f, cfg);
            if tile_jitter >= cfg.jitter_thr {
                return HorCallResult {
                    verdict: HorVerdict::Unresolved,
                    founder_bp: None,
                    multiplicity: None,
                    tile_bp: None,
                    family_score: family_score_best,
                    family_size: family_size_best,
                    jitter: tile_jitter,
                    reason: format!(
                        "unresolved:variable_length_tandem:tile_jitter={}/around_tile={}@±{:.0}%;\
                         family_fit_would_say:founder={},k={},tile={};peaks={{{}}}",
                        tile_jitter,
                        tile,
                        cfg.jitter_tol * 100.0,
                        m_f,
                        mult,
                        tile,
                        fam_str,
                    ),
                };
            }
            return HorCallResult {
                verdict: HorVerdict::Hor,
                founder_bp: Some(m_f),
                multiplicity: Some(mult),
                tile_bp: Some(tile),
                family_score: family_score_best,
                family_size: family_size_best,
                jitter: tile_jitter,
                reason: format!(
                    "hor:founder={},k={},tile={},family_size={}/{},family_share={:.2},tile_founder_ratio={:.2},tile_jitter={};peaks={{{}}}",
                    m_f,
                    mult,
                    tile,
                    family_size_best,
                    peaks.len(),
                    if total_score > 0.0 { family_score_best / total_score } else { 0.0 },
                    tile_founder_ratio,
                    tile_jitter,
                    fam_str,
                ),
            };
        }
    }

    // §2b alt path — jitter triggers Unresolved if no family fit
    if jitter >= cfg.jitter_thr {
        return HorCallResult {
            verdict: HorVerdict::Unresolved,
            founder_bp: None,
            multiplicity: None,
            tile_bp: None,
            family_score: 0.0,
            family_size: 0,
            jitter,
            reason: format!(
                "unresolved:jitter={}/{},top3={};no_family_fit",
                jitter,
                peaks.len().min(cfg.n_jitter),
                describe_top3(peaks),
            ),
        };
    }

    // §2d — tandem vs unresolved by dominance.
    let s1 = peaks[0].score;
    let s2 = peaks[1].score;
    if s2 > 0.0 && s1 / s2 >= cfg.dominance {
        return HorCallResult {
            verdict: HorVerdict::Tandem,
            founder_bp: Some(top1.period),
            multiplicity: Some(1),
            tile_bp: Some(top1.period),
            family_score: s1,
            family_size: 1,
            jitter,
            reason: format!(
                "tandem:monomer={},score={:.6},dominance={:.2};top3={}",
                top1.period,
                s1,
                s1 / s2,
                describe_top3(peaks),
            ),
        };
    }
    HorCallResult {
        verdict: HorVerdict::Unresolved,
        founder_bp: None,
        multiplicity: None,
        tile_bp: None,
        family_score: family_score_best,
        family_size: family_size_best,
        jitter,
        reason: format!(
            "unresolved:top3={};dominance={:.2}(<{}),family={}/{}",
            describe_top3(peaks),
            if s2 > 0.0 { s1 / s2 } else { 0.0 },
            cfg.dominance,
            family_size_best,
            cfg.min_family_size,
        ),
    }
}

/// Find the best HOR founder among kite peaks. Candidates iterated in
/// score-descending order; the first to pass ALL the pre/post checks
/// wins. Returns
/// `(best_m_f, family_score, family_size, best_tile_k, best_tile_period)`
/// or all-zero/None if no valid candidate.
///
/// Checks:
/// - **Pre-filter**: `founder_score ≥ min_founder_top1_share × top1_score`
///   — rejects spurious small-period founders (TRC_1__15071287 m_f=60).
/// - **Family size** ≥ `min_family_size`.
/// - **Top-K-in-family**: the top-`require_top_k_in_family` peaks (by
///   score) must all belong to the family — rejects cases where a small
///   founder happens to fit a family of low-score peaks while the top-
///   scoring peaks point to a different periodicity (TRC_1__15071287
///   top-3 = 178, 356, 238; 238 isn't a 178-multiple).
fn find_best_family(
    peaks: &[crate::kite::KitePeak],
    cfg: &HorCallConfig,
) -> (Option<usize>, f64, usize, usize, usize, f64, f64) {
    if peaks.is_empty() {
        return (None, 0.0, 0, 0, 0, 0.0, 0.0);
    }
    let top1_score = peaks[0].score;
    let mut ordered: Vec<&crate::kite::KitePeak> = peaks.iter().collect();
    ordered.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let n_top = cfg.require_top_k_in_family.min(peaks.len());
    let top_k: Vec<usize> = ordered.iter().take(n_top).map(|p| p.period).collect();

    for p in ordered.iter() {
        let m_f = p.period;
        if m_f < cfg.lo_period {
            continue;
        }
        // Pre-filter on founder score.
        if top1_score > 0.0 && p.score < cfg.min_founder_top1_share * top1_score {
            continue;
        }
        let (fam_score, fam_size, best_tile_k, best_tile_period, best_tile_score) =
            family_metrics(peaks, m_f, cfg);
        if fam_size < cfg.min_family_size {
            continue;
        }
        // Top-K-in-family check.
        let all_top_in = top_k
            .iter()
            .all(|tp| best_multiplicity(*tp, m_f, cfg).is_some());
        if !all_top_in {
            continue;
        }
        return (
            Some(m_f),
            fam_score,
            fam_size,
            best_tile_k,
            best_tile_period,
            best_tile_score,
            p.score,
        );
    }
    (None, 0.0, 0, 0, 0, 0.0, 0.0)
}

/// For a candidate founder `m_f`, return
/// `(family_score, family_size, best_tile_k, best_tile_period, best_tile_score)`.
///
/// `best_tile_*` identifies the highest-scoring family member with k ≥ 2
/// — the natural HOR-unit. `best_tile_score` is `-1.0` when no k≥2 peak
/// is in the family.
fn family_metrics(
    peaks: &[crate::kite::KitePeak],
    m_f: usize,
    cfg: &HorCallConfig,
) -> (f64, usize, usize, usize, f64) {
    let mut score = 0.0;
    let mut size = 0;
    let mut best_tile_k: usize = 1;
    let mut best_tile_period: usize = m_f;
    let mut best_tile_score: f64 = -1.0;
    for p in peaks {
        if let Some(k) = best_multiplicity(p.period, m_f, cfg) {
            score += p.score;
            size += 1;
            if k >= 2 && p.score > best_tile_score {
                best_tile_score = p.score;
                best_tile_k = k;
                best_tile_period = p.period;
            }
        }
    }
    (score, size, best_tile_k, best_tile_period, best_tile_score)
}

/// Round `p / m_f` to the nearest integer in `[1, qmax]`, return Some(k)
/// if `|p − k·m_f|` is within `tol_bp` or `tol_rel · (k·m_f)`. Else None.
fn best_multiplicity(p: usize, m_f: usize, cfg: &HorCallConfig) -> Option<usize> {
    if m_f == 0 {
        return None;
    }
    let k = ((p as f64) / (m_f as f64)).round() as usize;
    if k == 0 || k > cfg.qmax {
        return None;
    }
    let expected = k * m_f;
    let diff = p.abs_diff(expected);
    let tol = cfg.tol_bp.max((cfg.tol_rel * expected as f64) as usize);
    if diff <= tol {
        Some(k)
    } else {
        None
    }
}

fn compute_jitter(peaks: &[crate::kite::KitePeak], cfg: &HorCallConfig) -> usize {
    if peaks.is_empty() {
        return 0;
    }
    let n = peaks.len().min(cfg.n_jitter);
    let top1 = peaks[0].period;
    peaks_within_band(&peaks[..n], top1, cfg.jitter_tol)
}

/// Count peaks whose period is within `±tol × center` of `center`.
fn peaks_within_band(peaks: &[crate::kite::KitePeak], center: usize, tol: f64) -> usize {
    let band = tol * center as f64;
    peaks
        .iter()
        .filter(|p| {
            let diff = (p.period as f64 - center as f64).abs();
            diff <= band
        })
        .count()
}

fn describe_family(peaks: &[crate::kite::KitePeak], m_f: usize, cfg: &HorCallConfig) -> String {
    let mut parts = Vec::new();
    for p in peaks {
        if let Some(k) = best_multiplicity(p.period, m_f, cfg) {
            parts.push(format!("{}={}*{}(s={:.4})", p.period, k, m_f, p.score));
        }
    }
    parts.join(",")
}

fn describe_top3(peaks: &[crate::kite::KitePeak]) -> String {
    let mut parts = Vec::new();
    for p in peaks.iter().take(3) {
        parts.push(format!("{}(s={:.4})", p.period, p.score));
    }
    parts.join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kite::KitePeak;

    fn peak(period: usize, score: f64) -> KitePeak {
        KitePeak {
            period,
            peak_height: 0.0,
            score,
            score2: 0.0,
            score2_norm: 0.0,
            background: 0.0,
        }
    }

    fn result_with(peaks: Vec<KitePeak>) -> KiteResult {
        KiteResult {
            array_id: "test".into(),
            length_bp: 100_000,
            peaks,
            profile: None,
            background: None,
        }
    }

    #[test]
    fn empty_peaks_is_no_signal() {
        let r = classify(&result_with(vec![]), &HorCallConfig::default());
        assert_eq!(r.verdict, HorVerdict::NoSignal);
    }

    #[test]
    fn single_peak_is_tandem() {
        let r = classify(
            &result_with(vec![peak(178, 0.5)]),
            &HorCallConfig::default(),
        );
        assert_eq!(r.verdict, HorVerdict::Tandem);
        assert_eq!(r.founder_bp, Some(178));
        assert_eq!(r.multiplicity, Some(1));
    }

    #[test]
    fn clean_hor_family_fires() {
        // TRC_1__14926944 shape: peaks at 178, 356, 532, 712, 888 in
        // a clean 178-family.
        let peaks = vec![
            peak(888, 0.30),
            peak(178, 0.20),
            peak(356, 0.18),
            peak(712, 0.10),
            peak(532, 0.08),
        ];
        let r = classify(&result_with(peaks), &HorCallConfig::default());
        assert_eq!(r.verdict, HorVerdict::Hor);
        assert_eq!(r.founder_bp, Some(178));
        assert_eq!(r.multiplicity, Some(5));
        assert_eq!(r.tile_bp, Some(888)); // actual peak period
    }

    #[test]
    fn dominant_single_peak_is_tandem() {
        // top1 much larger than top2 → tandem
        let peaks = vec![peak(178, 0.50), peak(800, 0.02)];
        let r = classify(&result_with(peaks), &HorCallConfig::default());
        assert_eq!(r.verdict, HorVerdict::Tandem);
        assert_eq!(r.founder_bp, Some(178));
    }

    #[test]
    fn variable_length_tandem_is_unresolved() {
        // 8 peaks clustered within ±15 % of top1 (TRC_4 shape) with no
        // clean family fit.
        let peaks = vec![
            peak(10461, 0.20),
            peak(10050, 0.18),
            peak(9384, 0.15),
            peak(9037, 0.13),
            peak(9792, 0.11),
            peak(10113, 0.09),
            peak(10262, 0.08),
            peak(9387, 0.06),
        ];
        let r = classify(&result_with(peaks), &HorCallConfig::default());
        assert_eq!(r.verdict, HorVerdict::Unresolved);
        assert!(r.jitter >= 4);
    }

    #[test]
    fn submono_family_picks_smallest_founder() {
        // horsubmono shape: peaks at 53, 212(=4·53), 636(=12·53).
        // Founder=53, tile=636, k=12.
        let peaks = vec![peak(636, 0.30), peak(53, 0.20), peak(212, 0.15)];
        let r = classify(&result_with(peaks), &HorCallConfig::default());
        assert_eq!(r.verdict, HorVerdict::Hor);
        assert_eq!(r.founder_bp, Some(53));
        assert_eq!(r.multiplicity, Some(12));
        assert_eq!(r.tile_bp, Some(636));
    }
}

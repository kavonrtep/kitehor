//! Algorithm: candidate selection → per-record window plan → per-window
//! in-process kite → density + spatial + phase contrasts → per-record
//! decision. Mirrors `tools/rule_proto/tandem_validate.py` (spec v5).

use super::io::decision_label;
use super::Config;
use crate::kite::{analyze as kite_analyze, KiteConfig};
use crate::sequence::ArrayRecord;
use ahash::AHashMap;
use rayon::prelude::*;
use std::cmp::Ordering;

/// One row of the kite peaks TSV as the candidate picker needs it.
/// Owned here so the module is self-contained once `subrepeat` is
/// retired in commit 2.
#[derive(Debug, Clone, Copy)]
pub struct PeakRow {
    pub rank: u32,
    pub period: usize,
    pub score2_norm: f64,
}

/// One row of `verdicts.tsv` — kept for **every** verdict, not just HOR.
#[derive(Debug, Clone)]
pub struct VerdictRow {
    pub case_id: String,
    pub verdict: String,
    pub founder: Option<f64>,
    pub tile: Option<f64>,
    pub multiplicity: Option<i64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateKind {
    Founder,
    Other,
}

#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    pub period: f64,
    /// `-1` for founder, else the kite rank (1-based).
    pub rank: i32,
    pub score2_norm: f64,
    pub kind: CandidateKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Decision {
    LocalizedSubrepeat,
    ConfirmsHost,
    Ambiguous,
    NoSignal,
    NoCandidates,
    NoWindows,
    SkipK2,
    NoHost,
    NoVerdict,
}

#[derive(Debug, Clone)]
pub struct Row {
    pub record_id: String,
    pub verdict: String,
    pub host_period: Option<f64>,
    pub multiplicity: Option<i64>,
    pub window_bp: Option<usize>,
    pub n_candidates: usize,
    /// Semicolon-joined per-candidate diagnostic:
    /// `{kind}/{period}:d={density}:sc={sc}:pc={pc}:{decision}`.
    pub candidates: String,
    pub best_candidate_period: Option<f64>,
    pub best_candidate_kind: Option<CandidateKind>,
    pub density: Option<f64>,
    pub spatial_contrast: Option<f64>,
    pub phase_contrast: Option<f64>,
    pub n_windows_total: usize,
    pub n_windows_present: usize,
    pub decision_hint: Decision,
    /// `{kind}:{decision}` for real outcomes, bare skip-string otherwise.
    pub reason: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CandDecision {
    Localized,
    Ambiguous,
    Uniform,
    NoSignal,
}

fn cand_decision_label(d: CandDecision) -> &'static str {
    match d {
        CandDecision::Localized => "localized",
        CandDecision::Ambiguous => "ambiguous",
        CandDecision::Uniform => "uniform",
        CandDecision::NoSignal => "no_signal",
    }
}

fn kind_label(k: CandidateKind) -> &'static str {
    match k {
        CandidateKind::Founder => "founder",
        CandidateKind::Other => "other",
    }
}

/// Pick the host period for a verdict row. Returns `None` when no host
/// can be derived (kite peaks empty AND verdict provides none).
///
/// Mirrors the Python prototype's fall-through: try the
/// verdict-specific source first (`tile` for HOR, `founder` for
/// `simple_tr`), then fall back to the strongest kite peak. The
/// fall-back fires for `unresolved` verdicts AND for HOR/simple_tr rows
/// where the verdict-specific column is missing or zero.
pub fn pick_host_period(verdict: &VerdictRow, kite_peaks: &[PeakRow]) -> Option<f64> {
    let from_verdict = match verdict.verdict.as_str() {
        "hor" => verdict.tile.filter(|t| *t > 0.0),
        "simple_tr" => verdict.founder.filter(|f| *f > 0.0),
        _ => None,
    };
    if from_verdict.is_some() {
        return from_verdict;
    }
    kite_peaks
        .iter()
        .max_by(|a, b| {
            a.score2_norm
                .partial_cmp(&b.score2_norm)
                .unwrap_or(Ordering::Equal)
        })
        .map(|p| p.period as f64)
}

/// Reason a record was skipped before candidate scoring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Gate {
    Proceed,
    SkipK2,
    NoHost,
}

/// Return the candidate list (+ a gate). `Gate::Proceed` means scan;
/// `SkipK2` / `NoHost` mean stop and emit the corresponding skip row.
pub fn pick_candidates(
    verdict: &str,
    founder: Option<f64>,
    multiplicity: i64,
    host_period: f64,
    kite_peaks: &[PeakRow],
    cfg: &Config,
) -> (Vec<Candidate>, Gate) {
    if host_period <= 0.0 || kite_peaks.is_empty() {
        return (Vec::new(), Gate::NoHost);
    }
    if verdict == "hor" && multiplicity == 2 {
        return (Vec::new(), Gate::SkipK2);
    }

    let cand_max = host_period * cfg.host_inside_ratio;
    let record_max_score = kite_peaks
        .iter()
        .map(|p| p.score2_norm)
        .fold(0.0_f64, f64::max);
    let eff_floor = cfg
        .cand_score_floor
        .max(cfg.cand_rel_score_floor * record_max_score);

    let mut out: Vec<Candidate> = Vec::new();

    if verdict == "hor" && multiplicity >= 3 {
        if let Some(f) = founder.filter(|f| *f > 0.0) {
            // Founder candidate. At k=3, `founder == host/3` exactly, so
            // a 1% slack admits the borderline.
            if f <= cand_max * 1.01 {
                let score: f64 = kite_peaks
                    .iter()
                    .filter(|p| ((p.period as f64) - f).abs() / f <= cfg.founder_tol)
                    .map(|p| p.score2_norm)
                    .sum();
                out.push(Candidate {
                    period: f,
                    rank: -1,
                    score2_norm: score,
                    kind: CandidateKind::Founder,
                });
            }
            // Non-founder candidates exclude HOR-ladder rungs (m·founder
            // for m = 1..k-1; the m=1 case excludes the founder itself
            // from this loop since it's already added above).
            let explained: Vec<f64> = (1..multiplicity).map(|m| (m as f64) * f).collect();
            out.extend(pick_other_candidates(
                kite_peaks, cand_max, eff_floor, cfg, &explained,
            ));
            return (out, Gate::Proceed);
        }
    }

    // simple_tr / unresolved / hor-without-founder
    out.extend(pick_other_candidates(
        kite_peaks,
        cand_max,
        eff_floor,
        cfg,
        &[],
    ));
    (out, Gate::Proceed)
}

fn pick_other_candidates(
    kite_peaks: &[PeakRow],
    cand_max: f64,
    eff_floor: f64,
    cfg: &Config,
    explained: &[f64],
) -> Vec<Candidate> {
    // Sort peaks by score2_norm descending (stable for ties to mirror
    // pandas's sort_values stability).
    let mut sorted: Vec<&PeakRow> = kite_peaks
        .iter()
        .filter(|p| {
            (p.period as f64) >= cfg.cand_min_period as f64
                && (p.period as f64) < cand_max
                && p.score2_norm >= eff_floor
        })
        .collect();
    sorted.sort_by(|a, b| {
        b.score2_norm
            .partial_cmp(&a.score2_norm)
            .unwrap_or(Ordering::Equal)
    });

    let mut kept: Vec<Candidate> = Vec::new();
    for p in sorted {
        let per = p.period as f64;
        if explained
            .iter()
            .any(|x| *x > 0.0 && (per - *x).abs() / *x <= cfg.founder_tol)
        {
            continue;
        }
        kept.push(Candidate {
            period: per,
            rank: p.rank as i32,
            score2_norm: p.score2_norm,
            kind: CandidateKind::Other,
        });
        if kept.len() >= cfg.cand_top_n {
            break;
        }
    }
    kept
}

/// Plan sliding windows for a record. Window size is capped at the host
/// period to preserve within-host phase resolution.
pub fn plan_windows(
    seq_len: usize,
    host_period: f64,
    max_cand: f64,
    cfg: &Config,
) -> Vec<(usize, usize)> {
    let w_f = (host_period * cfg.window_host_frac)
        .max(max_cand * cfg.window_cand_mult)
        .max(cfg.min_window_bp as f64);
    let w_f = w_f.min(host_period); // hard cap at host
    let w = w_f.round() as usize;
    if w == 0 {
        return Vec::new();
    }
    if w >= seq_len {
        return vec![(0, seq_len)];
    }
    let step = (w / 4).max(1);
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut s = 0usize;
    while s + w <= seq_len {
        out.push((s, s + w));
        s += step;
    }
    if let Some(last) = out.last().copied() {
        if last.1 < seq_len {
            let start = seq_len.saturating_sub(w);
            out.push((start, seq_len));
        }
    }
    out
}

/// Place `value` into one of `n_bins` equal-width bins of `[lo, lo + span)`.
/// Clamps out-of-range values.
pub fn bin_index(value: f64, lo: f64, span: f64, n_bins: usize) -> usize {
    if span <= 0.0 || n_bins == 0 {
        return 0;
    }
    let raw = ((value - lo) / span * n_bins as f64) as i64;
    if raw < 0 {
        0
    } else if raw as usize >= n_bins {
        n_bins - 1
    } else {
        raw as usize
    }
}

fn kite_on_window(win_id: String, seq_slice: &[u8]) -> Vec<PeakRow> {
    if seq_slice.is_empty() {
        return Vec::new();
    }
    let rec = ArrayRecord::from_raw(win_id, seq_slice);
    let cfg = KiteConfig::default();
    let result = kite_analyze(&rec, &cfg);
    result
        .peaks
        .into_iter()
        .enumerate()
        .map(|(i, p)| PeakRow {
            rank: (i + 1) as u32,
            period: p.period,
            score2_norm: p.score2_norm,
        })
        .collect()
}

/// Per-window per-candidate presence test. Two criteria gated by
/// candidate kind:
/// * **Founder** (loose): sum of peaks within `±tol` of candidate
///   period is `≥ presence_rel_floor × top_score`. The founder is
///   rank-2/3 of kite in a clean HOR (the tile is rank-1), so a strict
///   "top must be founder" check would never fire.
/// * **Other** (strict): the window's top period is the candidate AND
///   `top_score ≥ window_score_floor`. The score floor distinguishes
///   real heterogeneity (some windows have weak tops) from uniform
///   tandem (all windows have strong tops).
fn is_present(cand: &Candidate, w_peaks: &[PeakRow], top_score: f64, cfg: &Config) -> bool {
    if w_peaks.is_empty() || top_score <= 0.0 {
        return false;
    }
    match cand.kind {
        CandidateKind::Founder => {
            let f_score: f64 = w_peaks
                .iter()
                .filter(|p| {
                    ((p.period as f64) - cand.period).abs() / cand.period <= cfg.period_match_tol
                })
                .map(|p| p.score2_norm)
                .sum();
            (f_score / top_score) >= cfg.presence_rel_floor
        }
        CandidateKind::Other => {
            let Some(top) = w_peaks.iter().max_by(|a, b| {
                a.score2_norm
                    .partial_cmp(&b.score2_norm)
                    .unwrap_or(Ordering::Equal)
            }) else {
                return false;
            };
            let top_period = top.period as f64;
            ((top_period - cand.period).abs() / cand.period <= cfg.period_match_tol)
                && (top.score2_norm >= cfg.window_score_floor)
        }
    }
}

/// Per-candidate metrics + decision.
#[derive(Debug, Clone)]
struct CandidateScore {
    candidate: Candidate,
    density: f64,
    spatial_contrast: f64,
    phase_contrast: Option<f64>,
    n_present: usize,
    decision: CandDecision,
}

/// Per-record window context handed to [`classify_candidate`]. Bundles
/// the data that doesn't vary by candidate so the classifier has a
/// small signature.
struct WindowCtx<'a> {
    wins: &'a [(usize, usize)],
    win_peaks: &'a [Option<Vec<PeakRow>>],
    win_top_scores: &'a [f64],
    seq_lo: usize,
    seq_span: usize,
    host: f64,
    window_bp: usize,
}

fn classify_candidate(cand: &Candidate, ctx: &WindowCtx<'_>, cfg: &Config) -> CandidateScore {
    let wins = ctx.wins;
    let n_total = wins.len();
    let mut present: Vec<u8> = Vec::with_capacity(n_total);
    for ((_, _), (tab, &top_score)) in wins
        .iter()
        .zip(ctx.win_peaks.iter().zip(ctx.win_top_scores.iter()))
    {
        let p = match tab {
            Some(t) => is_present(cand, t, top_score, cfg) as u8,
            None => 0,
        };
        present.push(p);
    }
    let n_present = present.iter().map(|p| *p as usize).sum::<usize>();
    let density = if n_total > 0 {
        n_present as f64 / n_total as f64
    } else {
        0.0
    };

    // Spatial contrast — 10 array-position bins, max-min over present-fractions.
    let mut spatial_total = vec![0_u32; cfg.n_bins];
    let mut spatial_pres = vec![0_u32; cfg.n_bins];
    for ((s, e), pres) in wins.iter().zip(present.iter()) {
        let mid = (*s as f64 + *e as f64) / 2.0;
        let b = bin_index(mid, ctx.seq_lo as f64, ctx.seq_span as f64, cfg.n_bins);
        spatial_total[b] += 1;
        spatial_pres[b] += *pres as u32;
    }
    let spatial_fracs: Vec<f64> = (0..cfg.n_bins)
        .filter(|i| spatial_total[*i] > 0)
        .map(|i| spatial_pres[i] as f64 / spatial_total[i] as f64)
        .collect();
    let spatial_contrast = if spatial_fracs.len() >= 2 {
        spatial_fracs
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            - spatial_fracs.iter().cloned().fold(f64::INFINITY, f64::min)
    } else {
        0.0
    };

    // Phase contrast — only meaningful when window < ~host (so midpoints
    // can hit different phase positions within one host cycle).
    let phase_contrast = if ctx.host > 0.0 && (ctx.window_bp as f64) < ctx.host * 0.95 {
        let cyc = ctx.host.round().max(1.0) as usize;
        let bin_width = ((ctx.host / cfg.n_bins as f64).round().max(1.0)) as usize;
        let mut phase_total = vec![0_u32; cfg.n_bins];
        let mut phase_pres = vec![0_u32; cfg.n_bins];
        for ((s, e), pres) in wins.iter().zip(present.iter()) {
            let mid = (s + e) / 2;
            let mut b = (mid % cyc) / bin_width;
            if b >= cfg.n_bins {
                b = cfg.n_bins - 1;
            }
            phase_total[b] += 1;
            phase_pres[b] += *pres as u32;
        }
        let fracs: Vec<f64> = (0..cfg.n_bins)
            .filter(|i| phase_total[*i] > 0)
            .map(|i| phase_pres[i] as f64 / phase_total[i] as f64)
            .collect();
        if fracs.len() >= 2 {
            Some(
                fracs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                    - fracs.iter().cloned().fold(f64::INFINITY, f64::min),
            )
        } else {
            None
        }
    } else {
        None
    };

    let is_loc = density <= cfg.density_dup_max
        || spatial_contrast >= cfg.contrast_dup_min
        || phase_contrast.is_some_and(|pc| pc >= cfg.contrast_dup_min);
    let is_uni = density >= cfg.density_hor_min
        && spatial_contrast <= cfg.contrast_hor_max
        && phase_contrast.is_none_or(|pc| pc <= cfg.contrast_hor_max);

    let decision = if n_present < cfg.min_present_windows {
        CandDecision::NoSignal
    } else if is_loc {
        CandDecision::Localized
    } else if is_uni {
        CandDecision::Uniform
    } else {
        CandDecision::Ambiguous
    };

    CandidateScore {
        candidate: *cand,
        density,
        spatial_contrast,
        phase_contrast,
        n_present,
        decision,
    }
}

fn decision_rank(d: CandDecision) -> u8 {
    match d {
        CandDecision::Localized => 0,
        CandDecision::Ambiguous => 1,
        CandDecision::Uniform => 2,
        CandDecision::NoSignal => 3,
    }
}

fn round_n(x: f64, n: i32) -> f64 {
    let p = 10_f64.powi(n);
    (x * p).round() / p
}

fn format_candidate_diag(cs: &CandidateScore) -> String {
    let pc = cs
        .phase_contrast
        .map(|x| format!("{:.3}", x))
        .unwrap_or_else(|| "nan".to_string());
    format!(
        "{}/{:.0}:d={:.3}:sc={:.3}:pc={}:{}",
        kind_label(cs.candidate.kind),
        cs.candidate.period,
        cs.density,
        cs.spatial_contrast,
        pc,
        cand_decision_label(cs.decision),
    )
}

fn empty_skip_row(
    rec_id: &str,
    verdict: &str,
    host: Option<f64>,
    mult: Option<i64>,
    d: Decision,
) -> Row {
    Row {
        record_id: rec_id.to_string(),
        verdict: verdict.to_string(),
        host_period: host.map(|h| round_n(h, 2)),
        multiplicity: mult,
        window_bp: None,
        n_candidates: 0,
        candidates: String::new(),
        best_candidate_period: None,
        best_candidate_kind: None,
        density: None,
        spatial_contrast: None,
        phase_contrast: None,
        n_windows_total: 0,
        n_windows_present: 0,
        decision_hint: d,
        reason: decision_label(d).to_string(),
    }
}

/// Scan one record. Pure function — no I/O.
pub fn scan_record(
    rec_id: &str,
    seq: &[u8],
    verdict: Option<&VerdictRow>,
    kite_peaks: &[PeakRow],
    cfg: &Config,
) -> Row {
    let Some(v) = verdict else {
        return empty_skip_row(rec_id, "", None, None, Decision::NoVerdict);
    };
    let host = match pick_host_period(v, kite_peaks) {
        Some(h) => h,
        None => return empty_skip_row(rec_id, &v.verdict, None, v.multiplicity, Decision::NoHost),
    };
    let mult = v.multiplicity.unwrap_or(1).max(1);
    let (cands, gate) = pick_candidates(&v.verdict, v.founder, mult, host, kite_peaks, cfg);
    match gate {
        Gate::SkipK2 => {
            return empty_skip_row(rec_id, &v.verdict, Some(host), Some(mult), Decision::SkipK2);
        }
        Gate::NoHost => {
            return empty_skip_row(rec_id, &v.verdict, Some(host), Some(mult), Decision::NoHost);
        }
        Gate::Proceed => {}
    }
    if cands.is_empty() {
        return empty_skip_row(
            rec_id,
            &v.verdict,
            Some(host),
            Some(mult),
            Decision::NoCandidates,
        );
    }
    let max_cand = cands.iter().map(|c| c.period).fold(0.0_f64, f64::max);
    let wins = plan_windows(seq.len(), host, max_cand, cfg);
    if wins.is_empty() {
        return empty_skip_row(
            rec_id,
            &v.verdict,
            Some(host),
            Some(mult),
            Decision::NoWindows,
        );
    }
    let window_bp = wins[0].1 - wins[0].0;

    // Per-window kite (parallel across windows of this record).
    let kite_per_window: Vec<Option<Vec<PeakRow>>> = wins
        .par_iter()
        .map(|(s, e)| {
            let wid = format!("{rec_id}__W_{s}_{e}");
            let peaks = kite_on_window(wid, &seq[*s..*e]);
            if peaks.is_empty() {
                None
            } else {
                Some(peaks)
            }
        })
        .collect();
    let win_top_scores: Vec<f64> = kite_per_window
        .iter()
        .map(|tab| match tab {
            Some(t) => t.iter().map(|p| p.score2_norm).fold(0.0_f64, f64::max),
            None => 0.0,
        })
        .collect();

    let seq_lo = wins[0].0;
    let seq_hi = wins.last().unwrap().1;
    let seq_span = seq_hi.saturating_sub(seq_lo).max(1);
    let ctx = WindowCtx {
        wins: &wins,
        win_peaks: &kite_per_window,
        win_top_scores: &win_top_scores,
        seq_lo,
        seq_span,
        host,
        window_bp,
    };

    let mut scores: Vec<CandidateScore> = cands
        .iter()
        .map(|c| classify_candidate(c, &ctx, cfg))
        .collect();

    // Best candidate: localized > ambiguous > uniform > no_signal,
    // secondary key = -(spatial_contrast + phase_contrast.unwrap_or(0)).
    scores.sort_by(|a, b| {
        let ka = (
            decision_rank(a.decision),
            -(a.spatial_contrast + a.phase_contrast.unwrap_or(0.0)),
        );
        let kb = (
            decision_rank(b.decision),
            -(b.spatial_contrast + b.phase_contrast.unwrap_or(0.0)),
        );
        ka.0.cmp(&kb.0)
            .then_with(|| ka.1.partial_cmp(&kb.1).unwrap_or(Ordering::Equal))
    });
    let best = scores[0].clone();

    let (record_decision, reason_decision) = match best.decision {
        CandDecision::Localized => (Decision::LocalizedSubrepeat, "localized"),
        CandDecision::Uniform => (Decision::ConfirmsHost, "uniform"),
        CandDecision::Ambiguous => (Decision::Ambiguous, "ambiguous"),
        CandDecision::NoSignal => (Decision::NoSignal, "no_signal"),
    };
    let reason = format!("{}:{}", kind_label(best.candidate.kind), reason_decision);

    let candidates_str = scores
        .iter()
        .map(format_candidate_diag)
        .collect::<Vec<_>>()
        .join(";");

    Row {
        record_id: rec_id.to_string(),
        verdict: v.verdict.clone(),
        host_period: Some(round_n(host, 2)),
        multiplicity: Some(mult),
        window_bp: Some(window_bp),
        n_candidates: scores.len(),
        candidates: candidates_str,
        best_candidate_period: Some(round_n(best.candidate.period, 1)),
        best_candidate_kind: Some(best.candidate.kind),
        density: Some(round_n(best.density, 4)),
        spatial_contrast: Some(round_n(best.spatial_contrast, 4)),
        phase_contrast: best.phase_contrast.map(|x| round_n(x, 4)),
        n_windows_total: wins.len(),
        n_windows_present: best.n_present,
        decision_hint: record_decision,
        reason,
    }
}

/// Scan all records in parallel. Per-window kite calls within a record
/// are also parallelized via [`scan_record`].
pub fn scan_records(
    records: &[(String, Vec<u8>)],
    verdicts: &[VerdictRow],
    kite_peaks: &AHashMap<String, Vec<PeakRow>>,
    cfg: &Config,
) -> Vec<Row> {
    let by_id: AHashMap<&str, &VerdictRow> =
        verdicts.iter().map(|v| (v.case_id.as_str(), v)).collect();
    let empty: Vec<PeakRow> = Vec::new();
    records
        .par_iter()
        .map(|(rec_id, seq)| {
            let v = by_id.get(rec_id.as_str()).copied();
            let peaks = kite_peaks.get(rec_id).unwrap_or(&empty);
            scan_record(rec_id, seq, v, peaks, cfg)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak(rank: u32, period: usize, s2n: f64) -> PeakRow {
        PeakRow {
            rank,
            period,
            score2_norm: s2n,
        }
    }

    fn cfg() -> Config {
        Config::default()
    }

    #[test]
    fn picks_no_candidates_for_clean_tandem() {
        // simple_tr at monomer=171; kite shows monomer + 2× + 3× harmonics.
        // All harmonics sit at or above host, so cand_max=host/3=57 leaves
        // nothing. Result: empty candidate list, no_candidates path.
        let peaks = vec![peak(1, 171, 0.94), peak(2, 342, 0.40), peak(3, 513, 0.20)];
        let (c, gate) = pick_candidates("simple_tr", Some(171.0), 1, 171.0, &peaks, &cfg());
        assert_eq!(gate, Gate::Proceed);
        assert!(c.is_empty(), "expected no candidates, got {:?}", c);
    }

    #[test]
    fn picks_founder_for_clean_hor_k3() {
        // HOR k=3: tile=510, founder=170. Kite shows tile (top) and founder.
        // Founder must be added as kind=Founder; no Other candidates (founder
        // is excluded by explained=[170] rung).
        let peaks = vec![peak(1, 510, 0.80), peak(2, 170, 0.30)];
        let (c, gate) = pick_candidates("hor", Some(170.0), 3, 510.0, &peaks, &cfg());
        assert_eq!(gate, Gate::Proceed);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].kind, CandidateKind::Founder);
        assert!((c[0].period - 170.0).abs() < 1e-9);
    }

    #[test]
    fn excludes_hor_ladder_harmonics() {
        // k=3 HOR with founder=503, tile=1509. Kite peaks at 503, 1006, 1512.
        // Founder added; 1006 = 2×503 must be excluded by the explained-rung
        // filter. 1512 > host/3=503 → filtered by cand_max anyway.
        let peaks = vec![peak(1, 1512, 0.70), peak(2, 503, 0.40), peak(3, 1006, 0.30)];
        let (c, gate) = pick_candidates("hor", Some(503.0), 3, 1509.0, &peaks, &cfg());
        assert_eq!(gate, Gate::Proceed);
        assert_eq!(c.len(), 1, "only founder should remain; got {:?}", c);
        assert_eq!(c[0].kind, CandidateKind::Founder);
    }

    #[test]
    fn respects_cand_max() {
        // host=900 → cand_max=300. A peak at 250 (< 300) is kept, peak at
        // 450 (> 300) is dropped.
        let peaks = vec![peak(1, 900, 0.80), peak(2, 250, 0.30), peak(3, 450, 0.30)];
        let (c, gate) = pick_candidates("simple_tr", Some(900.0), 1, 900.0, &peaks, &cfg());
        assert_eq!(gate, Gate::Proceed);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].period as usize, 250);
    }

    #[test]
    fn applies_rel_score_floor() {
        // record_max=0.94, rel_floor=0.03 ⇒ eff_floor ≈ 0.0282.
        // Peak with score 0.025 is dropped, score 0.030 is kept.
        let peaks = vec![peak(1, 600, 0.94), peak(2, 150, 0.030), peak(3, 170, 0.025)];
        let (c, gate) = pick_candidates("simple_tr", Some(600.0), 1, 600.0, &peaks, &cfg());
        assert_eq!(gate, Gate::Proceed);
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].period as usize, 150);
    }

    #[test]
    fn skip_k2_for_hor_k2() {
        let peaks = vec![peak(1, 340, 0.80), peak(2, 170, 0.30)];
        let (c, gate) = pick_candidates("hor", Some(170.0), 2, 340.0, &peaks, &cfg());
        assert_eq!(gate, Gate::SkipK2);
        assert!(c.is_empty());
    }

    #[test]
    fn no_host_when_peaks_empty_and_unresolved() {
        let (c, gate) = pick_candidates("unresolved", None, 1, 0.0, &[], &cfg());
        assert_eq!(gate, Gate::NoHost);
        assert!(c.is_empty());
    }

    #[test]
    fn pick_host_period_hor_uses_tile() {
        let v = VerdictRow {
            case_id: "x".into(),
            verdict: "hor".into(),
            founder: Some(170.0),
            tile: Some(510.0),
            multiplicity: Some(3),
        };
        assert_eq!(pick_host_period(&v, &[]), Some(510.0));
    }

    #[test]
    fn pick_host_period_simple_tr_uses_founder() {
        let v = VerdictRow {
            case_id: "x".into(),
            verdict: "simple_tr".into(),
            founder: Some(171.0),
            tile: None,
            multiplicity: Some(1),
        };
        assert_eq!(pick_host_period(&v, &[]), Some(171.0));
    }

    #[test]
    fn pick_host_period_unresolved_falls_back_to_kite_top() {
        let v = VerdictRow {
            case_id: "x".into(),
            verdict: "unresolved".into(),
            founder: None,
            tile: None,
            multiplicity: None,
        };
        let peaks = vec![peak(1, 401, 0.10), peak(2, 800, 0.50)];
        assert_eq!(pick_host_period(&v, &peaks), Some(800.0));
    }

    #[test]
    fn pick_host_period_simple_tr_without_founder_falls_back_to_kite() {
        // simple_tr with `founder=NA` (e.g., kite-classify k=1 output)
        // still has a valid host — fall back to the kite top.
        let v = VerdictRow {
            case_id: "tandem_pure".into(),
            verdict: "simple_tr".into(),
            founder: None,
            tile: None,
            multiplicity: Some(1),
        };
        let peaks = vec![peak(1, 60, 1.0)];
        assert_eq!(pick_host_period(&v, &peaks), Some(60.0));
    }

    #[test]
    fn plan_windows_caps_at_host() {
        // host=600 → host_frac floor = 200, cand_mult = 3*180 = 540, min=200
        // ⇒ w = 540, capped at 600. step=540/4=135.
        let w = plan_windows(3000, 600.0, 180.0, &cfg());
        assert!(w.iter().all(|(s, e)| e - s == 540));
        assert_eq!(w[1].0 - w[0].0, 135);
    }

    #[test]
    fn plan_windows_min_window_floor() {
        // host=180 → host_frac=60, cand_mult=3*30=90, min=200 ⇒ w = 200,
        // capped at host=180 ⇒ 180.
        let w = plan_windows(3000, 180.0, 30.0, &cfg());
        assert!(w.iter().all(|(s, e)| e - s == 180));
    }

    #[test]
    fn plan_windows_flush_tail() {
        // Tail window padded out when last regular window doesn't reach end.
        let w = plan_windows(1050, 600.0, 50.0, &cfg());
        // w = max(200, 150, 200) = 200, capped at 600 ⇒ 200. step=50.
        // 1050 / 50 = 21, but final flush is appended only if last.end < seq.
        assert_eq!(w.last().map(|p| p.1), Some(1050));
    }

    #[test]
    fn bin_index_clamps() {
        assert_eq!(bin_index(-5.0, 0.0, 100.0, 10), 0);
        assert_eq!(bin_index(100.0, 0.0, 100.0, 10), 9);
        assert_eq!(bin_index(55.0, 0.0, 100.0, 10), 5);
    }

    #[test]
    fn is_present_founder_loose() {
        let c = Candidate {
            period: 170.0,
            rank: -1,
            score2_norm: 0.0,
            kind: CandidateKind::Founder,
        };
        // Window top = tile=510 at 0.80; founder=170 at 0.20. Ratio = 0.25 ≥ 0.2 → present.
        let w = vec![peak(1, 510, 0.80), peak(2, 170, 0.20)];
        assert!(is_present(&c, &w, 0.80, &cfg()));
        // Now founder=0.10 → ratio 0.125 < 0.2 → absent.
        let w2 = vec![peak(1, 510, 0.80), peak(2, 170, 0.10)];
        assert!(!is_present(&c, &w2, 0.80, &cfg()));
    }

    #[test]
    fn is_present_other_strict() {
        let c = Candidate {
            period: 170.0,
            rank: 1,
            score2_norm: 0.0,
            kind: CandidateKind::Other,
        };
        // Top IS candidate and top_score ≥ 0.3 → present.
        let w = vec![peak(1, 170, 0.40)];
        assert!(is_present(&c, &w, 0.40, &cfg()));
        // Top IS candidate but score < 0.3 → absent.
        let w2 = vec![peak(1, 170, 0.20)];
        assert!(!is_present(&c, &w2, 0.20, &cfg()));
        // Top is NOT candidate → absent.
        let w3 = vec![peak(1, 510, 0.50), peak(2, 170, 0.40)];
        assert!(!is_present(&c, &w3, 0.50, &cfg()));
    }
}

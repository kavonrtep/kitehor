//! Per-record subrepeat scan: pick candidates → window sequence →
//! per-window in-process kite → smooth → block → flag.

use super::Config;
use crate::kite::{analyze as kite_analyze, KiteConfig};
use crate::sequence::ArrayRecord;

/// One kite-peak row, as needed by the candidate picker.
#[derive(Debug, Clone, Copy)]
pub struct PeakRow {
    pub rank: u32,
    pub period: usize,
    pub score2_norm: f64,
}

#[derive(Debug, Clone)]
pub struct SummaryRow {
    pub record_id: String,
    pub length_bp: usize,
    /// "NA" or int.
    pub host_period_bp: String,
    pub subrepeat_period_bp: String,
    pub subrepeat_flag: String, // yes / no / none
    pub reason: String,
    pub n_windows_total: usize,
    pub n_windows_sub: usize,
    pub n_windows_non_sub: usize,
    pub n_subrepeat_blocks: usize,
    pub subrepeat_coverage_bp: usize,
    /// Rounded to 2 decimals.
    pub subrepeat_coverage_pct: f64,
    /// "NA" or ";"-joined `s-e` pairs.
    pub blocks: String,
}

#[derive(Debug, Clone)]
pub struct WindowRow {
    pub record_id: String,
    pub window_start: usize,
    pub window_end: usize,
    /// -1 when no peaks.
    pub top_period: i64,
    pub top_score2_norm: f64,
    pub class_raw: String, // sub / non_sub
    pub class_: String,    // post-smoothing
}

/// Pick `(sub_candidate, host_candidate)` periods. Returns `None` if no
/// qualifying pair exists.
///
/// **Audit note**: the prototype docstring says "longest" host period
/// but the code picks **strongest-scored** in the top-N host pool. We
/// match the code (not the docstring).
fn pick_candidates(peaks: &[PeakRow], cfg: &Config) -> Option<(f64, f64)> {
    if peaks.is_empty() {
        return None;
    }
    // `nlargest(top_n, score2_norm)` — sort by score desc and take first N.
    let mut by_score: Vec<&PeakRow> = peaks.iter().collect();
    by_score.sort_by(|a, b| {
        b.score2_norm
            .partial_cmp(&a.score2_norm)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sub_pool: Vec<&PeakRow> = by_score
        .iter()
        .take(cfg.top_n_sub)
        .filter(|p| p.score2_norm >= cfg.sub_floor)
        .copied()
        .collect();
    if sub_pool.is_empty() {
        return None;
    }
    let sub_p = sub_pool
        .iter()
        .map(|p| p.period as f64)
        .fold(f64::INFINITY, f64::min);
    let host_pool: Vec<&PeakRow> = by_score
        .iter()
        .take(cfg.top_n_host)
        .filter(|p| (p.period as f64) >= (cfg.host_sub_ratio_min as f64) * sub_p)
        .copied()
        .collect();
    if host_pool.is_empty() {
        return None;
    }
    let host_p = host_pool
        .iter()
        .max_by(|a, b| {
            a.score2_norm
                .partial_cmp(&b.score2_norm)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|p| p.period as f64)
        .unwrap();
    Some((sub_p, host_p))
}

fn windows_for_record(l: usize, sub_cand: f64, cfg: &Config) -> Vec<(usize, usize)> {
    let w = (cfg.window_mult_sub as f64 * sub_cand).round() as usize;
    let w = w.max(cfg.min_window_bp);
    if w >= l {
        return vec![(0, l)];
    }
    let step = (w / cfg.step_frac).max(1);
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut s = 0usize;
    while s + w <= l {
        out.push((s, s + w));
        s += step;
    }
    // Final flush window if tail uncovered.
    if let Some(last) = out.last() {
        if last.1 < l {
            let start = l.saturating_sub(w);
            out.push((start, l));
        }
    }
    out
}

fn classify_window_period(
    period: i64,
    score: f64,
    sub_cand: f64,
    tol: f64,
    score_floor: f64,
) -> &'static str {
    if sub_cand > 0.0 && score >= score_floor && period >= 0 {
        let p = period as f64;
        if (p - sub_cand).abs() / sub_cand <= tol {
            return "sub";
        }
    }
    "non_sub"
}

/// Morphological smoothing — fixed-point absorb of short runs into
/// longer neighbours. Tie-break for interior runs: previous neighbour
/// wins when `prev_len >= next_len`. Matches the prototype exactly.
fn smooth_runs(windows: &[(usize, usize, String)], min_run: usize) -> Vec<(usize, usize, String)> {
    if min_run <= 1 || windows.is_empty() {
        return windows.to_vec();
    }
    let mut out: Vec<(usize, usize, String)> = windows.to_vec();
    loop {
        // Compute runs of identical class labels.
        let mut runs: Vec<(usize, usize, String)> = Vec::new(); // (start_idx, end_idx_excl, cls)
        let mut i = 0usize;
        while i < out.len() {
            let mut j = i;
            while j < out.len() && out[j].2 == out[i].2 {
                j += 1;
            }
            runs.push((i, j, out[i].2.clone()));
            i = j;
        }
        if runs.len() <= 1 {
            break;
        }
        // Find shortest run with length < min_run.
        let mut short_idx: Option<usize> = None;
        for (k, (i, j, _)) in runs.iter().enumerate() {
            let len = j - i;
            if len < min_run {
                if let Some(s) = short_idx {
                    let cur_len = runs[s].1 - runs[s].0;
                    if len < cur_len {
                        short_idx = Some(k);
                    }
                } else {
                    short_idx = Some(k);
                }
            }
        }
        let Some(k) = short_idx else {
            break;
        };
        let (i, j, _) = (runs[k].0, runs[k].1, runs[k].2.clone());
        let new_cls = if k == 0 {
            runs[1].2.clone()
        } else if k == runs.len() - 1 {
            runs[runs.len() - 2].2.clone()
        } else {
            let prev_len = runs[k - 1].1 - runs[k - 1].0;
            let next_len = runs[k + 1].1 - runs[k + 1].0;
            if prev_len >= next_len {
                runs[k - 1].2.clone()
            } else {
                runs[k + 1].2.clone()
            }
        };
        for slot in out.iter_mut().take(j).skip(i) {
            slot.2 = new_cls.clone();
        }
    }
    out
}

fn blocks_from_windows(windows: &[(usize, usize, String)]) -> Vec<(usize, usize)> {
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut cur: Option<(usize, usize)> = None;
    for (s, e, cls) in windows {
        if cls == "sub" {
            cur = Some(match cur {
                None => (*s, *e),
                Some((cs, ce)) => (cs, ce.max(*e)),
            });
        } else if let Some(b) = cur.take() {
            blocks.push(b);
        }
    }
    if let Some(b) = cur {
        blocks.push(b);
    }
    blocks
}

/// Round to the same precision the `kite-periodicity --out-peaks` TSV
/// produces (`%.10f`). The prototype shells out to kite per window,
/// reads the TSV back through pandas (which loses bits beyond the
/// 10-decimal cap), and re-emits. To match byte-for-byte we mirror
/// that cap by formatting and re-parsing.
fn truncate_score(x: f64) -> f64 {
    format!("{:.10}", x).parse::<f64>().unwrap_or(x)
}

/// Run kite on `seq[start..end]` with the prototype-compatible synthetic
/// window-id (`{rec_id}__w{s}_{e}`), then return `(top_period,
/// top_score2_norm)` of the top peak by score2_norm. `(-1, 0.0)` when
/// kite returns no peaks. The score is rounded to `%.10f` precision
/// (see `truncate_score`).
fn kite_on_window(rec_id: &str, seq: &[u8], start: usize, end: usize) -> (i64, f64) {
    if end <= start {
        return (-1, 0.0);
    }
    let win_id = format!("{rec_id}__w{start}_{end}");
    let win_seq = &seq[start..end];
    let rec = ArrayRecord::from_raw(win_id, win_seq);
    let cfg = KiteConfig::default();
    let result = kite_analyze(&rec, &cfg, false);
    if result.peaks.is_empty() {
        return (-1, 0.0);
    }
    // **Rank-1 by `score`** (kite already sorts peaks score-desc). The
    // prototype filters `rank == 1` from the TSV; this is its
    // equivalent. Note that the rank-1 peak is NOT necessarily the one
    // with the highest `score2_norm` — kite computes
    // `score2 = score * log2(position)`, which reranks longer-period
    // peaks higher even when their bare `score` is lower.
    let top = &result.peaks[0];
    (top.period as i64, truncate_score(top.score2_norm))
}

pub fn scan_record(
    rec_id: &str,
    seq: &[u8],
    peaks: &[PeakRow],
    cfg: &Config,
) -> (SummaryRow, Vec<WindowRow>) {
    let length_bp = seq.len();
    let candidates = pick_candidates(peaks, cfg);
    let Some((sub_cand, host_cand)) = candidates else {
        let sum = SummaryRow {
            record_id: rec_id.to_string(),
            length_bp,
            host_period_bp: "NA".into(),
            subrepeat_period_bp: "NA".into(),
            subrepeat_flag: "none".into(),
            reason: "no_candidate_pair".into(),
            n_windows_total: 0,
            n_windows_sub: 0,
            n_windows_non_sub: 0,
            n_subrepeat_blocks: 0,
            subrepeat_coverage_bp: 0,
            subrepeat_coverage_pct: 0.0,
            blocks: "NA".into(),
        };
        return (sum, Vec::new());
    };
    let windows = windows_for_record(length_bp, sub_cand, cfg);
    let mut classified: Vec<(usize, usize, String)> = Vec::with_capacity(windows.len());
    let mut raw_meta: Vec<(usize, usize, i64, f64)> = Vec::with_capacity(windows.len());
    for &(s, e) in &windows {
        let (period, score) = kite_on_window(rec_id, seq, s, e);
        let cls = classify_window_period(period, score, sub_cand, cfg.tol, cfg.window_score_floor);
        classified.push((s, e, cls.to_string()));
        raw_meta.push((s, e, period, score));
    }
    let smoothed = smooth_runs(&classified, cfg.min_run);
    let mut window_rows: Vec<WindowRow> = Vec::with_capacity(windows.len());
    for (raw, smoothed) in raw_meta.iter().zip(smoothed.iter()) {
        window_rows.push(WindowRow {
            record_id: rec_id.to_string(),
            window_start: raw.0,
            window_end: raw.1,
            top_period: raw.2,
            top_score2_norm: raw.3,
            class_raw: classified
                .iter()
                .find(|c| c.0 == raw.0 && c.1 == raw.1)
                .map(|c| c.2.clone())
                .unwrap_or_default(),
            class_: smoothed.2.clone(),
        });
    }
    let n_sub = smoothed.iter().filter(|w| w.2 == "sub").count();
    let n_non_sub = smoothed.iter().filter(|w| w.2 == "non_sub").count();
    let blocks = blocks_from_windows(&smoothed);
    let flag = if !blocks.is_empty() && n_non_sub > 0 {
        "yes"
    } else {
        "no"
    };
    let reason = if flag == "yes" {
        "blocks+non_sub"
    } else if blocks.is_empty() {
        "no_blocks"
    } else {
        "no_non_sub_windows"
    };
    let cov_bp: usize = blocks.iter().map(|(s, e)| e - s).sum();
    let cov_pct = if length_bp > 0 {
        ((100.0 * cov_bp as f64 / length_bp as f64) * 100.0).round() / 100.0
    } else {
        0.0
    };
    let blocks_str = if blocks.is_empty() {
        "NA".to_string()
    } else {
        blocks
            .iter()
            .map(|(s, e)| format!("{s}-{e}"))
            .collect::<Vec<_>>()
            .join(";")
    };
    let sum = SummaryRow {
        record_id: rec_id.to_string(),
        length_bp,
        host_period_bp: (host_cand.round() as i64).to_string(),
        subrepeat_period_bp: (sub_cand.round() as i64).to_string(),
        subrepeat_flag: flag.into(),
        reason: reason.into(),
        n_windows_total: windows.len(),
        n_windows_sub: n_sub,
        n_windows_non_sub: n_non_sub,
        n_subrepeat_blocks: blocks.len(),
        subrepeat_coverage_bp: cov_bp,
        subrepeat_coverage_pct: cov_pct,
        blocks: blocks_str,
    };
    (sum, window_rows)
}

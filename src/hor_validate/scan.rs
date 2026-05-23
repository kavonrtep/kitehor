//! Per-record scan: within-tile peaks + spatial density windows + phase
//! folding + decision. Mirrors `hor_within_tile_check.py` exactly.

use super::Config;
use crate::kite::{analyze as kite_analyze, KiteConfig};
use crate::sequence::ArrayRecord;
use crate::subrepeat::scan::PeakRow as KitePeakRow;
use ahash::AHashMap;

#[derive(Debug, Clone)]
pub struct HorVerdict {
    pub case_id: String,
    pub founder: f64,
    pub tile: f64,
    pub multiplicity: Option<u32>,
}

#[derive(Debug, Clone)]
pub struct ValidationRow {
    pub record_id: String,
    pub global_founder_bp: f64,
    pub global_tile_bp: i64,
    pub global_founder_score: f64,
    pub global_tile_score: f64,
    /// NaN when tile_score == 0.
    pub global_founder_tile_ratio: f64,
    /// "NA" or integer string.
    pub within_top_period: String,
    pub within_top_score: f64,
    pub within_founder_score: f64,
    pub within_founder_top_ratio: f64,
    pub decision_hint: String,
    /// "NA" or float (rounded to 4 decimals).
    pub founder_density: String,
    pub phase_contrast: String,
    pub density_n_windows: u64,
    pub density_hint: String,
    pub skip_reason: String,
}

fn sum_near(peaks: &[KitePeakRow], target: f64, tol: f64) -> f64 {
    if target <= 0.0 {
        return 0.0;
    }
    peaks
        .iter()
        .filter(|p| ((p.period as f64) - target).abs() / target <= tol)
        .map(|p| p.score2_norm)
        .sum()
}

fn max_score(peaks: &[KitePeakRow]) -> f64 {
    peaks
        .iter()
        .map(|p| p.score2_norm)
        .fold(0.0f64, |a, b| a.max(b))
}

fn round_n(x: f64, n: i32) -> f64 {
    let p = 10f64.powi(n);
    (x * p).round() / p
}

fn kite_on_window(win_id: String, seq_slice: &[u8]) -> Vec<KitePeakRow> {
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
        .map(|(i, p)| KitePeakRow {
            rank: (i + 1) as u32,
            period: p.period,
            score2_norm: p.score2_norm,
        })
        .collect()
}

/// Run the validation over `verdicts`.
pub fn run(
    records: &[(String, Vec<u8>)],
    verdicts: &[HorVerdict],
    global: &AHashMap<String, Vec<KitePeakRow>>,
    cfg: &Config,
) -> Vec<ValidationRow> {
    let rec_by_id: AHashMap<&str, &[u8]> = records
        .iter()
        .map(|(id, s)| (id.as_str(), s.as_slice()))
        .collect();
    let empty_peaks: Vec<KitePeakRow> = Vec::new();
    let mut rows: Vec<ValidationRow> = Vec::with_capacity(verdicts.len());

    for v in verdicts {
        let rec_id = v.case_id.clone();
        let founder = v.founder;
        let tile_f = v.tile;
        let tile_int = tile_f.round() as i64;

        let mut row = ValidationRow {
            record_id: rec_id.clone(),
            global_founder_bp: round_n(founder, 2),
            global_tile_bp: tile_int,
            global_founder_score: 0.0,
            global_tile_score: 0.0,
            global_founder_tile_ratio: f64::NAN,
            within_top_period: "NA".into(),
            within_top_score: 0.0,
            within_founder_score: 0.0,
            within_founder_top_ratio: f64::NAN,
            decision_hint: "NA".into(),
            founder_density: "NA".into(),
            phase_contrast: "NA".into(),
            density_n_windows: 0,
            density_hint: "NA".into(),
            skip_reason: "".into(),
        };

        let g = global.get(&rec_id).unwrap_or(&empty_peaks);
        row.global_founder_score = sum_near(g, founder, cfg.period_match_tol);
        row.global_tile_score = sum_near(g, tile_f, cfg.period_match_tol);
        if row.global_tile_score > 0.0 {
            row.global_founder_tile_ratio = row.global_founder_score / row.global_tile_score;
        }

        let Some(seq) = rec_by_id.get(rec_id.as_str()) else {
            rows.push(row);
            continue;
        };

        // Within-tile (first tile) check.
        let mut skip_reason = "".to_string();
        if tile_int < (cfg.min_window_bp as i64) || tile_int > (cfg.max_tile_bp as i64) {
            skip_reason = "tile_out_of_range".into();
        } else if (tile_int as usize) > seq.len() {
            skip_reason = "tile_exceeds_array".into();
        }
        if !skip_reason.is_empty() {
            row.skip_reason = skip_reason;
            apply_density_or_skip(&mut row, v, seq, global, cfg);
            rows.push(row);
            continue;
        }

        let win_id = format!("{rec_id}__TILE0_{tile_int}");
        let win_peaks = kite_on_window(win_id, &seq[..tile_int as usize]);

        if win_peaks.is_empty() {
            row.skip_reason = "kite_no_peaks_in_window".into();
            apply_density_or_skip(&mut row, v, seq, global, cfg);
            rows.push(row);
            continue;
        }

        let top = win_peaks
            .iter()
            .max_by(|a, b| {
                a.score2_norm
                    .partial_cmp(&b.score2_norm)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        row.within_top_period = top.period.to_string();
        row.within_top_score = top.score2_norm;
        row.within_founder_score = sum_near(&win_peaks, founder, cfg.period_match_tol);
        if row.within_top_score > 0.0 {
            row.within_founder_top_ratio = row.within_founder_score / row.within_top_score;
        }
        row.decision_hint = if row.within_founder_top_ratio.is_nan() {
            "NA".into()
        } else {
            let r = row.within_founder_top_ratio;
            if r >= 0.5 {
                "strongly_confirms_hor".into()
            } else if r >= 0.2 {
                "weakly_confirms_hor".into()
            } else if r >= 0.05 {
                "ambiguous".into()
            } else {
                "suggests_within_monomer_duplication".into()
            }
        };

        apply_density_or_skip(&mut row, v, seq, global, cfg);
        rows.push(row);
    }

    rows
}

/// Compute the spatial density + phase contrast for one row. Mutates
/// `row` in place. Handles the `k <= MIN_K_FOR_DENSITY - 1` skip.
fn apply_density_or_skip(
    row: &mut ValidationRow,
    v: &HorVerdict,
    seq: &[u8],
    _global: &AHashMap<String, Vec<KitePeakRow>>,
    cfg: &Config,
) {
    let founder = v.founder;
    let tile = v.tile;
    if let Some(k) = v.multiplicity {
        if k < cfg.min_k_for_density {
            row.founder_density = "NA".into();
            row.density_n_windows = 0;
            row.phase_contrast = "NA".into();
            row.density_hint = format!("k_too_low_for_test(k={k})");
            return;
        }
    }
    let window_bp_f = (tile / cfg.density_window_tile_frac as f64)
        .max(cfg.min_founder_mult as f64 * founder)
        .max(cfg.min_density_window_bp as f64);
    let window_bp = window_bp_f.round() as usize;
    if window_bp > seq.len() {
        return;
    }
    let base_step = (founder.round() as usize).max(1);
    let max_step = ((seq.len() - window_bp) / cfg.max_density_windows + 1).max(1);
    let step = base_step.max(max_step);

    let mut n_total: u64 = 0;
    let mut n_founder: u64 = 0;
    let mut phase_total = vec![0u64; cfg.phase_fold_bins];
    let mut phase_present = vec![0u64; cfg.phase_fold_bins];
    let bin_width = ((tile / cfg.phase_fold_bins as f64).round() as usize).max(1);
    let tile_round = tile.round() as usize;

    let mut s = 0usize;
    while s + window_bp <= seq.len() {
        n_total += 1;
        let win_mid = (s + s + window_bp) / 2;
        let mut phase_bin = if tile_round > 0 {
            (win_mid % tile_round) / bin_width
        } else {
            0
        };
        if phase_bin >= cfg.phase_fold_bins {
            phase_bin = cfg.phase_fold_bins - 1;
        }
        phase_total[phase_bin] += 1;

        let win_id = format!("{}__SP_{}_{}", v.case_id, s, s + window_bp);
        let win_peaks = kite_on_window(win_id, &seq[s..s + window_bp]);
        if !win_peaks.is_empty() {
            let top_score = max_score(&win_peaks);
            let f_score = sum_near(&win_peaks, founder, cfg.period_match_tol);
            if top_score > 0.0 && f_score / top_score >= cfg.density_rel_floor {
                n_founder += 1;
                phase_present[phase_bin] += 1;
            }
        }
        s += step;
    }
    if n_total == 0 {
        return;
    }
    let density = n_founder as f64 / n_total as f64;
    row.founder_density = fmt_round4(density);
    row.density_n_windows = n_total;
    let frac_per_bin: Vec<f64> = (0..cfg.phase_fold_bins)
        .filter(|i| phase_total[*i] > 0)
        .map(|i| phase_present[i] as f64 / phase_total[i] as f64)
        .collect();
    if frac_per_bin.len() >= 2 {
        let max_f = frac_per_bin
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_f = frac_per_bin.iter().cloned().fold(f64::INFINITY, f64::min);
        let contrast = max_f - min_f;
        row.phase_contrast = fmt_round4(contrast);
        let is_dup = density <= cfg.density_dup_max || contrast >= cfg.phase_contrast_dup_min;
        let is_hor = density >= cfg.density_hor_min && contrast <= cfg.phase_contrast_hor_max;
        row.density_hint = if is_dup {
            "localized_duplication".into()
        } else if is_hor {
            "spatially_confirms_hor".into()
        } else {
            "ambiguous".into()
        };
    } else {
        row.phase_contrast = "NA".into();
        row.density_hint = "insufficient_phase_bins".into();
    }
}

fn fmt_round4(x: f64) -> String {
    // The prototype calls `round(x, 4)` and emits via pandas's
    // `to_csv(..., float_format="%.6g")` — i.e. round first, then format.
    let rounded = round_n(x, 4);
    crate::rule_classify::io::fmt_g(6, rounded)
}

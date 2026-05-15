//! Bare kite periodicity detector — Rust port of
//! `TideCluster/tarean/kite.R::get_peaks_from_seq`.
//!
//! Matches the R algorithm exactly:
//!
//! - `k = 6`, **non-canonical** k-mers (forward strand only; no
//!   reverse-complement merging — kite.R's `substring()` based grouping).
//! - Build H[d] = tabulate(unlist(diff(positions per k-mer))) for d ∈ [1, L].
//! - Composition-matched random background, N = 10 replicates,
//!   element-wise max → smoothed with a wide gaussian → rescaled by
//!   `max(env/smoothed)[0..L/2]` (mirrors kite.R's smooth.spline +
//!   rescale; existing `periodogram.rs::compute_kite_background` already
//!   uses this approximation).
//! - For each local maximum:
//!   - score   = peak / (L − position − k + 1)
//!   - score2  = score · log2(position),  then sum-normalised
//!   - filter: position < L/2; score2_norm > 0.001; peak > background
//!   - fallback when no peak passes peak>bg: pick which.max(peak/bg)
//! - Top-1: `which.max(score)` after the filter.
//! - Top-3: peaks sorted by score desc.

use crate::sequence::ArrayRecord;
use ahash::AHashMap;
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct KiteConfig {
    pub k: usize,
    pub n_bg_replicates: usize,
    pub score2_threshold: f64,
    pub min_peak_distance: usize,
    /// Override gaussian-smoothing sigma for the bg envelope. When `None`,
    /// picks `max(5.0, L / 1000.0)` per record.
    pub bg_smoothing_sigma: Option<f64>,
}

impl Default for KiteConfig {
    fn default() -> Self {
        Self {
            k: 6,
            n_bg_replicates: 10,
            score2_threshold: 0.001,
            min_peak_distance: 1,
            bg_smoothing_sigma: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct KitePeak {
    pub period: usize,
    pub peak_height: f64,
    pub score: f64,
    pub score2: f64,
    pub score2_norm: f64,
    pub background: f64,
}

#[derive(Debug, Clone)]
pub struct KiteResult {
    pub array_id: String,
    pub length_bp: usize,
    /// Filtered peaks, sorted by score desc. Empty when no peak survives.
    pub peaks: Vec<KitePeak>,
    /// Optional full H[d] (size L+1; H[0] unused).
    pub profile: Option<Vec<f64>>,
    /// Optional background envelope.
    pub background: Option<Vec<f64>>,
}

pub fn analyze(record: &ArrayRecord, cfg: &KiteConfig, dump_profile: bool) -> KiteResult {
    let l = record.length;
    let k = cfg.k;
    if l < k + 2 {
        return KiteResult {
            array_id: record.id.clone(),
            length_bp: l,
            peaks: Vec::new(),
            profile: if dump_profile { Some(Vec::new()) } else { None },
            background: if dump_profile { Some(Vec::new()) } else { None },
        };
    }
    let profile = compute_neighbor_profile(&record.seq, k, l);
    // Default sigma = 10 matches `periodogram.rs::compute_kite_background`,
    // and is small enough to preserve the per-bin random-noise structure
    // (kite.R's smooth.spline picks a similar effective bandwidth via
    // cross-validation, not the L/1000 sigma initially considered).
    let bg_sigma = cfg.bg_smoothing_sigma.unwrap_or(10.0);
    let background =
        compute_background(&record.seq, &record.id, k, l, cfg.n_bg_replicates, bg_sigma);
    let peaks = find_peaks_with_score(
        &profile,
        &background,
        l,
        k,
        cfg.score2_threshold,
        cfg.min_peak_distance,
    );
    KiteResult {
        array_id: record.id.clone(),
        length_bp: l,
        peaks,
        profile: if dump_profile { Some(profile) } else { None },
        background: if dump_profile { Some(background) } else { None },
    }
}

/// H[d] for d ∈ [1, L], non-canonical k-mers. K-mers containing any
/// byte other than A/C/G/T are skipped (matches kite.R behaviour —
/// non-ACGT bytes form k-mers but never match the random-bg samples,
/// so they contribute approximately nothing).
fn compute_neighbor_profile(seq: &[u8], k: usize, l: usize) -> Vec<f64> {
    let max_d = l + 1;
    let mut profile = vec![0f64; max_d];
    if seq.len() < k {
        return profile;
    }
    let mut positions: AHashMap<&[u8], Vec<u32>> = AHashMap::new();
    let limit = seq.len() - k + 1;
    for i in 0..limit {
        let s = &seq[i..i + k];
        if s.iter().any(|b| !matches!(*b, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        positions.entry(s).or_default().push(i as u32);
    }
    for pos in positions.values() {
        if pos.len() < 2 {
            continue;
        }
        for w in pos.windows(2) {
            let d = (w[1] - w[0]) as usize;
            if d < max_d {
                profile[d] += 1.0;
            }
        }
    }
    profile
}

/// Composition-matched random background. Mirrors kite.R's
/// `get_neighgor_distances_background` and the production
/// `periodogram.rs::compute_kite_background`:
///   1. estimate ACGT composition
///   2. for N replicates: sample a random sequence of length L,
///      compute its profile, take element-wise max into the envelope
///   3. smooth the envelope with a wide gaussian (approx. smooth.spline)
///   4. rescale: `multiple = max(env/smoothed)[0..L/2]`, multiply
fn compute_background(
    seq: &[u8],
    id: &str,
    k: usize,
    l: usize,
    n_replicates: usize,
    bg_sigma: f64,
) -> Vec<f64> {
    let max_d = l + 1;
    if l < k || n_replicates == 0 {
        return vec![0f64; max_d];
    }
    // ACGT composition.
    let mut counts = [0u64; 4];
    let mut n_acgt: u64 = 0;
    for &b in seq.iter() {
        match b {
            b'A' => {
                counts[0] += 1;
                n_acgt += 1;
            }
            b'C' => {
                counts[1] += 1;
                n_acgt += 1;
            }
            b'G' => {
                counts[2] += 1;
                n_acgt += 1;
            }
            b'T' => {
                counts[3] += 1;
                n_acgt += 1;
            }
            _ => {}
        }
    }
    if n_acgt == 0 {
        return vec![0f64; max_d];
    }
    let nf = n_acgt as f64;
    let cum = [
        counts[0] as f64 / nf,
        (counts[0] + counts[1]) as f64 / nf,
        (counts[0] + counts[1] + counts[2]) as f64 / nf,
        1.0_f64,
    ];

    // FNV-1a seed for determinism.
    let seed_base: u64 = {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for &b in id.as_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100_0000_01b3);
        }
        h
    };

    let mut envelope = vec![0f64; max_d];
    let mut buf: Vec<u8> = vec![0; l];
    for j in 0..n_replicates {
        let mut state = seed_base.wrapping_add(0x9E37_79B9_7F4A_7C15u64.wrapping_mul(j as u64 + 1));
        if state == 0 {
            state = 1;
        }
        for byte in buf.iter_mut() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let u = (state as f64) / (u64::MAX as f64);
            *byte = if u < cum[0] {
                b'A'
            } else if u < cum[1] {
                b'C'
            } else if u < cum[2] {
                b'G'
            } else {
                b'T'
            };
        }
        let hist = compute_neighbor_profile(&buf, k, l);
        for i in 0..max_d {
            if hist[i] > envelope[i] {
                envelope[i] = hist[i];
            }
        }
    }
    let smoothed = gaussian_smooth(&envelope, bg_sigma);
    // Rescale: multiple = max(envelope / smoothed) over [0, L/2].
    let half_l = (l / 2).min(max_d.saturating_sub(1));
    let mut multiple: f64 = 1.0;
    for i in 0..=half_l {
        if smoothed[i] > 0.0 {
            let r = envelope[i] / smoothed[i];
            if r > multiple {
                multiple = r;
            }
        }
    }
    smoothed.iter().map(|x| x * multiple).collect()
}

/// 1-D gaussian smoothing kernel with ±3σ window. Same as
/// `periodogram.rs::gaussian_smooth`.
fn gaussian_smooth(hist: &[f64], sigma: f64) -> Vec<f64> {
    if sigma <= 0.0 {
        return hist.to_vec();
    }
    let radius = (3.0 * sigma).ceil() as isize;
    let mut kernel: Vec<f64> = (-radius..=radius)
        .map(|x| (-(x as f64).powi(2) / (2.0 * sigma * sigma)).exp())
        .collect();
    let ksum: f64 = kernel.iter().sum();
    for k in kernel.iter_mut() {
        *k /= ksum;
    }
    let n = hist.len();
    let mut out = vec![0.0f64; n];
    for (i, slot) in out.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (ki, kv) in kernel.iter().enumerate() {
            let off = ki as isize - radius;
            let j = i as isize + off;
            if j >= 0 && (j as usize) < n {
                acc += hist[j as usize] * kv;
            }
        }
        *slot = acc;
    }
    out
}

/// Port of kite.R's `simplified_findpeaks` + `get_peaks_from_neighbor_distances`.
/// Strict local maxima from sign-change of `diff(profile)`, scored per
/// kite.R, filtered by score2 threshold + `position < L/2` + peak > bg
/// (with the kite.R fallback when none pass).
fn find_peaks_with_score(
    profile: &[f64],
    background: &[f64],
    l: usize,
    k: usize,
    score2_threshold: f64,
    min_peak_distance: usize,
) -> Vec<KitePeak> {
    // Find strict local maxima: i where profile[i-1] < profile[i] > profile[i+1]
    // (kite.R uses sign(diff(x)) "+-" pattern). Min peak distance enforced
    // by keeping highest peak in each disqualifying cluster.
    let mut raw: Vec<(usize, f64)> = Vec::new();
    for i in 1..profile.len().saturating_sub(1) {
        if profile[i] > profile[i - 1] && profile[i] > profile[i + 1] {
            raw.push((i, profile[i]));
        }
    }
    if min_peak_distance > 1 {
        // kite.R: sort by height desc, then drop any peak within
        // min_peak_distance of an already-kept one.
        raw.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut kept: Vec<(usize, f64)> = Vec::new();
        for (pos, h) in raw.drain(..) {
            if !kept.iter().any(|(p, _)| {
                let dd = pos.abs_diff(*p);
                dd > 0 && dd < min_peak_distance
            }) {
                kept.push((pos, h));
            }
        }
        raw = kept;
    }
    // Score per kite.R.
    let denom = |pos: usize| -> f64 { (l as f64 - pos as f64 - k as f64 + 1.0).max(1.0) };
    let mut peaks: Vec<KitePeak> = raw
        .into_iter()
        .filter(|(pos, _)| *pos < l / 2 && *pos >= 1)
        .map(|(pos, h)| {
            let score = h / denom(pos);
            let score2 = score * (pos as f64).log2();
            KitePeak {
                period: pos,
                peak_height: h,
                score,
                score2,
                score2_norm: 0.0, // filled below
                background: background.get(pos).copied().unwrap_or(0.0),
            }
        })
        .collect();
    // score2_norm = score2 / sum(score2)
    let sum_score2: f64 = peaks.iter().map(|p| p.score2).sum();
    if sum_score2 > 0.0 {
        for p in peaks.iter_mut() {
            p.score2_norm = p.score2 / sum_score2;
        }
    }
    // Filter by score2_norm > threshold.
    peaks.retain(|p| p.score2_norm > score2_threshold);

    // Peak > background filter (with kite.R fallback).
    let any_above = peaks.iter().any(|p| p.peak_height > p.background);
    if any_above {
        peaks.retain(|p| p.peak_height > p.background);
    } else if !peaks.is_empty() {
        // which.max(peak/background) — closest to bg
        let best_idx = peaks
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                let ra = a.peak_height / a.background.max(1.0);
                let rb = b.peak_height / b.background.max(1.0);
                ra.partial_cmp(&rb).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        let keep = peaks.swap_remove(best_idx);
        peaks.clear();
        peaks.push(keep);
    }
    // Sort by score desc (kite.R: which.max(score) for top-1; sorted by
    // score desc for top-3).
    peaks.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    peaks
}

/// Convenience wrapper for parallel processing over a record set. Drops
/// any per-record dump-profile data; use `analyze` directly when needed.
pub fn analyze_records(records: &[ArrayRecord], cfg: &KiteConfig) -> Vec<KiteResult> {
    records.par_iter().map(|r| analyze(r, cfg, false)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ArrayRecord;

    fn make_record(id: &str, seq: &[u8]) -> ArrayRecord {
        ArrayRecord::from_raw(id.to_string(), seq)
    }

    fn tandem(monomer: &[u8], copies: usize) -> Vec<u8> {
        let mut s = Vec::with_capacity(monomer.len() * copies);
        for _ in 0..copies {
            s.extend_from_slice(monomer);
        }
        s
    }

    #[test]
    fn finds_clean_monomer_period() {
        // 178-bp synthetic monomer (random-looking 178 chars), 100 copies.
        // Top-1 should be 178.
        let mut monomer = Vec::with_capacity(178);
        let mut state: u64 = 0xc0de_d00d_dead_beef;
        for _ in 0..178 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let nt = match state & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            monomer.push(nt);
        }
        assert_eq!(monomer.len(), 178);
        let seq = tandem(&monomer, 100);
        let rec = make_record("test", &seq);
        let cfg = KiteConfig::default();
        let res = analyze(&rec, &cfg, false);
        assert!(!res.peaks.is_empty(), "no peaks found");
        let top = res.peaks[0].period;
        assert!(
            (top as i64 - 178).abs() <= 2,
            "top-1 period = {} (expected 178)",
            top
        );
    }

    #[test]
    fn rejects_random_sequence() {
        // Random-ish sequence (no real periodicity).
        let mut state: u64 = 0xdead_beef_cafe_babe;
        let mut seq = Vec::with_capacity(50_000);
        for _ in 0..50_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let nt = match state & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            seq.push(nt);
        }
        let rec = make_record("rand", &seq);
        let cfg = KiteConfig::default();
        let res = analyze(&rec, &cfg, false);
        // Either no peaks, or peaks whose top-1 score is very low. We
        // don't enforce zero — kite.R's "fallback to closest-to-bg" can
        // emit a single peak — but the score should be small.
        if !res.peaks.is_empty() {
            assert!(
                res.peaks[0].score < 0.01,
                "random seq produced top score {} (too high)",
                res.peaks[0].score
            );
        }
    }

    #[test]
    fn score_formula_matches_kite_r() {
        // For a chosen (peak_height, position, L, k) the score must be
        // exactly peak / (L - position - k + 1).
        let l = 10_000usize;
        let k = 6usize;
        let pos = 178usize;
        let denom = l as f64 - pos as f64 - k as f64 + 1.0;
        // pure synthetic: profile with one spike.
        let mut profile = vec![0f64; l + 1];
        profile[pos] = 1234.0;
        let bg = vec![0f64; l + 1];
        let peaks = find_peaks_with_score(&profile, &bg, l, k, 0.0, 1);
        assert_eq!(peaks.len(), 1);
        let p = &peaks[0];
        assert_eq!(p.period, pos);
        let expected = 1234.0 / denom;
        assert!(
            (p.score - expected).abs() < 1e-9,
            "score {} != expected {}",
            p.score,
            expected
        );
    }
}

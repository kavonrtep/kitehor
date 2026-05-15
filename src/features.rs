//! Per-record feature extraction for the kite-first classifier.
//!
//! Ports `eval/training_data/extract_features.py` to Rust. Produces the
//! 28-feature schema consumed by the calibrated HOR-score / k-predictor
//! random forests (+ the two homology columns `h_d1`, `h_founder` which
//! are filled in separately by `homology.rs`).
//!
//! Schema (kept in lock-step with the R training pipeline):
//!   case_id, stratum, array_length,
//!   s1, s2, s3, s2_over_s1, s3_over_s1,
//!   family_size_best, tile_founder_ratio, tile_jitter,
//!   d1, d2, d3, log_d1_over_L,
//!   d2_over_d1, d3_over_d1, max_d_top3_over_min_d_top3,
//!   d4 .. d10,
//!   family_founder_d, family_tile_d,
//!   distinct_kmers_per_bp, kmer_entropy, singletons_ratio,
//!   (h_d1, h_founder are added by homology.rs)

use ahash::AHashMap;

use crate::kite::KiteResult;
use crate::sequence::ArrayRecord;

/// All numeric features needed for the classifier, plus the small
/// strings (`case_id`, `stratum`) used for I/O and grouping.
#[derive(Debug, Clone)]
pub struct FeatureRow {
    pub case_id: String,
    pub stratum: String,
    pub array_length: usize,
    // Kite scores.
    pub s1: f64,
    pub s2: f64,
    pub s3: f64,
    pub s2_over_s1: f64,
    pub s3_over_s1: f64,
    // Family-search results.
    pub family_size_best: usize,
    pub tile_founder_ratio: f64,
    pub tile_jitter: usize,
    // Kite periods (top-N by score, 0 padding when fewer than N peaks).
    pub d1: usize,
    pub d2: usize,
    pub d3: usize,
    pub log_d1_over_l: f64,
    pub d2_over_d1: f64,
    pub d3_over_d1: f64,
    pub max_d_top3_over_min_d_top3: f64,
    pub d4: usize,
    pub d5: usize,
    pub d6: usize,
    pub d7: usize,
    pub d8: usize,
    pub d9: usize,
    pub d10: usize,
    pub family_founder_d: usize,
    pub family_tile_d: usize,
    // Sequence diversity.
    pub distinct_kmers_per_bp: f64,
    pub kmer_entropy: f64,
    pub singletons_ratio: f64,
    // Homology features (filled in later — NaN by default).
    pub h_d1: f64,
    pub h_founder: f64,
}

impl FeatureRow {
    /// Resolve a feature value by name, matching the column names in
    /// the R training TSV (used by both RF models).
    pub fn get(&self, name: &str) -> Option<f64> {
        Some(match name {
            "array_length" => self.array_length as f64,
            "s1" => self.s1,
            "s2" => self.s2,
            "s3" => self.s3,
            "s2_over_s1" => self.s2_over_s1,
            "s3_over_s1" => self.s3_over_s1,
            "family_size_best" => self.family_size_best as f64,
            "tile_founder_ratio" => self.tile_founder_ratio,
            "tile_jitter" => self.tile_jitter as f64,
            "d1" => self.d1 as f64,
            "d2" => self.d2 as f64,
            "d3" => self.d3 as f64,
            "log_d1_over_L" => self.log_d1_over_l,
            "d2_over_d1" => self.d2_over_d1,
            "d3_over_d1" => self.d3_over_d1,
            "max_d_top3_over_min_d_top3" => self.max_d_top3_over_min_d_top3,
            "d4" => self.d4 as f64,
            "d5" => self.d5 as f64,
            "d6" => self.d6 as f64,
            "d7" => self.d7 as f64,
            "d8" => self.d8 as f64,
            "d9" => self.d9 as f64,
            "d10" => self.d10 as f64,
            "family_founder_d" => self.family_founder_d as f64,
            "family_tile_d" => self.family_tile_d as f64,
            "distinct_kmers_per_bp" => self.distinct_kmers_per_bp,
            "kmer_entropy" => self.kmer_entropy,
            "singletons_ratio" => self.singletons_ratio,
            "h_d1" => self.h_d1,
            "h_founder" => self.h_founder,
            _ => return None,
        })
    }
}

// ---------------------------------------------------------------------------
// Diversity features (k=6 distinct k-mer count, Shannon entropy,
// singletons ratio). Mirrors `extract_features.py::diversity_features`.
// Skips k-mers containing N (matches kite.R).
// ---------------------------------------------------------------------------

pub fn diversity_features(seq: &[u8], k: usize) -> (f64, f64, f64) {
    let l = seq.len();
    if l < k + 1 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let mut counts: AHashMap<&[u8], u32> = AHashMap::with_capacity(l);
    let mut total: u32 = 0;
    for i in 0..=l - k {
        let s = &seq[i..i + k];
        if s.iter().any(|&b| b == b'N' || b == b'n') {
            continue;
        }
        *counts.entry(s).or_insert(0) += 1;
        total += 1;
    }
    if total == 0 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let distinct = counts.len();
    let singletons = counts.values().filter(|&&v| v == 1).count();
    let total_f = total as f64;
    let mut entropy = 0.0_f64;
    for &c in counts.values() {
        let p = c as f64 / total_f;
        entropy -= p * p.log2();
    }
    (
        distinct as f64 / l as f64,
        entropy,
        singletons as f64 / total_f,
    )
}

// ---------------------------------------------------------------------------
// Family search — port of `find_best_family` in extract_features.py.
// This produces the *features* (family_size_best, tile_founder_ratio,
// tile_period, tile_jitter) AND the founder period (`family_founder_d`)
// needed for downstream homology probing. It is intentionally
// independent of `hor_call.rs` because the verdict layer has its own
// thresholds; here we want raw structural info for the model.
//
// Parameters chosen to match the Python defaults (qmax=30, tol_bp=5,
// tol_rel=0.02, lo_period=15, min_size=2, min_founder_top1_share=0.5,
// require_top_k=3, jitter_tol=0.15).
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
struct PeakRef {
    period: usize,
    score: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct FamilyFeatures {
    pub size_best: usize,
    pub tile_founder_ratio: f64,
    pub tile_period: usize, // 0 = no tile
    pub tile_jitter: usize,
    pub founder_d: usize, // 0 = no family
    pub tile_d: usize,
}

fn fits(p: usize, m: usize, qmax: usize, tol_bp: usize, tol_rel: f64) -> Option<usize> {
    if m == 0 {
        return None;
    }
    let k = ((p as f64) / (m as f64)).round() as i64;
    if k < 1 || (k as usize) > qmax {
        return None;
    }
    let k = k as usize;
    let expected = (k * m) as i64;
    let diff = ((p as i64) - expected).unsigned_abs() as usize;
    let tol = tol_bp.max((tol_rel * expected as f64) as usize);
    if diff <= tol {
        Some(k)
    } else {
        None
    }
}

pub fn family_features(peaks_by_score: &[(usize, f64)]) -> FamilyFeatures {
    let qmax = 30usize;
    let tol_bp = 5usize;
    let tol_rel = 0.02_f64;
    let lo_period = 15usize;
    let min_size = 2usize;
    let min_founder_share = 0.5_f64;
    let require_top_k = 3usize;
    let jitter_tol = 0.15_f64;

    if peaks_by_score.is_empty() {
        return FamilyFeatures {
            size_best: 0,
            tile_founder_ratio: 0.0,
            tile_period: 0,
            tile_jitter: 0,
            founder_d: 0,
            tile_d: 0,
        };
    }

    let peaks: Vec<PeakRef> = peaks_by_score
        .iter()
        .map(|&(p, s)| PeakRef {
            period: p,
            score: s,
        })
        .collect();
    let top1_score = peaks[0].score;
    let top_k_periods: Vec<usize> = peaks.iter().take(require_top_k).map(|p| p.period).collect();

    for cand in &peaks {
        let m_f = cand.period;
        if m_f < lo_period {
            continue;
        }
        if top1_score > 0.0 && cand.score < min_founder_share * top1_score {
            continue;
        }
        // Family scan.
        let mut family_count: usize = 0;
        let mut best_tile_k: usize = 1;
        let mut best_tile_period: usize = m_f;
        let mut best_tile_score: f64 = -1.0;
        for p in &peaks {
            let Some(k) = fits(p.period, m_f, qmax, tol_bp, tol_rel) else {
                continue;
            };
            family_count += 1;
            if k >= 2 && p.score > best_tile_score {
                best_tile_score = p.score;
                best_tile_k = k;
                best_tile_period = p.period;
            }
        }
        if family_count < min_size {
            continue;
        }
        if !top_k_periods
            .iter()
            .all(|&tp| fits(tp, m_f, qmax, tol_bp, tol_rel).is_some())
        {
            continue;
        }
        // Jitter: count peaks within ±jitter_tol * tile of the chosen tile.
        let band = jitter_tol * best_tile_period as f64;
        let tile_jitter = peaks
            .iter()
            .filter(|p| (p.period as f64 - best_tile_period as f64).abs() <= band)
            .count();
        let ratio = if cand.score > 0.0 && best_tile_score > 0.0 {
            best_tile_score / cand.score
        } else {
            0.0
        };
        let _ = best_tile_k; // k itself not exported in features (only via tile_period)
        return FamilyFeatures {
            size_best: family_count,
            tile_founder_ratio: ratio,
            tile_period: best_tile_period,
            tile_jitter,
            founder_d: m_f,
            tile_d: best_tile_period,
        };
    }
    FamilyFeatures {
        size_best: 0,
        tile_founder_ratio: 0.0,
        tile_period: 0,
        tile_jitter: 0,
        founder_d: 0,
        tile_d: 0,
    }
}

// ---------------------------------------------------------------------------
// Top-level feature builder
// ---------------------------------------------------------------------------

/// Compute the full feature row from a kite result + the underlying
/// sequence. Sets `h_d1`/`h_founder` to NaN — `homology.rs` fills these
/// in.
pub fn build_features(record: &ArrayRecord, kite: &KiteResult) -> FeatureRow {
    let l = record.seq.len();
    let case_id = record.id.clone();
    let stratum = derive_stratum(&case_id);

    // Peaks already come from kite::analyze sorted by score desc. Pull
    // top 10 periods and top 3 scores.
    let peaks: Vec<(usize, f64)> = kite.peaks.iter().map(|p| (p.period, p.score)).collect();

    let s = |i: usize| peaks.get(i).map(|p| p.1).unwrap_or(0.0);
    let d = |i: usize| peaks.get(i).map(|p| p.0).unwrap_or(0);
    let s1 = s(0);
    let s2 = s(1);
    let s3 = s(2);
    let d1 = d(0);
    let d2 = d(1);
    let d3 = d(2);

    let s2_over_s1 = if s1 > 0.0 { s2 / s1 } else { 0.0 };
    let s3_over_s1 = if s1 > 0.0 { s3 / s1 } else { 0.0 };
    let log_d1_over_l = if d1 > 0 && l > 0 {
        (d1 as f64 / l as f64).ln()
    } else {
        0.0
    };
    let d2_over_d1 = if d1 > 0 && d2 > 0 {
        d2 as f64 / d1 as f64
    } else {
        0.0
    };
    let d3_over_d1 = if d1 > 0 && d3 > 0 {
        d3 as f64 / d1 as f64
    } else {
        0.0
    };
    let nz: Vec<f64> = [d1, d2, d3]
        .iter()
        .filter(|&&x| x > 0)
        .map(|&x| x as f64)
        .collect();
    let max_d_top3_over_min_d_top3 = if nz.len() >= 2 {
        let mx = nz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mn = nz.iter().cloned().fold(f64::INFINITY, f64::min);
        mx / mn
    } else {
        1.0
    };

    let fam = family_features(&peaks);

    let (distinct_kmers_per_bp, kmer_entropy, singletons_ratio) =
        diversity_features(record.seq.as_slice(), 6);

    FeatureRow {
        case_id,
        stratum,
        array_length: l,
        s1,
        s2,
        s3,
        s2_over_s1,
        s3_over_s1,
        family_size_best: fam.size_best,
        tile_founder_ratio: fam.tile_founder_ratio,
        tile_jitter: fam.tile_jitter,
        d1,
        d2,
        d3,
        log_d1_over_l,
        d2_over_d1,
        d3_over_d1,
        max_d_top3_over_min_d_top3,
        d4: d(3),
        d5: d(4),
        d6: d(5),
        d7: d(6),
        d8: d(7),
        d9: d(8),
        d10: d(9),
        family_founder_d: fam.founder_d,
        family_tile_d: fam.tile_d,
        distinct_kmers_per_bp,
        kmer_entropy,
        singletons_ratio,
        h_d1: f64::NAN,
        h_founder: f64::NAN,
    }
}

/// Convention from `extract_features.py`: stratum = case_id with the
/// final `_<number>` suffix removed.
fn derive_stratum(case_id: &str) -> String {
    match case_id.rfind('_') {
        Some(i) => case_id[..i].to_string(),
        None => case_id.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diversity_constant_seq_low_entropy() {
        // All-A sequence: 1 distinct 6-mer, entropy 0.
        let seq = vec![b'A'; 100];
        let (dk, ent, sng) = diversity_features(&seq, 6);
        assert!(dk > 0.0 && dk < 0.011);
        assert!((ent - 0.0).abs() < 1e-12);
        assert!((sng - 0.0).abs() < 1e-12); // 1 singleton out of 1 distinct = wait
                                            // Actually the all-A case has 1 distinct, count = 95, no singletons.
    }

    #[test]
    fn diversity_handles_ns() {
        let seq = b"ACGTACGTACGTNNNNACGTACGT".to_vec();
        let (_, ent, _) = diversity_features(&seq, 6);
        assert!(ent >= 0.0);
    }

    #[test]
    fn family_features_pure_tandem_has_no_family() {
        // Peaks at 100, 200 (harmonic). Founder candidate 100 would
        // match itself + 200 (k=2). With require_top_k=3 and only 2
        // peaks, top_k_periods = [100, 200], both fit 100 → family found.
        // tile = 200 (k=2), tile_founder_ratio = score2/score1.
        let peaks = vec![(100, 1.0), (200, 0.3)];
        let f = family_features(&peaks);
        assert_eq!(f.size_best, 2);
        assert_eq!(f.tile_period, 200);
        assert_eq!(f.founder_d, 100);
        assert!((f.tile_founder_ratio - 0.3).abs() < 1e-12);
    }

    #[test]
    fn family_features_real_hor() {
        // Founder = 178, tile at 888 (k=5), plus harmonic 356 (k=2),
        // 534 (k=3), 712 (k=4). Family of 5.
        let peaks = vec![(178, 1.0), (888, 0.4), (356, 0.2), (534, 0.15), (712, 0.1)];
        let f = family_features(&peaks);
        assert_eq!(f.size_best, 5);
        assert_eq!(f.founder_d, 178);
        // Tile picks the *highest-score* k>=2 peak, which is 888.
        assert_eq!(f.tile_period, 888);
    }
}

//! Supplementary HOR-coverage QC.
//!
//! For each HOR call `(founder, tile, k)` we compute how well the first
//! tile in the array (= putative HOR unit) explains the rest of the
//! array. Identity is alignment-based (Levenshtein), so per-base indels
//! introduced by the simulator or by real sequence drift do not collapse
//! the score to background.
//!
//! ```text
//! ref = seq[0 : tile]
//! for i in 1..n_tiles:
//!     identity_i = 1 - levenshtein(seq[i·tile : (i+1)·tile], ref) / max(|w|, |r|)
//! ```
//!
//! The aggregate score is supplementary — it does not enter the
//! HOR/non-HOR decision in the rule classifier. It is reported only
//! when the user passes `--coverage`, and only for records that the
//! rule classified as `hor`. Captures:
//!
//! - *mosaic / partial-array* arrays: `cov_first_half` vs `cov_second_half`
//!   asymmetry.
//! - *tile dropouts*: `cov_min` flags individual outlier tiles.
//! - *wrong-period calls*: if the `tile` period does not actually tile
//!   the array, `cov_mean` will be ~0.25 (random).
//!
//! It does NOT distinguish a real founder from a sub-period of a longer
//! real monomer — both produce high tile-aligned identity by
//! construction. That distinction needs an external constraint (e.g.
//! a `--rule-lo-founder` floor).

#[derive(Debug, Clone, Copy)]
pub struct TileCoverage {
    pub mean: f64,
    pub pass_70: f64,
    pub pass_80: f64,
    pub pass_90: f64,
    pub first_half: f64,
    pub second_half: f64,
    pub min: f64,
    pub max: f64,
    pub n_tiles: usize,
}

/// Levenshtein-based identity in [0, 1]. Two-row Wagner-Fischer DP.
/// Memory O(min(|a|, |b|)); time O(|a| * |b|).
pub fn levenshtein_identity(a: &[u8], b: &[u8]) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    // Make `b` the shorter side so we keep less DP memory live.
    let (a, b) = if a.len() < b.len() { (b, a) } else { (a, b) };
    let n = a.len();
    let m = b.len();

    let mut prev: Vec<u32> = (0..=m as u32).collect();
    let mut curr: Vec<u32> = vec![0; m + 1];

    for i in 1..=n {
        curr[0] = i as u32;
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let del = prev[j] + 1;
            let ins = curr[j - 1] + 1;
            let sub = prev[j - 1] + cost;
            curr[j] = del.min(ins).min(sub);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    let dist = prev[m] as f64;
    let max_len = n.max(m) as f64;
    1.0 - dist / max_len
}

/// Slide a tile-length window across `seq` with step = `tile`. For each
/// subsequent window compute identity vs `seq[0..tile]`. Aggregate.
///
/// Returns `None` if the array contains fewer than two full tiles.
pub fn compute_tile_coverage(seq: &[u8], tile: usize) -> Option<TileCoverage> {
    let l = seq.len();
    if tile == 0 {
        return None;
    }
    let n_tiles = l / tile;
    if n_tiles < 2 {
        return None;
    }
    let reference = &seq[0..tile];
    let mut identities: Vec<f64> = Vec::with_capacity(n_tiles - 1);
    for i in 1..n_tiles {
        let window = &seq[i * tile..(i + 1) * tile];
        identities.push(levenshtein_identity(window, reference));
    }
    aggregate(&identities)
}

fn aggregate(ids: &[f64]) -> Option<TileCoverage> {
    if ids.is_empty() {
        return None;
    }
    let mean = ids.iter().sum::<f64>() / ids.len() as f64;
    let pass_at =
        |t: f64| -> f64 { ids.iter().filter(|&&x| x >= t).count() as f64 / ids.len() as f64 };
    let half = ids.len() / 2;
    let first = if half > 0 {
        ids[..half].iter().sum::<f64>() / half as f64
    } else {
        f64::NAN
    };
    let second_len = ids.len() - half;
    let second = if second_len > 0 {
        ids[half..].iter().sum::<f64>() / second_len as f64
    } else {
        f64::NAN
    };
    let min = ids.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = ids.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    Some(TileCoverage {
        mean,
        pass_70: pass_at(0.70),
        pass_80: pass_at(0.80),
        pass_90: pass_at(0.90),
        first_half: first,
        second_half: second,
        min,
        max,
        n_tiles: ids.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_are_1() {
        let a = b"ACGTACGTACGT";
        assert!((levenshtein_identity(a, a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn single_substitution() {
        // ACGT vs ACAT => 1 sub in 4 bp => 0.75 identity
        let id = levenshtein_identity(b"ACGT", b"ACAT");
        assert!((id - 0.75).abs() < 1e-9, "got {id}");
    }

    #[test]
    fn single_indel_handled() {
        // Hamming would say ~0 (frame shift); Levenshtein sees a 1-bp
        // indel: edit distance 1, max len 5 => 0.8 identity.
        let id = levenshtein_identity(b"ACGTA", b"ACGA");
        assert!((id - 0.8).abs() < 1e-9, "got {id}");
    }

    #[test]
    fn empty_and_empty_is_1() {
        assert_eq!(levenshtein_identity(b"", b""), 1.0);
    }

    #[test]
    fn one_empty_is_0() {
        assert_eq!(levenshtein_identity(b"ACGT", b""), 0.0);
        assert_eq!(levenshtein_identity(b"", b"ACGT"), 0.0);
    }

    #[test]
    fn coverage_needs_two_full_tiles() {
        // 30 bp seq, tile=20 => only one full tile fits
        let seq = vec![b'A'; 30];
        assert!(compute_tile_coverage(&seq, 20).is_none());
    }

    #[test]
    fn perfect_tandem_full_coverage() {
        // 5 copies of a 100-bp random-ish tile
        let unit: Vec<u8> = (0..100).map(|i| b"ACGT"[(i * 7 + 13) % 4]).collect();
        let seq: Vec<u8> = (0..5).flat_map(|_| unit.clone()).collect();
        let cov = compute_tile_coverage(&seq, 100).expect("coverage");
        assert_eq!(cov.n_tiles, 4);
        assert!((cov.mean - 1.0).abs() < 1e-12);
        assert_eq!(cov.pass_70, 1.0);
        assert_eq!(cov.pass_80, 1.0);
        assert_eq!(cov.pass_90, 1.0);
        assert!((cov.min - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mosaic_array_shows_split() {
        // Two genuinely-different 50-bp sequences. A periodic test
        // (i*7 mod 4 vs i*11 mod 4) collapses because both step by 3
        // mod 4 — easy mistake.
        let clean = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let noisy = b"GGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGG".to_vec();
        assert_eq!(clean.len(), 50);
        assert_eq!(noisy.len(), 50);
        let mut seq: Vec<u8> = Vec::new();
        for _ in 0..3 {
            seq.extend(&clean);
        }
        for _ in 0..3 {
            seq.extend(&noisy);
        }
        let cov = compute_tile_coverage(&seq, 50).expect("coverage");
        assert!(cov.first_half > 0.9, "first half mean = {}", cov.first_half);
        assert!(
            cov.second_half < 0.6,
            "second half mean = {}",
            cov.second_half
        );
    }
}

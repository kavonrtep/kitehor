//! Internal analysis blocks for same-width mixed detection
//! (`docs/new/detect_m7_plan.md` Q1–Q4 + Q8).
//!
//! M7.1 lays down the load-bearing pieces:
//!
//! - [`build_blocks`] — chops `n_rows` into deterministic blocks.
//!   Adaptive size: `max(min_segment_rows, ceil(n_rows / max_segments_per_array))`.
//!   Phase-shift boundaries are folded in as extra splits.
//!   For HOR/irregular_HOR arrays with a known HOR-unit width,
//!   block boundaries snap to full HOR units so partial units at
//!   block edges can't drive the consensus-identity test.
//! - [`block_consensuses`] — per-block majority-vote consensus at a
//!   caller-chosen comparison width (caller picks `hor_length_bp`
//!   if ≥ `min_complete_units_per_block` complete units fit per
//!   block, else `base_width_bp`).
//! - [`pairwise_identity`] — all-pairs Hamming-with-N-skip identity,
//!   with `min_identity_coverage` filtering. Returns the valid
//!   `IdentityPair` set.
//! - [`pick_medoid`] — chooses the reference block (highest sum of
//!   pairwise identities; ties broken by smallest index).
//!
//! **`AnalysisBlock` is an internal type and never written to
//! `segments.tsv`.** Reported segments (`detect::types::Segment`)
//! are a separate, biology-meaningful concept emitted only for
//! phase-shift boundaries and mixed-class sub-blocks (M7.2).
//!
//! Class behaviour is unchanged at M7.1 — these helpers compute
//! the data; the override that uses it ships in M7.2.

use crate::detect::config::DetectorConfig;
use crate::detect::consensus::consensus_on_slice;

/// One internal analysis block. Row coordinates are at the
/// comparison width's wrap (the same width used for the per-block
/// consensus); `bp` coordinates are derived from those rows times
/// the comparison width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnalysisBlock {
    pub start_row: usize,
    pub end_row: usize,
}

impl AnalysisBlock {
    pub fn n_rows(&self) -> usize {
        self.end_row.saturating_sub(self.start_row)
    }
}

/// Per-pair identity outcome retained for use by the M7.2 override.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdentityPair {
    pub i: usize,
    pub j: usize,
    /// Hamming identity over non-N positions: `matches / (matches + mismatches)`.
    /// Range [0, 1]. Higher = more similar.
    pub identity: f64,
    /// Fraction of consensus length that contributed (non-N in
    /// both consensuses).
    pub coverage: f64,
}

/// Build the analysis blocks for one array.
///
/// `n_rows` is the row count at the comparison width. `unit_rows`
/// is the HOR-unit length in rows at the comparison width (i.e.,
/// `hor_length_bp / comparison_width`) — pass `None` for simple_TR
/// or for HOR arrays where unit-level alignment isn't meaningful.
///
/// `extra_split_rows` lets the caller force splits at phase-shift
/// boundaries; the builder folds these into the deterministic
/// adaptive grid.
///
/// Blocks below `cfg.min_segment_rows` are merged into a
/// neighbour. With unit-aligned mode (`unit_rows = Some(u)` where
/// `u >= 1`), block sizes are coerced to multiples of `u` and the
/// trailing partial unit (if any) is dropped from the final block.
pub fn build_blocks(
    n_rows: usize,
    unit_rows: Option<usize>,
    extra_split_rows: &[usize],
    cfg: &DetectorConfig,
) -> Vec<AnalysisBlock> {
    if n_rows < cfg.min_segment_rows.max(2) {
        // Single block — too few rows to subdivide meaningfully.
        return vec![AnalysisBlock { start_row: 0, end_row: n_rows }];
    }

    // Target block size (rows), adaptive to keep ≤ max_segments_per_array.
    let target_rows = {
        let div_up = (n_rows + cfg.max_segments_per_array - 1) / cfg.max_segments_per_array;
        cfg.min_segment_rows.max(div_up)
    };
    let target_rows = match unit_rows {
        Some(u) if u >= 1 => {
            // Snap target up to a multiple of u — at least one full unit.
            let multiplier = (target_rows + u - 1) / u;
            (multiplier.max(1)) * u
        }
        _ => target_rows,
    };

    // Build deterministic boundary set (in rows). Always include 0
    // and the cap (n_rows or its unit-aligned truncation).
    let cap = match unit_rows {
        Some(u) if u >= 1 => (n_rows / u) * u,
        _ => n_rows,
    };
    if cap < cfg.min_segment_rows {
        return vec![AnalysisBlock { start_row: 0, end_row: n_rows }];
    }

    let mut bounds: Vec<usize> = Vec::new();
    let mut r = 0usize;
    while r + target_rows <= cap {
        bounds.push(r);
        r += target_rows;
    }
    bounds.push(cap);

    // Fold in phase-shift breakpoints (also snap to unit_rows if set).
    for &split in extra_split_rows {
        let aligned = match unit_rows {
            Some(u) if u >= 1 => (split / u) * u,
            _ => split,
        };
        if aligned > 0 && aligned < cap && !bounds.contains(&aligned) {
            bounds.push(aligned);
        }
    }
    bounds.sort();
    bounds.dedup();

    // Materialise blocks; merge anything below min_segment_rows
    // into its successor (or predecessor if it's the last block).
    let mut blocks: Vec<AnalysisBlock> = bounds
        .windows(2)
        .map(|w| AnalysisBlock { start_row: w[0], end_row: w[1] })
        .collect();
    let min_rows = cfg.min_segment_rows;
    let mut i = 0;
    while i < blocks.len() {
        if blocks[i].n_rows() < min_rows && blocks.len() > 1 {
            if i + 1 < blocks.len() {
                // Merge into successor.
                blocks[i + 1].start_row = blocks[i].start_row;
                blocks.remove(i);
                // Don't advance — re-check the merged block.
            } else {
                // Last block too small — merge into predecessor.
                blocks[i - 1].end_row = blocks[i].end_row;
                blocks.remove(i);
                break;
            }
        } else {
            i += 1;
        }
    }
    if blocks.is_empty() {
        blocks.push(AnalysisBlock { start_row: 0, end_row: cap });
    }
    blocks
}

/// Compute per-block consensus at the supplied comparison width.
/// `block.start_row` and `block.end_row` are interpreted at
/// `comparison_width`'s wrap (rows of length `comparison_width`).
///
/// Returns one `Option<Vec<u8>>` per input block — `None` when the
/// block's `[start_row, end_row)` is too narrow to call
/// `consensus_on_slice` (fewer than 2 rows in range).
pub fn block_consensuses(
    seq: &[u8],
    comparison_width: usize,
    blocks: &[AnalysisBlock],
) -> Vec<Option<Vec<u8>>> {
    blocks
        .iter()
        .map(|b| consensus_on_slice(seq, comparison_width, b.start_row, b.end_row))
        .collect()
}

/// Hamming identity ignoring positions where either consensus has
/// `N`. Returns `None` when the consensuses are different lengths
/// or when no informative positions remain.
///
/// The returned tuple is `(identity, coverage)`:
///   - `identity` = matches / (matches + mismatches) over non-N pairs
///   - `coverage` = (matches + mismatches) / consensus_length
pub fn hamming_identity_n_skip(a: &[u8], b: &[u8]) -> Option<(f64, f64)> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    for (x, y) in a.iter().zip(b.iter()) {
        if *x == b'N' || *y == b'N' {
            continue;
        }
        if x == y {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }
    let informative = matches + mismatches;
    if informative == 0 {
        return None;
    }
    let identity = matches as f64 / informative as f64;
    let coverage = informative as f64 / a.len() as f64;
    Some((identity, coverage))
}

/// All-pairs identity over the supplied consensuses, with
/// coverage filtering. Pairs whose coverage is below
/// `cfg.min_identity_coverage` are dropped — they aren't
/// informative enough to drive the mixed override.
///
/// `consensuses` is expected in the same order as
/// `AnalysisBlock` instances; the `i` / `j` fields of returned
/// pairs index into that slice.
pub fn pairwise_identity(
    consensuses: &[Option<Vec<u8>>],
    cfg: &DetectorConfig,
) -> Vec<IdentityPair> {
    let mut pairs: Vec<IdentityPair> = Vec::new();
    for i in 0..consensuses.len() {
        let Some(ci) = &consensuses[i] else { continue };
        for j in (i + 1)..consensuses.len() {
            let Some(cj) = &consensuses[j] else { continue };
            let Some((identity, coverage)) = hamming_identity_n_skip(ci, cj) else {
                continue;
            };
            if coverage < cfg.min_identity_coverage {
                continue;
            }
            pairs.push(IdentityPair { i, j, identity, coverage });
        }
    }
    pairs
}

/// Pick the reference (medoid) block: the index with the highest
/// sum of pairwise identities to the other admitted blocks.
/// Ties broken by smallest index.
///
/// `n_blocks` is the consensuses slice length so we can iterate
/// over indices even when some pairs were filtered out. Returns
/// `None` when no pair survived the coverage filter.
pub fn pick_medoid(n_blocks: usize, pairs: &[IdentityPair]) -> Option<usize> {
    if pairs.is_empty() {
        return None;
    }
    let mut sums = vec![0.0f64; n_blocks];
    let mut counts = vec![0usize; n_blocks];
    for p in pairs {
        sums[p.i] += p.identity;
        sums[p.j] += p.identity;
        counts[p.i] += 1;
        counts[p.j] += 1;
    }
    // Pick max (sum, smallest index) over blocks that participated
    // in at least one valid pair.
    let mut best: Option<(usize, f64)> = None;
    for (idx, (s, &c)) in sums.iter().zip(counts.iter()).enumerate() {
        if c == 0 {
            continue;
        }
        match best {
            None => best = Some((idx, *s)),
            Some((_, prev_s)) if *s > prev_s => best = Some((idx, *s)),
            // ties: keep the smaller index (already true since we iterate ascending)
            _ => {}
        }
    }
    best.map(|(i, _)| i)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_cfg() -> DetectorConfig {
        DetectorConfig::default()
    }

    #[test]
    fn build_blocks_single_for_short_array() {
        // Below min_segment_rows → no subdivision.
        let cfg = default_cfg();
        let blocks = build_blocks(10, None, &[], &cfg);
        assert_eq!(blocks, vec![AnalysisBlock { start_row: 0, end_row: 10 }]);
    }

    #[test]
    fn build_blocks_respects_max_segments_per_array_on_huge() {
        // 10 000 rows, max_segments=32, min_segment_rows=20 →
        // target ≈ ceil(10000/32) = 313 → 32 blocks of 313 rows
        // (last block absorbs the remainder).
        let mut cfg = default_cfg();
        cfg.max_segments_per_array = 32;
        cfg.min_segment_rows = 20;
        let blocks = build_blocks(10_000, None, &[], &cfg);
        assert!(
            blocks.len() <= cfg.max_segments_per_array,
            "{} blocks exceeds cap {}",
            blocks.len(),
            cfg.max_segments_per_array,
        );
        // Each block ≥ min_segment_rows.
        for b in &blocks {
            assert!(b.n_rows() >= cfg.min_segment_rows, "block too small: {:?}", b);
        }
        // Spans the full array.
        assert_eq!(blocks.first().unwrap().start_row, 0);
        assert_eq!(blocks.last().unwrap().end_row, 10_000);
    }

    #[test]
    fn build_blocks_hor_unit_alignment_drops_partial_units_at_end() {
        // 105 rows, unit_rows=12 → cap = (105/12)*12 = 96; target
        // gets snapped up to a multiple of 12.
        let mut cfg = default_cfg();
        cfg.max_segments_per_array = 4;
        cfg.min_segment_rows = 12;
        let blocks = build_blocks(105, Some(12), &[], &cfg);
        for b in &blocks {
            assert_eq!(b.n_rows() % 12, 0, "block not unit-aligned: {:?}", b);
        }
        // Cap drops trailing partial unit (rows 96..105 excluded).
        assert_eq!(blocks.last().unwrap().end_row, 96);
    }

    #[test]
    fn build_blocks_folds_in_extra_splits() {
        // 200 rows, target ≈ 100 → 2 blocks; phase shift at row 50
        // should split the first block into [0,50) + [50,100).
        let mut cfg = default_cfg();
        cfg.max_segments_per_array = 2;
        cfg.min_segment_rows = 20;
        let blocks = build_blocks(200, None, &[50], &cfg);
        let starts: Vec<usize> = blocks.iter().map(|b| b.start_row).collect();
        assert!(starts.contains(&50), "phase-shift split at 50 missing: {starts:?}");
    }

    #[test]
    fn build_blocks_merges_undersized_block_into_neighbour() {
        // 200 rows, target=100; phase shift at row 195 would create
        // a tiny trailing block (5 rows) below min_segment_rows=20.
        // The trailing block should merge into its predecessor.
        let mut cfg = default_cfg();
        cfg.max_segments_per_array = 2;
        cfg.min_segment_rows = 20;
        let blocks = build_blocks(200, None, &[195], &cfg);
        // None of the surviving blocks should be < min_segment_rows.
        for b in &blocks {
            assert!(b.n_rows() >= cfg.min_segment_rows, "tiny block survived: {:?}", b);
        }
    }

    #[test]
    fn hamming_identity_n_skip_basic() {
        let a = b"ACGTACGT";
        let b = b"ACGTACGT";
        let (id, cov) = hamming_identity_n_skip(a, b).unwrap();
        assert!((id - 1.0).abs() < 1e-9);
        assert!((cov - 1.0).abs() < 1e-9);
    }

    #[test]
    fn hamming_identity_n_skip_drops_n_positions() {
        // 8 positions: 4 match, 2 mismatch, 2 N → identity = 4/6.
        let a = b"AAAANNAA";
        let b = b"AAGCNNAT";
        // positions: A-A (m), A-A (m), A-G (mm), A-C (mm), N-N (skip),
        // N-N (skip), A-A (m), A-T (mm) → 3 matches, 3 mismatches over
        // 6 informative.
        let (id, cov) = hamming_identity_n_skip(a, b).unwrap();
        assert!((id - 0.5).abs() < 1e-9, "id={id}");
        assert!((cov - 0.75).abs() < 1e-9, "cov={cov}");
    }

    #[test]
    fn hamming_identity_n_skip_returns_none_for_all_n() {
        let a = b"NNNN";
        let b = b"NNNN";
        assert!(hamming_identity_n_skip(a, b).is_none());
    }

    #[test]
    fn hamming_identity_n_skip_returns_none_for_length_mismatch() {
        let a = b"AAA";
        let b = b"AAAA";
        assert!(hamming_identity_n_skip(a, b).is_none());
    }

    #[test]
    fn pairwise_identity_skips_low_coverage() {
        // Build three "consensuses". Two have N-rate 90% so their
        // mutual coverage is ~10% — below default 0.70.
        let mut cfg = default_cfg();
        cfg.min_identity_coverage = 0.70;
        let n_heavy = vec![b'N'; 9].into_iter().chain([b'A']).collect::<Vec<u8>>();
        let n_heavy2 = vec![b'N'; 9].into_iter().chain([b'A']).collect::<Vec<u8>>();
        let clean = b"ACGTACGTAC".to_vec();
        let cs = vec![Some(n_heavy), Some(n_heavy2), Some(clean)];
        let pairs = pairwise_identity(&cs, &cfg);
        // n_heavy vs n_heavy2: coverage 1/10 = 0.1 → dropped.
        // n_heavy vs clean: coverage 1/10 → dropped.
        // (only fully-clean vs clean would survive, but we only have one clean)
        assert!(pairs.is_empty(), "expected all pairs filtered; got {pairs:?}");
    }

    #[test]
    fn pairwise_identity_admits_high_coverage_pairs() {
        let mut cfg = default_cfg();
        cfg.min_identity_coverage = 0.70;
        let a = b"ACGTACGTAC".to_vec();
        let b = b"ACGTACGTAG".to_vec();
        let c = b"TGCATGCATG".to_vec();
        let cs = vec![Some(a), Some(b), Some(c)];
        let pairs = pairwise_identity(&cs, &cfg);
        assert_eq!(pairs.len(), 3);
        // a vs b: 9/10 matches → identity 0.9.
        let ab = pairs.iter().find(|p| (p.i, p.j) == (0, 1)).unwrap();
        assert!((ab.identity - 0.9).abs() < 1e-9);
    }

    #[test]
    fn pick_medoid_returns_central_block() {
        // 3 blocks. 0 and 1 are similar (identity 0.95);
        // 2 is divergent from both (identity 0.50, 0.55).
        // Sums: blk 0 = 0.95 + 0.50 = 1.45;
        //       blk 1 = 0.95 + 0.55 = 1.50;
        //       blk 2 = 0.50 + 0.55 = 1.05.
        // Medoid = blk 1.
        let pairs = vec![
            IdentityPair { i: 0, j: 1, identity: 0.95, coverage: 1.0 },
            IdentityPair { i: 0, j: 2, identity: 0.50, coverage: 1.0 },
            IdentityPair { i: 1, j: 2, identity: 0.55, coverage: 1.0 },
        ];
        assert_eq!(pick_medoid(3, &pairs), Some(1));
    }

    #[test]
    fn pick_medoid_ties_break_by_smallest_index() {
        // 3 blocks, all pairwise identical 0.9 — all sums equal.
        let pairs = vec![
            IdentityPair { i: 0, j: 1, identity: 0.9, coverage: 1.0 },
            IdentityPair { i: 0, j: 2, identity: 0.9, coverage: 1.0 },
            IdentityPair { i: 1, j: 2, identity: 0.9, coverage: 1.0 },
        ];
        assert_eq!(pick_medoid(3, &pairs), Some(0));
    }

    #[test]
    fn pick_medoid_none_when_no_pairs() {
        assert_eq!(pick_medoid(3, &[]), None);
    }
}

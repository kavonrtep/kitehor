//! Reported segmentation (`detect_impl_plan.md §6.9`, A2; M7 plan
//! §3.1 + §A19).
//!
//! Two distinct sources of `Segment` rows in `segments.tsv`:
//!
//! 1. **Phase-shift segments** — when `n_phase_shifts > 0`, the
//!    array is split at the shift bp positions. Each segment
//!    inherits the whole-array class (light Q6 option). The new
//!    M7.2 columns (`consensus_identity_to_reference` /
//!    `_coverage`) stay `None` here.
//! 2. **Mixed analysis-block segments** (M7.2) — when the array's
//!    final class is `Mixed` after the consensus-identity override,
//!    emit one row per analysis block with the identity columns
//!    filled. The whole-array `phase_shift_positions` are folded
//!    into the block boundaries upstream (see
//!    `analysis_blocks::build_blocks`).
//!
//! These sources are mutually exclusive: when class=Mixed the
//! per-block segments take precedence. Clean single-family arrays
//! with no phase shifts emit no segment rows (`n_segments = 1`).

use crate::detect::analysis_blocks::{AnalysisBlock, IdentityPair};
use crate::detect::types::{Properties, Segment};

/// Per-block context used to emit mixed segments. Bundled so the
/// caller in `mod.rs::run_array_m4` can hand the whole structure
/// to `split()` without recomputing anything.
#[derive(Debug, Clone)]
pub struct MixedBlocksContext {
    pub blocks: Vec<AnalysisBlock>,
    /// Per-block consensus byte string, aligned 1:1 with `blocks`.
    /// `None` when a block was too short to compute consensus.
    pub consensuses: Vec<Option<Vec<u8>>>,
    /// All-pairs identities surviving the coverage gate. Order is
    /// `(i, j)` with `i < j`.
    pub pairs: Vec<IdentityPair>,
    /// Reference (medoid) block index into `blocks`/`consensuses`.
    pub reference_block: usize,
    /// The width at which block consensuses were computed.
    pub comparison_width: usize,
}

/// Emit reported segment rows from the whole-array `Properties`
/// and (optionally) a mixed-blocks context.
///
/// When `mixed` is supplied AND `props.class == Mixed`, emit one
/// row per analysis block with `consensus_identity_to_reference`
/// and `consensus_identity_coverage` filled. Otherwise fall back
/// to the phase-shift split.
pub fn split(props: &Properties, mixed: Option<&MixedBlocksContext>) -> Vec<Segment> {
    if let Some(ctx) = mixed {
        if matches!(props.class, crate::detect::types::Class::Mixed) {
            return split_mixed(props, ctx);
        }
    }
    split_by_phase_shifts(props)
}

fn split_by_phase_shifts(props: &Properties) -> Vec<Segment> {
    if props.n_phase_shifts == 0 || props.phase_shift_positions.is_empty() {
        return Vec::new();
    }
    let mut boundaries: Vec<usize> = vec![0];
    boundaries.extend_from_slice(&props.phase_shift_positions);
    boundaries.push(props.length_bp);
    boundaries.sort();
    boundaries.dedup();

    let mut out = Vec::with_capacity(boundaries.len().saturating_sub(1));
    for (i, win) in boundaries.windows(2).enumerate() {
        let start_bp = win[0];
        let end_bp = win[1];
        if end_bp <= start_bp {
            continue;
        }
        out.push(Segment {
            array_id: props.array_id.clone(),
            segment_id: i + 1,
            start_bp,
            end_bp,
            class: props.class,
            base_width_bp: props.base_width_bp,
            hor_k: props.hor_k,
            column_conservation: props.column_conservation,
            phase_separation: props.phase_separation,
            wobble_amplitude_bp: props.wobble_amplitude_bp,
            irregularity_score: props.irregularity_score,
            consensus_identity_to_reference: None,
            consensus_identity_coverage: None,
        });
    }
    out
}

fn split_mixed(props: &Properties, ctx: &MixedBlocksContext) -> Vec<Segment> {
    let cw = ctx.comparison_width;
    let mut out = Vec::with_capacity(ctx.blocks.len());
    for (i, blk) in ctx.blocks.iter().enumerate() {
        let (identity, coverage) = identity_to_reference(i, ctx);
        out.push(Segment {
            array_id: props.array_id.clone(),
            segment_id: i + 1,
            start_bp: blk.start_row * cw,
            end_bp: blk.end_row * cw,
            class: props.class,
            // For mixed arrays the whole-array base_width/k/etc. are
            // None (per M7 plan Q5 + review #1); segments inherit that.
            base_width_bp: None,
            hor_k: None,
            column_conservation: None,
            phase_separation: None,
            wobble_amplitude_bp: None,
            irregularity_score: None,
            consensus_identity_to_reference: identity,
            consensus_identity_coverage: coverage,
        });
    }
    out
}

fn identity_to_reference(
    block_idx: usize,
    ctx: &MixedBlocksContext,
) -> (Option<f64>, Option<f64>) {
    if block_idx == ctx.reference_block {
        return (Some(1.0), Some(1.0));
    }
    for p in &ctx.pairs {
        let other = if p.i == block_idx {
            p.j
        } else if p.j == block_idx {
            p.i
        } else {
            continue;
        };
        if other == ctx.reference_block {
            return (Some(p.identity), Some(p.coverage));
        }
    }
    (None, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::analysis_blocks::IdentityPair;
    use crate::detect::types::Class;

    fn props_with_shifts(length: usize, shifts: Vec<usize>) -> Properties {
        let mut p = Properties::placeholder("arr", length);
        p.class = Class::HOR;
        p.base_width_bp = Some(171);
        p.hor_k = Some(12);
        p.phase_shift_positions = shifts.clone();
        p.n_phase_shifts = shifts.len();
        p.n_segments = 1 + shifts.len();
        p
    }

    #[test]
    fn no_shifts_yields_no_segments() {
        let p = Properties::placeholder("arr", 1000);
        assert!(split(&p, None).is_empty());
    }

    #[test]
    fn one_shift_yields_two_segments() {
        let p = props_with_shifts(1000, vec![400]);
        let s = split(&p, None);
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].start_bp, 0);
        assert_eq!(s[0].end_bp, 400);
        assert_eq!(s[1].start_bp, 400);
        assert_eq!(s[1].end_bp, 1000);
        for seg in &s {
            assert_eq!(seg.class, Class::HOR);
            assert_eq!(seg.base_width_bp, Some(171));
            assert_eq!(seg.hor_k, Some(12));
            // Phase-shift segments leave the M7.2 columns blank.
            assert_eq!(seg.consensus_identity_to_reference, None);
            assert_eq!(seg.consensus_identity_coverage, None);
        }
    }

    #[test]
    fn two_shifts_yields_three_segments() {
        let p = props_with_shifts(1000, vec![300, 700]);
        let s = split(&p, None);
        assert_eq!(s.len(), 3);
        assert_eq!(s[0].end_bp, 300);
        assert_eq!(s[1].start_bp, 300);
        assert_eq!(s[1].end_bp, 700);
        assert_eq!(s[2].start_bp, 700);
        assert_eq!(s[2].end_bp, 1000);
    }

    #[test]
    fn out_of_order_shifts_get_sorted_and_deduped() {
        let p = props_with_shifts(1000, vec![700, 300, 700]);
        let s = split(&p, None);
        assert_eq!(s.len(), 3);
    }

    // M7.2: mixed-class arrays emit one segment per analysis block.
    #[test]
    fn mixed_emits_one_segment_per_block_with_identity() {
        let mut p = Properties::placeholder("arr", 4000);
        p.class = Class::Mixed;
        let ctx = MixedBlocksContext {
            blocks: vec![
                AnalysisBlock { start_row: 0,  end_row: 10 },
                AnalysisBlock { start_row: 10, end_row: 20 },
                AnalysisBlock { start_row: 20, end_row: 40 },
            ],
            consensuses: vec![Some(b"A".to_vec()), Some(b"A".to_vec()), Some(b"T".to_vec())],
            // Block 0 ↔ 1: high identity; block 0/1 ↔ 2: low identity.
            pairs: vec![
                IdentityPair { i: 0, j: 1, identity: 0.95, coverage: 1.0 },
                IdentityPair { i: 0, j: 2, identity: 0.50, coverage: 1.0 },
                IdentityPair { i: 1, j: 2, identity: 0.55, coverage: 1.0 },
            ],
            reference_block: 0,
            comparison_width: 100,
        };
        let s = split(&p, Some(&ctx));
        assert_eq!(s.len(), 3);
        // Reference block (index 0) → identity to itself is 1.0.
        assert_eq!(s[0].consensus_identity_to_reference, Some(1.0));
        // Block 1: identity to reference (block 0) = 0.95.
        assert!((s[1].consensus_identity_to_reference.unwrap() - 0.95).abs() < 1e-9);
        // Block 2: identity to reference = 0.50.
        assert!((s[2].consensus_identity_to_reference.unwrap() - 0.50).abs() < 1e-9);
        // bp coordinates derived from start_row × comparison_width.
        assert_eq!(s[0].start_bp, 0);
        assert_eq!(s[0].end_bp, 1000);
        assert_eq!(s[1].start_bp, 1000);
        assert_eq!(s[2].end_bp, 4000);
        // Mixed class on every row.
        for seg in &s {
            assert_eq!(seg.class, Class::Mixed);
            // No base_width/k/etc. for mixed segments (per A19 light Q6 + review #1).
            assert_eq!(seg.base_width_bp, None);
            assert_eq!(seg.hor_k, None);
        }
    }

    #[test]
    fn mixed_context_ignored_when_class_not_mixed() {
        // If class is HOR, the mixed context is ignored — phase-shift
        // path still wins. Verifies the mutex between the two paths.
        let p = props_with_shifts(1000, vec![400]);
        let ctx = MixedBlocksContext {
            blocks: vec![AnalysisBlock { start_row: 0, end_row: 10 }],
            consensuses: vec![Some(b"A".to_vec())],
            pairs: vec![],
            reference_block: 0,
            comparison_width: 100,
        };
        let s = split(&p, Some(&ctx));
        assert_eq!(s.len(), 2, "HOR path should still emit phase-shift segments");
        assert_eq!(s[0].class, Class::HOR);
    }
}

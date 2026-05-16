//! Phase-shift segmentation (`detect_impl_plan.md §6.9`, A2).
//!
//! When `n_phase_shifts > 0`, the array is split at the shift bp
//! positions and each segment is reported in `segments.tsv` with its
//! own class / base_width / k. For MVP M4, segments inherit the
//! whole-array class because all CI fixtures with phase shifts share
//! the same underlying repeat architecture on both sides of every
//! shift (`detect_impl_plan.md §A5`). A future iteration can
//! recompute per-segment widths once HMM-based segmentation lands
//! (OQ6).

use crate::detect::types::{Class, Properties, Segment};

/// Emit one segment per inter-shift region. `phase_shift_positions`
/// are sorted bp positions inside the array.
pub fn split(props: &Properties) -> Vec<Segment> {
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
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
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
        assert!(split(&p).is_empty());
    }

    #[test]
    fn one_shift_yields_two_segments() {
        let p = props_with_shifts(1000, vec![400]);
        let s = split(&p);
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].start_bp, 0);
        assert_eq!(s[0].end_bp, 400);
        assert_eq!(s[1].start_bp, 400);
        assert_eq!(s[1].end_bp, 1000);
        for seg in &s {
            assert_eq!(seg.class, Class::HOR);
            assert_eq!(seg.base_width_bp, Some(171));
            assert_eq!(seg.hor_k, Some(12));
        }
    }

    #[test]
    fn two_shifts_yields_three_segments() {
        let p = props_with_shifts(1000, vec![300, 700]);
        let s = split(&p);
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
        let s = split(&p);
        assert_eq!(s.len(), 3);
    }
}

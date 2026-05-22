//! Sigmoid logit confidence score (`detect_impl_plan.md §6.11`,
//! spec §9).
//!
//! ```text
//! HOR / irregular_HOR:
//!     logit = α · phase_separation
//!           + γ · log10(n_complete_copies + 1)
//!           + δ · mean_column_IC
//!           − ε · irregularity_score
//!           − ζ · wobble_amplitude / width
//!           − η · |mean_shift| / width
//!
//! simple_TR:
//!     logit = α · column_IC                 (using IC in place of phase_sep)
//!           + γ · log10(n_complete_copies + 1)
//!           + δ · column_IC
//!           − ε · irregularity_score
//!           − ζ · wobble_amplitude / width
//!           − η · |mean_shift| / width
//!
//! mixed / ambiguous (Review-2026-05-16 #6): logit reflects the
//! quality of the negative evidence rather than a class constant.
//!   `mixed`: weak positive bias when multiple high-score periods
//!            are present (we're more sure the array is composite
//!            when several independent periods voted). Bias scales
//!            with phase_separation if it was computed before the
//!            classifier gave up.
//!   `ambiguous`: small negative bias; tightens further when no
//!            width supported the rescue floor (truly no signal).
//! ```
//!
//! **This is a heuristic score, not a calibrated probability.**
//! Weights live in `DetectorConfig.confidence_weights` so M6
//! calibration can adjust them without touching code. Default values
//! produce ≥ 0.85 for clean HOR / simple_TR cases on the CI corpus
//! and ≤ 0.5 for the negative-control fixture. Downstream consumers
//! should treat this as a per-class signal-quality estimate, not
//! P(class is correct).

use crate::detect::config::DetectorConfig;
use crate::detect::types::{Class, Properties};

pub fn compute(props: &Properties, cfg: &DetectorConfig) -> f64 {
    let w = &cfg.confidence_weights;
    let width = props.base_width_bp.unwrap_or(1).max(1) as f64;
    let n_copies = props.n_complete_copies.unwrap_or(0) as f64;
    let copies_term = w.gamma * (n_copies + 1.0).log10();
    let ic = props.column_conservation.unwrap_or(0.0);
    let irreg = props.irregularity_score.unwrap_or(0.0);
    let wob = props.wobble_amplitude_bp.unwrap_or(0.0).abs() / width;
    let shift = props.mean_shift_bp.unwrap_or(0.0).abs() / width;

    let logit = match props.class {
        Class::HOR | Class::IrregularHOR => {
            let phase = props.phase_separation.unwrap_or(0.0);
            w.alpha * phase + copies_term + w.delta * ic
                - w.epsilon * irreg
                - w.zeta * wob
                - w.eta * shift
        }
        Class::SimpleTR => {
            // For simple_TR, `phase_separation` is small by
            // construction. Substitute column_IC for both the
            // "phase" and "IC" terms so high IC + high copy count
            // still produces high confidence.
            w.alpha * ic + copies_term + w.delta * ic
                - w.epsilon * irreg
                - w.zeta * wob
                - w.eta * shift
        }
        Class::Mixed => {
            // Review-2026-05-16 #6: derive from evidence rather
            // than the previous class-constant -1.0. After the
            // finding-#1 fix, mixed retains `n_complete_copies`
            // (the structural-evidence count from the decision)
            // and `n_phase_shifts` (an independent indicator of
            // architectural inconsistency). Use both: more copies
            // and more shifts make us more confident the array
            // really is structurally composite.
            let shifts_term = (props.n_phase_shifts as f64).ln_1p();
            -0.8 + copies_term + 0.3 * shifts_term
        }
        Class::Ambiguous => {
            // No recovered signal by definition — keep the bias
            // negative, but let it vary slightly with copy count
            // so two ambiguous calls don't collide on a single
            // class-constant confidence value.
            -1.5 + 0.3 * copies_term
        }
    };
    sigmoid(logit)
}

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::types::Properties;

    fn make(class: Class, ic: f64, phase_sep: f64, n_copies: usize) -> Properties {
        let mut p = Properties::placeholder("x", 1_000_000);
        p.class = class;
        p.column_conservation = Some(ic);
        p.phase_separation = Some(phase_sep);
        p.n_complete_copies = Some(n_copies);
        p.base_width_bp = Some(171);
        p.wobble_amplitude_bp = Some(0.0);
        p.mean_shift_bp = Some(0.0);
        p.irregularity_score = Some(0.0);
        p
    }

    #[test]
    fn clean_hor_scores_high() {
        let cfg = DetectorConfig::default();
        let p = make(Class::HOR, 1.5, 0.30, 200);
        let c = compute(&p, &cfg);
        assert!(
            c > 0.85,
            "expected confidence > 0.85 for clean HOR; got {c}"
        );
    }

    #[test]
    fn clean_simple_tr_scores_high() {
        let cfg = DetectorConfig::default();
        let p = make(Class::SimpleTR, 1.8, 0.0, 1000);
        let c = compute(&p, &cfg);
        assert!(
            c > 0.85,
            "expected confidence > 0.85 for clean simple TR; got {c}"
        );
    }

    #[test]
    fn ambiguous_scores_low() {
        let cfg = DetectorConfig::default();
        let p = make(Class::Ambiguous, 0.0, 0.0, 0);
        let c = compute(&p, &cfg);
        assert!(c < 0.5, "expected confidence < 0.5 for ambiguous; got {c}");
    }

    // Review-2026-05-16 #6: mixed/ambiguous reflect evidence
    // rather than a class constant. Each call should vary with
    // its recovered properties.
    #[test]
    fn mixed_with_more_copies_scores_higher() {
        let cfg = DetectorConfig::default();
        let mixed_high = make(Class::Mixed, 0.0, 0.0, 500);
        let mixed_low = make(Class::Mixed, 0.0, 0.0, 5);
        assert!(
            compute(&mixed_high, &cfg) > compute(&mixed_low, &cfg),
            "mixed with more complete copies should score higher"
        );
    }

    #[test]
    fn mixed_outscores_signal_free_ambiguous() {
        let cfg = DetectorConfig::default();
        let mixed = make(Class::Mixed, 0.0, 0.0, 200);
        let ambiguous = make(Class::Ambiguous, 0.0, 0.0, 0);
        assert!(
            compute(&mixed, &cfg) > compute(&ambiguous, &cfg),
            "mixed (some structure detected) should outscore signal-free ambiguous"
        );
    }

    #[test]
    fn high_wobble_penalises_score() {
        let cfg = DetectorConfig::default();
        let mut p = make(Class::HOR, 1.5, 0.30, 200);
        let clean = compute(&p, &cfg);
        p.wobble_amplitude_bp = Some(20.0); // huge wobble for w=171
        let dirty = compute(&p, &cfg);
        assert!(dirty < clean, "high wobble should reduce confidence");
    }
}

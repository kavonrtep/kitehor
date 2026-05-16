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
//! mixed / ambiguous: small negative bias → confidence ≈ 0.27.
//! ```
//!
//! Weights live in `DetectorConfig.confidence_weights` so M6
//! calibration can adjust them without touching code. Default values
//! produce ≥ 0.85 for clean HOR / simple_TR cases on the CI corpus
//! and ≤ 0.5 for the negative-control fixture.

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
        Class::Mixed | Class::Ambiguous => -1.0,
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
        assert!(c > 0.85, "expected confidence > 0.85 for clean HOR; got {c}");
    }

    #[test]
    fn clean_simple_tr_scores_high() {
        let cfg = DetectorConfig::default();
        let p = make(Class::SimpleTR, 1.8, 0.0, 1000);
        let c = compute(&p, &cfg);
        assert!(c > 0.85, "expected confidence > 0.85 for clean simple TR; got {c}");
    }

    #[test]
    fn ambiguous_scores_low() {
        let cfg = DetectorConfig::default();
        let p = make(Class::Ambiguous, 0.0, 0.0, 0);
        let c = compute(&p, &cfg);
        assert!(c < 0.5, "expected confidence < 0.5 for ambiguous; got {c}");
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

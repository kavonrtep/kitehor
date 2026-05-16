//! Array-level classification per `detect_spec.md §8` and
//! `detect_impl_plan.md §6.10`.
//!
//! The regime A/B/C decision is the load-bearing piece: HOR is
//! called only when **both** the base period and the HOR-unit
//! period are statistically valid. Regime-C inputs (high inter-slot
//! divergence) fall through to `simple_TR` at the HOR-unit width.

use crate::detect::autocorr;
use crate::detect::config::DetectorConfig;
use crate::detect::embed;
use crate::detect::phase;
use crate::detect::types::{Class, PeriodCandidate, WidthFeatures};
use crate::detect::wrap;
use crate::sequence::ArrayRecord;
use std::collections::HashMap;

/// Output of `decide_array`: a class call + the supporting evidence.
#[derive(Debug, Clone)]
pub struct ArrayDecision {
    pub class: Class,
    pub base_width_bp: Option<usize>,
    pub hor_k: Option<usize>,
    pub hor_length_bp: Option<usize>,
    pub n_complete_copies: Option<usize>,
    pub column_conservation: Option<f64>,
    pub phase_separation: Option<f64>,
    pub inter_monomer_identity: Option<f64>,
    pub reason: String,
}

#[derive(Debug, Clone)]
struct Candidate {
    class: Class,
    base_width_bp: usize,
    hor_k: Option<usize>,
    hor_length_bp: Option<usize>,
    column_conservation: f64,
    phase_separation: f64,
    r_lag1: f64,
    n_complete_copies: usize,
    reason: String,
    /// For regime-C-derived simple_TRs, the original HOR base width
    /// that was collapsing. None for plain simple_TR candidates.
    /// Used to suppress harmonics of the underlying base width.
    underlying_base: Option<usize>,
}

/// Main entry point. Iterates input periods by score, evaluates each
/// width as an HOR candidate / simple_TR / regime-C fall-through,
/// then combines candidates.
pub fn decide_array(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    width_features: &[WidthFeatures],
    cfg: &DetectorConfig,
) -> ArrayDecision {
    let by_width: HashMap<usize, &WidthFeatures> =
        width_features.iter().map(|w| (w.width_bp, w)).collect();

    // Sort input periods by score (desc); tie-break smaller-width-first
    // so we prefer base over harmonic at equal score.
    let mut sorted: Vec<&PeriodCandidate> = pers.iter().collect();
    sorted.sort_by(|a, b| {
        b.period_score
            .partial_cmp(&a.period_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.period_bp.cmp(&b.period_bp))
    });

    if width_features.is_empty() {
        return ambiguous("no widths available");
    }

    // No width clears `ic_threshold_min` → array has no detectable
    // repeat structure.
    let any_supported = width_features
        .iter()
        .any(|w| w.column_ic.map(|ic| ic >= cfg.ic_threshold_min).unwrap_or(false));
    if !any_supported {
        return ambiguous("no width achieves ic_threshold_min");
    }

    let mut hor_calls_raw: Vec<Candidate> = Vec::new();
    let mut simple_calls_raw: Vec<Candidate> = Vec::new();

    let bg = wrap::Background::compute(&arr.seq);
    for p in &sorted {
        let w = match by_width.get(&p.period_bp) {
            Some(w) if w.rows >= cfg.min_rows_per_width => *w,
            _ => continue,
        };
        let ic = w.column_ic.unwrap_or(0.0);

        // Recompute R(k) here so we have the full curve for primitive
        // correction + phase separation. `width_features` only stores
        // the summary.
        let embs = embed::embed_rows(&arr.seq, w.width_bp, cfg);
        let summary = autocorr::compute(&embs, cfg.max_hor_k);
        let r_k = summary.r_k.clone();
        let r1 = summary.r_lag1.unwrap_or(0.0);
        let best_lag_raw = summary.best_lag.unwrap_or(1);
        let k_corrected =
            phase::primitive_correct(&r_k, best_lag_raw, cfg.primitive_correction_delta);
        let phase_sep = phase::phase_separation(&r_k, k_corrected);
        let n_units = w.rows / k_corrected.max(1);

        // ---- HOR candidate (regime B) ----
        if k_corrected >= 2
            && phase_sep >= cfg.phase_separation_threshold
            && n_units >= cfg.min_hor_units
        {
            let unit_width = w.width_bp * k_corrected;
            let unit_ic = by_width
                .get(&unit_width)
                .and_then(|u| u.column_ic)
                .unwrap_or_else(|| recompute_ic(&arr.seq, unit_width, &bg, cfg));
            let base_ic_ok = ic >= cfg.ic_threshold_hor_base;
            let unit_ic_ok = unit_ic >= cfg.ic_threshold_hor_unit;
            let r1_ok = r1 >= cfg.regime_c_r1_threshold;
            if base_ic_ok && unit_ic_ok && r1_ok {
                hor_calls_raw.push(Candidate {
                    class: Class::HOR,
                    base_width_bp: w.width_bp,
                    hor_k: Some(k_corrected),
                    hor_length_bp: Some(unit_width),
                    column_conservation: ic,
                    phase_separation: phase_sep,
                    r_lag1: r1,
                    n_complete_copies: n_units,
                    reason: format!(
                        "regime B — base width {} bp + HOR-unit width {} bp both valid (k={})",
                        w.width_bp, unit_width, k_corrected
                    ),
                    underlying_base: None,
                });
                continue;
            }
            // Regime C: HOR collapsed — base period not statistically
            // valid (low R(1) OR low base IC), but the HOR-unit width
            // is a strong simple_TR. Either trigger qualifies.
            let regime_c =
                (!r1_ok || !base_ic_ok) && unit_ic >= cfg.ic_threshold_simple_tr;
            if regime_c {
                simple_calls_raw.push(Candidate {
                    class: Class::SimpleTR,
                    base_width_bp: unit_width,
                    hor_k: None,
                    hor_length_bp: None,
                    column_conservation: unit_ic,
                    phase_separation: 0.0,
                    r_lag1: r1,
                    n_complete_copies: n_units,
                    reason: format!(
                        "regime C — HOR with k={} collapsed; simple_TR at HOR-unit width {} bp (R(1)={:.3}, base IC={:.3})",
                        k_corrected, unit_width, r1, ic
                    ),
                    underlying_base: Some(w.width_bp),
                });
                continue;
            }
        }

        // ---- Simple TR candidate ----
        if ic >= cfg.ic_threshold_simple_tr && phase_sep < cfg.phase_separation_threshold {
            // Regime-A tag: uniformly high R(k) curve (R(1) ≈ R(best_lag))
            // suggests a collapsed HOR (div ≈ 0). The data is genuinely
            // indistinguishable from a plain simple_TR, so we always
            // call simple_TR — we just flag the "regime A" possibility
            // in the reason for downstream consumers.
            let r_best = summary.best_lag_score.unwrap_or(0.0);
            let is_regime_a = best_lag_raw >= 2
                && r1 >= cfg.regime_a_r1_floor
                && (r_best - r1).abs() < cfg.regime_a_r_curve_flatness;
            let reason = if is_regime_a {
                format!(
                    "regime A — uniformly high R(k); simple_TR at base width {} bp (IC {:.3})",
                    w.width_bp, ic
                )
            } else {
                format!(
                    "simple_TR — base width {} bp; column IC {:.3}, phase_sep {:.3}",
                    w.width_bp, ic, phase_sep
                )
            };
            simple_calls_raw.push(Candidate {
                class: Class::SimpleTR,
                base_width_bp: w.width_bp,
                hor_k: None,
                hor_length_bp: None,
                column_conservation: ic,
                phase_separation: phase_sep,
                r_lag1: r1,
                n_complete_copies: w.rows,
                reason,
                underlying_base: None,
            });
            continue;
        }
    }

    // ---- Dedup candidates ----
    //
    // Two cleanups, both prevent the same array from being called
    // `mixed` just because its repeat structure produces multiple
    // candidate widths along the same harmonic chain:
    //
    //   (1) Within each list, sort by base_width ascending and drop
    //       any candidate whose base_width is a strict multiple of an
    //       earlier candidate's base_width — primitive correction at
    //       the candidate level.
    //   (2) Suppress simple_TR candidates whose base_width matches
    //       the `hor_length_bp` of an existing HOR call. Those are
    //       the unit-band signature of the same HOR, not an
    //       independent claim.
    let hor_calls = dedup_by_multiplicity(hor_calls_raw);
    let mut simple_calls = dedup_by_multiplicity(simple_calls_raw);
    // Suppress any simple_TR candidate whose base_width is a multiple
    // of an existing HOR's base_width (so the HOR-unit width 12·171
    // and the harmonics 2·171, 3·171 are all absorbed into the HOR
    // call rather than producing spurious mixed-class artefacts).
    simple_calls.retain(|c| {
        !hor_calls.iter().any(|h| {
            h.base_width_bp > 0 && c.base_width_bp % h.base_width_bp == 0
        })
    });
    // Regime-C cleanup: a regime-C simple_TR at HOR-unit width
    // suppresses every plain simple_TR whose width is a multiple of
    // the *underlying* HOR base. Otherwise harmonics of the collapsed
    // base (e.g. 3·170 = 510 for T07) fire as independent simple_TR
    // claims and we end up calling `mixed`.
    let regime_c_underlying_bases: std::collections::HashSet<usize> = simple_calls
        .iter()
        .filter_map(|c| c.underlying_base)
        .collect();
    simple_calls.retain(|c| {
        c.underlying_base.is_some()
            || !regime_c_underlying_bases
                .iter()
                .any(|&b| b > 0 && c.base_width_bp % b == 0)
    });

    // ---- Combine candidates into a class ----

    if hor_calls.len() >= 2 {
        let first = &hor_calls[0];
        let any_diff = hor_calls.iter().any(|c| {
            c.base_width_bp != first.base_width_bp || c.hor_k != first.hor_k
        });
        if any_diff {
            return mixed_decision(
                hor_calls.iter().map(|c| c.n_complete_copies).sum::<usize>(),
                "mixed — multiple HOR candidates with distinct (base_width, k)",
            );
        }
    }
    if simple_calls.len() >= 2 && hor_calls.is_empty() {
        let first = &simple_calls[0];
        let any_diff = simple_calls
            .iter()
            .any(|c| c.base_width_bp != first.base_width_bp);
        if any_diff {
            return mixed_decision(
                simple_calls.iter().map(|c| c.n_complete_copies).sum::<usize>(),
                "mixed — multiple simple_TR candidates with distinct base widths",
            );
        }
    }
    if !hor_calls.is_empty() && !simple_calls.is_empty() {
        let hor_base = hor_calls[0].base_width_bp;
        let any_simple_diff = simple_calls
            .iter()
            .any(|c| c.base_width_bp != hor_base);
        if any_simple_diff {
            return mixed_decision(
                None,
                "mixed — HOR and simple_TR candidates at different base widths",
            );
        }
    }

    // Single coherent HOR call wins.
    if let Some(c) = hor_calls.into_iter().next() {
        return ArrayDecision {
            class: c.class,
            base_width_bp: Some(c.base_width_bp),
            hor_k: c.hor_k,
            hor_length_bp: c.hor_length_bp,
            n_complete_copies: Some(c.n_complete_copies),
            column_conservation: Some(c.column_conservation),
            phase_separation: Some(c.phase_separation),
            inter_monomer_identity: Some(c.r_lag1),
            reason: c.reason,
        };
    }
    if let Some(c) = simple_calls.into_iter().next() {
        return ArrayDecision {
            class: c.class,
            base_width_bp: Some(c.base_width_bp),
            hor_k: None,
            hor_length_bp: None,
            n_complete_copies: Some(c.n_complete_copies),
            column_conservation: Some(c.column_conservation),
            phase_separation: Some(c.phase_separation),
            inter_monomer_identity: None,
            reason: c.reason,
        };
    }

    // Coexisting-periods fallback for `mixed`: when no candidate
    // passed individually but the upstream period generator emitted
    // two or more high-score periods that aren't multiples of each
    // other, trust that signal. This handles T13-style cases where
    // each block's HOR-12 signal is diluted by the other block's
    // rows at the whole-array level.
    let mut strong: Vec<usize> = pers
        .iter()
        .filter(|p| p.period_score >= cfg.strong_period_score)
        .map(|p| p.period_bp)
        .collect();
    strong.sort();
    strong.dedup();
    let mut strong_primary: Vec<usize> = Vec::new();
    for w in &strong {
        if !strong_primary.iter().any(|p| *w % p == 0) {
            strong_primary.push(*w);
        }
    }
    if strong_primary.len() >= 2 {
        return ArrayDecision {
            class: Class::Mixed,
            base_width_bp: None,
            hor_k: None,
            hor_length_bp: None,
            n_complete_copies: None,
            column_conservation: None,
            phase_separation: None,
            inter_monomer_identity: None,
            reason: format!(
                "mixed — upstream emitted {} high-score periods after multiplicity dedup ({:?})",
                strong_primary.len(),
                strong_primary
            ),
        };
    }

    ambiguous("no width passed HOR or simple_TR criteria")
}

fn ambiguous(reason: &str) -> ArrayDecision {
    ArrayDecision {
        class: Class::Ambiguous,
        base_width_bp: None,
        hor_k: None,
        hor_length_bp: None,
        n_complete_copies: None,
        column_conservation: None,
        phase_separation: None,
        inter_monomer_identity: None,
        reason: reason.to_string(),
    }
}

fn mixed_decision(n: impl Into<Option<usize>>, reason: &str) -> ArrayDecision {
    ArrayDecision {
        class: Class::Mixed,
        base_width_bp: None,
        hor_k: None,
        hor_length_bp: None,
        n_complete_copies: n.into(),
        column_conservation: None,
        phase_separation: None,
        inter_monomer_identity: None,
        reason: reason.to_string(),
    }
}

fn dedup_by_multiplicity(mut v: Vec<Candidate>) -> Vec<Candidate> {
    v.sort_by_key(|c| c.base_width_bp);
    let mut out: Vec<Candidate> = Vec::new();
    for c in v {
        let is_multiple = out
            .iter()
            .any(|prev| prev.base_width_bp > 0 && c.base_width_bp % prev.base_width_bp == 0);
        if !is_multiple {
            out.push(c);
        }
    }
    out
}

fn recompute_ic(
    seq: &[u8],
    width: usize,
    bg: &wrap::Background,
    cfg: &DetectorConfig,
) -> f64 {
    wrap::wrap_and_ic(seq, width, bg, cfg)
        .map(|s| s.mean_column_ic)
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    // Decision tests live in tests/detect_m4.rs — they need a real
    // synth → detect pipeline run, which is integration-test territory.
}

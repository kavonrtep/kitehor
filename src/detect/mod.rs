//! v2 line-width detector (`kitehor detect*`).
//!
//! Design contract: `docs/new/detect_impl_plan.md`. Built across
//! milestones M0..M6; this file scaffolds the public entry points
//! and the pipeline skeleton. Algorithm modules (`widths`, `wrap`,
//! `embed`, ...) are stubs at M0 and fill in milestone by milestone.

pub mod config;
pub mod io;
pub mod types;

// Algorithm modules — empty stubs at M0; populated M1..M5.
pub mod analysis_blocks;
pub mod autocorr;
pub mod classify;
pub mod confidence;
pub mod consensus;
pub mod edges;
pub mod embed;
pub mod irregularity;
pub mod phase;
pub mod segment;
pub mod shift;
pub mod viz;
pub mod widths;
pub mod wrap;

pub use config::DetectorConfig;
pub use consensus::ConsensusRecord;
pub use viz::VizFlags;
pub use types::{
    Class, ClassHint, PeriodCandidate, Properties, Segment, WidthFeatures,
    PROPERTIES_HEADER, SEGMENTS_HEADER, WIDTH_FEATURES_HEADER,
};

use crate::sequence::ArrayRecord;
use anyhow::Result;
use std::path::Path;

/// End-to-end pipeline for a single FASTA + period-TSV pair.
///
/// M0 behaviour: produces `properties.tsv` with one
/// `Properties::placeholder` row per FASTA record. `segments.tsv` and
/// `width_features.tsv` are header-only. No detection logic yet —
/// that arrives in M1+.
pub fn run_one(
    fasta: &Path,
    periods: &Path,
    out_prefix: &Path,
    cfg: &DetectorConfig,
    viz_flags: &VizFlags,
    allow_missing_periods: bool,
    allow_extra_periods: bool,
) -> Result<DetectorReport> {
    cfg.validate()?;
    let arrays = io::load_arrays(fasta)?;
    let periods_by_id = io::load_periods(periods)?;
    let default_stem = fasta
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string());
    let paired = io::join_arrays_with_periods(
        arrays,
        periods_by_id,
        default_stem.as_deref(),
        allow_missing_periods,
        allow_extra_periods,
    )?;

    let mut properties: Vec<Properties> = Vec::with_capacity(paired.len());
    let mut segments: Vec<Segment> = Vec::new();
    let mut width_features: Vec<WidthFeatures> = Vec::new();
    let mut consensus_records: Vec<ConsensusRecord> = Vec::new();

    for (arr, pers) in &paired {
        let (props, mut widths, mixed_ctx) = run_array_m4(arr, pers, cfg);
        // DH2 + M7.2: emit segments. The mixed context (when present)
        // wins over phase-shift segments because the array class
        // already became Mixed inside run_array_m4.
        let mut new_segments = segment::split(&props, mixed_ctx.as_ref());
        segments.append(&mut new_segments);
        // M5: build consensus + viz when we have a chosen base width.
        // Review-2026-05-16 #1: only emit for resolved classes so
        // Mixed/Ambiguous don't leak a heuristic-width consensus into
        // the output bundle.
        let resolved = matches!(
            props.class,
            Class::SimpleTR | Class::HOR | Class::IrregularHOR
        );
        if resolved {
            if let Some(base_w) = props.base_width_bp {
                if let Some(monomer) = consensus::consensus(&arr.seq, base_w) {
                    let hor_unit = props
                        .hor_length_bp
                        .and_then(|hu| consensus::consensus(&arr.seq, hu));
                    consensus_records.push(ConsensusRecord {
                        array_id: arr.id.clone(),
                        monomer,
                        hor_unit,
                        hor_k: props.hor_k,
                    });
                }
                if viz_flags.is_active() {
                    emit_viz(arr, base_w, cfg, viz_flags)?;
                }
            }
        }
        properties.push(props);
        width_features.append(&mut widths);
    }

    io::write_properties(out_prefix, &properties)?;
    io::write_segments(out_prefix, &segments)?;
    io::write_width_features(out_prefix, &width_features)?;
    io::write_diagnostics(out_prefix, &properties, &width_features, &segments)?;
    if !consensus_records.is_empty() {
        consensus::write_fasta(out_prefix, &consensus_records)?;
    }

    Ok(DetectorReport {
        n_arrays: paired.len(),
        n_segments: segments.len(),
        n_width_rows: width_features.len(),
    })
}

fn emit_viz(
    arr: &ArrayRecord,
    base_w: usize,
    cfg: &DetectorConfig,
    viz_flags: &VizFlags,
) -> Result<()> {
    let bg = wrap::Background::compute(&arr.seq);
    let stats = wrap::wrap_and_ic(&arr.seq, base_w, &bg, cfg);
    let n_rows = stats.as_ref().map(|s| s.n_rows).unwrap_or(0);
    let column_ic_vec: Option<Vec<f64>> = stats.as_ref().map(|s| s.column_ic.clone());
    let embs = embed::embed_rows(&arr.seq, base_w, cfg);
    let ac = autocorr::compute(&embs, cfg.max_hor_k);
    let r_k_vec = ac.r_k.clone();
    let edge = if n_rows >= 2 {
        edges::compute(&arr.seq, base_w, n_rows)
    } else {
        None
    };
    let column_edge_rate_vec: Option<Vec<f64>> =
        edge.as_ref().map(|e| e.column_edge_rate.clone());
    let shift_feats = if n_rows >= 2 {
        shift::compute(&arr.seq, base_w, n_rows, cfg)
    } else {
        None
    };
    let best_shift_vec: Vec<i32> = shift_feats.map(|s| s.best_shift).unwrap_or_default();
    let bundle = viz::VizBundle {
        array_id: &arr.id,
        width_bp: base_w,
        seq: &arr.seq,
        n_rows,
        column_ic: column_ic_vec.as_deref(),
        column_edge_rate: column_edge_rate_vec.as_deref(),
        r_k: Some(&r_k_vec),
        best_shift: Some(&best_shift_vec),
    };
    viz::export(viz_flags, &bundle)
}

// M0 placeholder removed; M4 path produces real classifications.

/// M7.2: compute the analysis-block consensus-identity context.
/// Returns `Some` when ≥ 2 blocks were buildable and at least one
/// pairwise comparison cleared the coverage gate (i.e., the
/// context is usable for the mixed override + per-block segments).
///
/// Comparison-width choice follows `docs/new/detect_m7_plan.md` Q2:
///   - prefer `hor_length_bp` consensus when each analysis block fits
///     at least `cfg.min_complete_units_per_block` complete HOR units;
///   - otherwise fall back to `base_width_bp` with unit-aligned blocks.
fn compute_mixed_blocks_context(
    arr: &ArrayRecord,
    props: &Properties,
    cfg: &DetectorConfig,
) -> Option<segment::MixedBlocksContext> {
    let base_w = props.base_width_bp?;

    // Compare at base_width_bp. The plan (Q2) preferred
    // `hor_length_bp` to catch families that share monomer
    // composition but differ in slot arrangement. Empirically,
    // the v2 corpus's mixed cases use DIFFERENT monomer templates,
    // so base_width-level Hamming captures them. The hor_length_bp
    // comparison was also catastrophically confused by undetected
    // small phase shifts, which the detector's shift::compute
    // misses for short offsets — at hor_length_bp width, a
    // sub-monomer shift offsets the row alignment by a full
    // column-shift and identity collapses to ~0.25 even for a
    // single coherent HOR. base_width comparison is invariant to
    // multiples-of-base-width drift, which is the only kind a
    // single-family array can exhibit. Re-evaluate when corpora
    // contain HOR families with shared monomers but distinct slot
    // arrangements (currently out of M7 scope).
    let n_rows_base = arr.seq.len() / base_w.max(1);
    // Emit unit-aligned blocks at base_width when k is known so
    // boundaries land on slot edges.
    let unit_rows_in_base = props.hor_k.unwrap_or(1).max(1);
    let unit_rows_in_cmp = if unit_rows_in_base > 1 {
        Some(unit_rows_in_base)
    } else {
        None
    };
    let comparison_width = base_w;
    let n_rows_cmp = n_rows_base;

    // Phase-shift splits map from bp → rows in the comparison-width grid.
    let extra_splits: Vec<usize> = props
        .phase_shift_positions
        .iter()
        .filter_map(|&bp| {
            if comparison_width == 0 { return None; }
            Some(bp / comparison_width)
        })
        .collect();

    let blocks = analysis_blocks::build_blocks(
        n_rows_cmp,
        unit_rows_in_cmp,
        &extra_splits,
        cfg,
    );
    if blocks.len() < 2 {
        return None;
    }
    let mut consensuses =
        analysis_blocks::block_consensuses(&arr.seq, comparison_width, &blocks);
    // M7.2 calibration (2026-05-17): drop blocks whose per-block
    // column IC is below `ic_threshold_hor_unit`. These are
    // "unstructured" blocks (e.g., the foreign sequence inside a
    // hor_insertion case), and their consensus would otherwise look
    // ~random to the identity test and force a false-mixed call.
    let bg = wrap::Background::compute(&arr.seq);
    // M7.2 calibration: filter at `ic_threshold_hor_base` (0.30 default,
    // stricter than `ic_threshold_hor_unit`). Catches insertion blocks
    // whose IC at the chosen comparison width sits in the 0.2–0.4 range
    // without filtering out the mostly-clean blocks of cross-width
    // mixed cases.
    for (i, blk) in blocks.iter().enumerate() {
        if consensuses[i].is_none() {
            continue;
        }
        let ic = block_column_ic(
            &arr.seq,
            comparison_width,
            blk.start_row,
            blk.end_row,
            &bg,
        );
        if ic < cfg.ic_threshold_hor_base {
            consensuses[i] = None;
        }
    }
    let pairs = analysis_blocks::pairwise_identity(&consensuses, cfg);
    let medoid = analysis_blocks::pick_medoid(consensuses.len(), &pairs)?;
    Some(segment::MixedBlocksContext {
        blocks,
        consensuses,
        pairs,
        reference_block: medoid,
        comparison_width,
    })
}

/// M7.2: per-block mean column IC at the comparison width. Returns
/// 0.0 when the slice is too narrow to compute. Uses the array-wide
/// background frequencies passed in.
fn block_column_ic(
    seq: &[u8],
    width: usize,
    start_row: usize,
    end_row: usize,
    bg: &wrap::Background,
) -> f64 {
    let n_rows_total = seq.len() / width.max(1);
    let end = end_row.min(n_rows_total);
    if start_row >= end || width == 0 {
        return 0.0;
    }
    let n_rows = end - start_row;
    let mut total_ic = 0.0;
    let mut counted_cols = 0;
    for c in 0..width {
        let mut counts = [0usize; 4];
        let mut n_acgt = 0usize;
        for r in start_row..end {
            match seq[r * width + c] {
                b'A' => { counts[0] += 1; n_acgt += 1; }
                b'C' => { counts[1] += 1; n_acgt += 1; }
                b'G' => { counts[2] += 1; n_acgt += 1; }
                b'T' => { counts[3] += 1; n_acgt += 1; }
                _ => {}
            }
        }
        if n_acgt == 0 {
            continue;
        }
        let mut col_ic = 0.0;
        for i in 0..4 {
            let p = counts[i] as f64 / n_acgt as f64;
            if p > 0.0 {
                col_ic += p * (p / bg.q[i]).log2();
            }
        }
        total_ic += col_ic;
        counted_cols += 1;
        // `n_rows` unused but kept for parity with the whole-array
        // computation; per-row scaling already in p.
        let _ = n_rows;
    }
    if counted_cols == 0 {
        0.0
    } else {
        total_ic / counted_cols as f64
    }
}

/// M7.2: returns the worst-pair identity (lowest, most divergent),
/// or `None` if no valid pair survived the coverage gate.
fn worst_pair_identity(ctx: &segment::MixedBlocksContext) -> Option<(f64, usize, usize)> {
    ctx.pairs
        .iter()
        .min_by(|a, b| {
            a.identity
                .partial_cmp(&b.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|p| (p.identity, p.i, p.j))
}


/// M4 per-array work: produces a real `class` + supporting fields by
/// running the classify module over `width_features`. Then layers M3.5
/// Pass-B phase-shift offset recovery on top, using the chosen
/// `base_width_bp` as the primary width for shift analysis.
///
/// M7.2: also computes the analysis-block consensus-identity context.
/// When the override fires (HOR/IrregularHOR + min-pair identity
/// ≤ `stratification_diff_threshold`), the class is rewritten to
/// `Mixed`, class-defining fields are cleared (per review #1), and
/// the context is returned so per-block segment rows can be emitted
/// by the caller.
fn run_array_m4(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    cfg: &DetectorConfig,
) -> (Properties, Vec<WidthFeatures>, Option<segment::MixedBlocksContext>) {
    let (mut props_m35, widths) = run_array_m3_5(arr, pers, cfg);
    // DH4: capture the pre-decision width BEFORE the classify-result
    // copy overwrites props_m35.base_width_bp.
    let pre_decision_width = props_m35.base_width_bp;
    let decision = classify::decide_array(arr, pers, &widths, cfg);

    // Review-2026-05-16 #1: decision is authoritative for every
    // class-defining field. For Mixed/Ambiguous, the classifier
    // intentionally returns no base_width / k / IC / phase_sep —
    // we must NOT silently fall back to the M3.5 heuristic width
    // (publishing it produced misleading TSV rows like
    // `mixed,base=100,k=NA`). Keep shift / wobble fields from
    // M3.5 since they're computed independently and remain
    // informative even when classification fails.
    let is_resolved = matches!(
        decision.class,
        Class::SimpleTR | Class::HOR | Class::IrregularHOR
    );
    props_m35.class = decision.class;
    if is_resolved {
        props_m35.base_width_bp = decision.base_width_bp.or(props_m35.base_width_bp);
        props_m35.hor_k = decision.hor_k;
        props_m35.hor_length_bp = decision.hor_length_bp;
        props_m35.n_complete_copies = decision.n_complete_copies;
        props_m35.column_conservation = decision.column_conservation;
        props_m35.phase_separation = decision.phase_separation;
        props_m35.inter_monomer_identity = decision.inter_monomer_identity;
    } else {
        props_m35.base_width_bp = None;
        props_m35.hor_k = None;
        props_m35.hor_length_bp = None;
        props_m35.n_complete_copies = decision.n_complete_copies;
        props_m35.column_conservation = None;
        props_m35.phase_separation = None;
        props_m35.inter_monomer_identity = None;
    }
    props_m35.reason = decision.reason;

    // M7.2: same-width mixed override. Compute analysis-block
    // consensus identities for HOR / irregular_HOR arrays; if any
    // valid pair is at or below `stratification_diff_threshold`,
    // rewrite class to Mixed and clear class-defining fields.
    //
    // simple_TR is intentionally NOT eligible (M7 plan Q5).
    //
    // Single-family-perturbation gates (M7 plan risks §): pure
    // column-Hamming is confused by phase shifts (alignment offset),
    // wobble (drifting alignment), and mean-shift drift (gradual
    // misalignment over the array). All three would fire the mixed
    // override for what is structurally one family with alignment
    // perturbations. Skip the override when any of these signals is
    // strong. Inversion is documented as accepted false-mixed in
    // M7 plan §Risks (strand-aware deferred to v2).
    let mut mixed_ctx: Option<segment::MixedBlocksContext> = None;
    // M7.2 calibration: rely on the tightened `stratification_diff_threshold`
    // (0.50) plus best-alignment Hamming identity to discriminate
    // genuine mixed (random match rate ≈ 0.25-0.36) from single-family
    // shift/wobble cases (best-alignment identity 0.6+). Earlier
    // single_family_perturbed gates over-fired on mixed cases where
    // M3.5 detected spurious phase shifts in the confused wrap.
    if matches!(props_m35.class, Class::HOR | Class::IrregularHOR) {
        if let Some(ctx) = compute_mixed_blocks_context(arr, &props_m35, cfg) {
            if let Some((min_id, i, j)) = worst_pair_identity(&ctx) {
                log::debug!(
                    target: "kitehor::detect::analysis_blocks",
                    "{}: blocks={} cmp_w={} pairs={} min_id={:.3} medoid={}",
                    arr.id,
                    ctx.blocks.len(),
                    ctx.comparison_width,
                    ctx.pairs.len(),
                    min_id,
                    ctx.reference_block,
                );
                if min_id <= cfg.stratification_diff_threshold {
                    let n_blocks = ctx.blocks.len();
                    let original_class = props_m35.class.as_str();
                    props_m35.class = Class::Mixed;
                    props_m35.base_width_bp = None;
                    props_m35.hor_k = None;
                    props_m35.hor_length_bp = None;
                    props_m35.column_conservation = None;
                    props_m35.phase_separation = None;
                    props_m35.inter_monomer_identity = None;
                    props_m35.n_segments = n_blocks;
                    props_m35.reason = format!(
                        "mixed — block-consensus identity {:.3} ≤ diff_threshold {:.3} \
                         (block {} vs {} of {} at width {}; was {})",
                        min_id,
                        cfg.stratification_diff_threshold,
                        i, j, n_blocks,
                        ctx.comparison_width,
                        original_class,
                    );
                    mixed_ctx = Some(ctx);
                }
            }
        }
    }

    // DH3: copy irregularity_score from the chosen base_width's
    // width_features row. Demote HOR → irregular_HOR when the score
    // exceeds the calibrated threshold.
    //
    // Review-2026-05-16 #4: smooth wobble inflates block-level IC
    // variance because per-block row alignment drifts with the
    // wobble phase, even though the architecture is still a
    // coherent HOR with a wobble property. Suppress the demotion
    // when wobble_amplitude_bp dominates: large amplitude relative
    // to the base width AND no detected phase_shifts (which would
    // signal genuine architectural inconsistency).
    if let Some(bw) = props_m35.base_width_bp {
        if let Some(w) = widths.iter().find(|w| w.width_bp == bw) {
            props_m35.irregularity_score = w.irregularity_score;
            if matches!(props_m35.class, Class::HOR) {
                if let Some(irr) = w.irregularity_score {
                    if irr >= cfg.irregularity_demote_threshold {
                        // Wobble-dominance guard: high wobble (≥ 5%
                        // of base width) + no phase shifts → keep
                        // the HOR call, surface wobble as the
                        // explanation rather than demoting.
                        let wobble_frac = props_m35
                            .wobble_amplitude_bp
                            .map(|w_amp| w_amp.abs() / bw.max(1) as f64)
                            .unwrap_or(0.0);
                        let wobble_dominates = wobble_frac >= 0.05
                            && props_m35.n_phase_shifts == 0;
                        if !wobble_dominates {
                            props_m35.class = Class::IrregularHOR;
                            props_m35.reason = format!(
                                "{} (irregular_HOR — block-level IC variance {:.3} ≥ {:.3})",
                                props_m35.reason, irr, cfg.irregularity_demote_threshold
                            );
                        } else {
                            props_m35.reason = format!(
                                "{} (irregularity {:.3} attributed to wobble {:.1} bp / {} bp = {:.1}%; HOR retained)",
                                props_m35.reason,
                                irr,
                                props_m35.wobble_amplitude_bp.unwrap_or(0.0),
                                bw,
                                100.0 * wobble_frac,
                            );
                        }
                    }
                }
            }
        }
    }
    // M6: compute confidence over the populated property fields.
    props_m35.confidence = Some(confidence::compute(&props_m35, cfg));

    // If classification picked a different base_width than M3.5's
    // heuristic, rerun Pass-A/B at the new base_width to keep the
    // phase-shift positions/offsets consistent.
    if let Some(target_w) = decision.base_width_bp {
        let already_done = pre_decision_width == Some(target_w);
        if !already_done {
            if let Some(w_features) = widths.iter().find(|w| w.width_bp == target_w) {
                if w_features.rows >= 2 {
                    if let Some(shift_feats) =
                        shift::compute(&arr.seq, target_w, w_features.rows, cfg)
                    {
                        let window_rows = (cfg.block_size_rows_min / 8).max(8);
                        let offsets = shift::recover_offsets_at_breakpoints(
                            &arr.seq,
                            target_w,
                            w_features.rows,
                            &shift_feats.breakpoints,
                            window_rows,
                        );
                        let positions: Vec<usize> = shift_feats
                            .breakpoints
                            .iter()
                            .map(|&b| (b + 1) * target_w)
                            .collect();
                        props_m35.n_phase_shifts = positions.len();
                        props_m35.phase_shift_positions = positions;
                        props_m35.phase_shift_offsets =
                            offsets.iter().map(|&v| v as i64).collect();
                        props_m35.n_segments = 1 + props_m35.n_phase_shifts;
                        props_m35.mean_shift_bp = Some(shift_feats.mean_shift_bp);
                        props_m35.wobble_amplitude_bp = Some(shift_feats.wobble_amplitude_bp);
                        props_m35.wobble_periodicity_bp = shift_feats.wobble_periodicity_bp;
                    }
                }
            }
        }
    }

    // M7.2: the phase-shift recovery above runs from `decision.base_width_bp`
    // (the original HOR decision), so it would overwrite `n_segments` even
    // for Mixed arrays. Reassert the per-block segment count.
    if let Some(ctx) = &mixed_ctx {
        props_m35.n_segments = ctx.blocks.len();
    }

    (props_m35, widths, mixed_ctx)
}

/// M3.5 per-array work: M3 + Pass-B phase-shift offset recovery.
///
/// The "primary width" is the input period with the highest score
/// among those that produced valid Pass-A stats. (A heuristic until
/// M4's `phase::pick_best_width` lands.) Pass B runs at that width
/// to populate `Properties.n_phase_shifts`,
/// `phase_shift_positions`, and `phase_shift_offsets`.
fn run_array_m3_5(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    cfg: &DetectorConfig,
) -> (Properties, Vec<WidthFeatures>) {
    let (mut props, widths) = run_array_m3(arr, pers, cfg);

    // Find the highest-scored input period that produced valid Pass A
    // (i.e. width_features row has n_phase_shifts populated).
    let mut periods_sorted: Vec<&PeriodCandidate> = pers.iter().collect();
    periods_sorted.sort_by(|a, b| {
        b.period_score
            .partial_cmp(&a.period_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let primary = periods_sorted.iter().find_map(|p| {
        widths.iter().find(|w| w.width_bp == p.period_bp && w.rows >= 2)
    });
    if let Some(w) = primary {
        let width = w.width_bp;
        let n_rows = w.rows;
        if let Some(shift_feats) = shift::compute(&arr.seq, width, n_rows, cfg) {
            let window_rows = (cfg.block_size_rows_min / 8).max(8);
            let offsets = shift::recover_offsets_at_breakpoints(
                &arr.seq,
                width,
                n_rows,
                &shift_feats.breakpoints,
                window_rows,
            );
            // bp position of breakpoint at best_shift index `b` is
            // (b + 1) * width — the start of the post-shift row.
            let positions: Vec<usize> = shift_feats
                .breakpoints
                .iter()
                .map(|&b| (b + 1) * width)
                .collect();
            props.n_phase_shifts = positions.len();
            props.phase_shift_positions = positions;
            props.phase_shift_offsets = offsets.iter().map(|&v| v as i64).collect();
            props.n_segments = 1 + props.n_phase_shifts;
            props.mean_shift_bp = Some(shift_feats.mean_shift_bp);
            props.wobble_amplitude_bp = Some(shift_feats.wobble_amplitude_bp);
            props.wobble_periodicity_bp = shift_feats.wobble_periodicity_bp;
            props.base_width_bp = Some(width);
        }
    }
    (props, widths)
}

/// Width-feature builder: per-width compute used by M3+. Returns the
/// fully populated `WidthFeatures` row (DH6: includes phase_separation,
/// irregularity_score, and class_hint).
fn build_width_features(
    arr: &ArrayRecord,
    width: usize,
    bg: &wrap::Background,
    cfg: &DetectorConfig,
) -> WidthFeatures {
    let stats = wrap::wrap_and_ic(&arr.seq, width, bg, cfg);
    let (rows, ic, fc) = match &stats {
        Some(s) => (s.n_rows, Some(s.mean_column_ic), Some(s.fraction_conserved)),
        None => (0usize, None, None),
    };
    let mut row = WidthFeatures {
        array_id: arr.id.clone(),
        width_bp: width,
        rows,
        column_ic: ic,
        fraction_conserved_columns: fc,
        row_lag1_similarity: None,
        best_lag: None,
        best_lag_score: None,
        phase_separation: None,
        vertical_edge_rate: None,
        column_edge_autocorr_k: None,
        column_edge_autocorr_score: None,
        mean_shift_bp: None,
        wobble_amplitude_bp: None,
        n_phase_shifts: 0,
        irregularity_score: None,
        class_hint: ClassHint::UnsupportedWidth,
    };
    if stats.is_none() {
        return row;
    }

    // R(k) → phase_separation, primitive-corrected best lag.
    let embs = embed::embed_rows(&arr.seq, width, cfg);
    let ac = autocorr::compute(&embs, cfg.max_hor_k);
    row.row_lag1_similarity = ac.r_lag1;
    row.best_lag_score = ac.best_lag_score;
    let best_k_raw = ac.best_lag.unwrap_or(1);
    let best_k = if !ac.r_k.is_empty() {
        phase::primitive_correct(&ac.r_k, best_k_raw, cfg.primitive_correction_delta)
    } else {
        best_k_raw
    };
    row.best_lag = Some(best_k);
    if best_k >= 2 && !ac.r_k.is_empty() {
        row.phase_separation = Some(phase::phase_separation(&ac.r_k, best_k));
    } else {
        row.phase_separation = Some(0.0);
    }

    // Edge field.
    if rows >= 2 {
        if let Some(e) = edges::compute(&arr.seq, width, rows) {
            row.vertical_edge_rate = Some(e.vertical_edge_rate);
            row.column_edge_autocorr_k = e.column_edge_autocorr_k;
            row.column_edge_autocorr_score = e.column_edge_autocorr_score;
        }
    }

    // Pass-A shift signal.
    if rows >= 2 {
        if let Some(s) = shift::compute(&arr.seq, width, rows, cfg) {
            row.mean_shift_bp = Some(s.mean_shift_bp);
            row.wobble_amplitude_bp = Some(s.wobble_amplitude_bp);
            row.n_phase_shifts = s.breakpoints.len();
        }
    }

    // Irregularity (block-variance metric, A4 / DH3).
    row.irregularity_score = irregularity::compute(&arr.seq, width, bg, cfg);

    // Class hint per-width: cheap heuristic that lets calibration
    // see what each width "looks like" without re-running classify.
    row.class_hint = class_hint_for(&row, cfg);
    row
}

fn class_hint_for(w: &WidthFeatures, cfg: &DetectorConfig) -> ClassHint {
    let ic = w.column_ic.unwrap_or(0.0);
    if ic < cfg.ic_threshold_min || w.rows < cfg.min_rows_per_width {
        return ClassHint::UnsupportedWidth;
    }
    let r1 = w.row_lag1_similarity.unwrap_or(0.0);
    let phase = w.phase_separation.unwrap_or(0.0);
    let best_lag = w.best_lag.unwrap_or(1);
    if best_lag >= 2
        && phase >= cfg.phase_separation_threshold
        && r1 >= cfg.regime_c_r1_threshold
        && ic >= cfg.ic_threshold_hor_base
    {
        return ClassHint::HORBaseWidth { k: best_lag };
    }
    if r1 >= 0.95 && ic >= cfg.ic_threshold_hor_unit {
        return ClassHint::HORUnitWidth;
    }
    if ic >= cfg.ic_threshold_simple_tr && phase < cfg.phase_separation_threshold {
        return ClassHint::SimpleTRBaseWidth;
    }
    ClassHint::UnsupportedWidth
}

/// M3 per-array work: M2 + edge field + Pass-A shift signal.
/// DH6: now populates phase_separation, irregularity_score, and
/// class_hint via `build_width_features`.
fn run_array_m3(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    cfg: &DetectorConfig,
) -> (Properties, Vec<WidthFeatures>) {
    let bg = wrap::Background::compute(&arr.seq);
    let widths = widths::expand(pers, cfg, arr.length);
    let out: Vec<WidthFeatures> = widths
        .into_iter()
        .map(|w| build_width_features(arr, w, &bg, cfg))
        .collect();
    (Properties::placeholder(&arr.id, arr.length), out)
}

#[derive(Debug, Clone)]
pub struct DetectorReport {
    pub n_arrays: usize,
    pub n_segments: usize,
    pub n_width_rows: usize,
}

/// Batch over a directory of FASTAs + matching periods TSVs.
/// `<fasta_dir>/<stem>.fa` pairs with `<periods_dir>/<stem>.periods.tsv`.
pub fn run_batch(
    fasta_dir: &Path,
    periods_dir: &Path,
    out_dir: &Path,
    cfg: &DetectorConfig,
    viz_flags: &VizFlags,
    allow_missing_periods: bool,
    allow_extra_periods: bool,
) -> Result<usize> {
    use rayon::prelude::*;
    std::fs::create_dir_all(out_dir)?;
    let pairs = discover_pairs(fasta_dir, periods_dir, allow_extra_periods)?;
    let n = pairs.len();
    pairs.par_iter().try_for_each(|(fa, pe, stem)| -> Result<()> {
        let prefix = out_dir.join(stem);
        // Per-file invocations always pass `allow_extra_periods=true`
        // here because the batch loop already checked symmetry at the
        // directory level (`discover_pairs`).
        run_one(
            fa,
            pe,
            &prefix,
            cfg,
            viz_flags,
            allow_missing_periods,
            /*allow_extra_periods=*/ true,
        )
        .map(|_| ())
    })?;
    Ok(n)
}

fn discover_pairs(
    fasta_dir: &Path,
    periods_dir: &Path,
    allow_extra_periods: bool,
) -> Result<Vec<(std::path::PathBuf, std::path::PathBuf, String)>> {
    use std::collections::BTreeSet;
    let mut fasta_stems: BTreeSet<String> = BTreeSet::new();
    let mut periods_stems: BTreeSet<String> = BTreeSet::new();

    let mut out = Vec::new();
    for entry in std::fs::read_dir(fasta_dir)? {
        let e = entry?;
        let p = e.path();
        if p.extension().and_then(|x| x.to_str()) != Some("fa") {
            continue;
        }
        let stem = p
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow::anyhow!("non-utf8 stem in {:?}", p))?
            .to_string();
        fasta_stems.insert(stem.clone());
        let periods_path = periods_dir.join(format!("{stem}.periods.tsv"));
        if !periods_path.exists() {
            anyhow::bail!(
                "FASTA {:?} has no matching periods TSV at {:?}",
                p, periods_path
            );
        }
        out.push((p, periods_path, stem));
    }

    // DH11: symmetric pairing — periods TSVs without a matching FASTA
    // are misspelled or stale. Fail unless explicitly allowed.
    if !allow_extra_periods {
        for entry in std::fs::read_dir(periods_dir)? {
            let e = entry?;
            let p = e.path();
            let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if !name.ends_with(".periods.tsv") {
                continue;
            }
            let stem = name.trim_end_matches(".periods.tsv").to_string();
            periods_stems.insert(stem.clone());
        }
        let extras: Vec<&String> = periods_stems.difference(&fasta_stems).collect();
        if !extras.is_empty() {
            anyhow::bail!(
                "periods directory {:?} contains {} unmatched files: {:?}; \
                 pass `--allow-extra-periods` to ignore them",
                periods_dir, extras.len(), extras
            );
        }
    }

    out.sort();
    Ok(out)
}

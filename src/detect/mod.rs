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
        let (props, mut widths) = run_array_m4(arr, pers, cfg);
        // DH2: emit one Segment row per inter-shift region when
        // n_phase_shifts > 0.
        let mut new_segments = segment::split(&props);
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

/// M4 per-array work: produces a real `class` + supporting fields by
/// running the classify module over `width_features`. Then layers M3.5
/// Pass-B phase-shift offset recovery on top, using the chosen
/// `base_width_bp` as the primary width for shift analysis.
fn run_array_m4(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    cfg: &DetectorConfig,
) -> (Properties, Vec<WidthFeatures>) {
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

    (props_m35, widths)
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

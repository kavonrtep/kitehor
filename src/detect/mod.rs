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
) -> Result<DetectorReport> {
    cfg.validate()?;
    let arrays = io::load_arrays(fasta)?;
    let periods_by_id = io::load_periods(periods)?;
    let default_stem = fasta
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string());
    let paired = io::join_arrays_with_periods(arrays, periods_by_id, default_stem.as_deref())?;

    let mut properties: Vec<Properties> = Vec::with_capacity(paired.len());
    let segments: Vec<Segment> = Vec::new();
    let mut width_features: Vec<WidthFeatures> = Vec::new();
    let mut consensus_records: Vec<ConsensusRecord> = Vec::new();

    for (arr, pers) in &paired {
        let (props, mut widths) = run_array_m4(arr, pers, cfg);
        // M5: build consensus + viz when we have a chosen base width.
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
        properties.push(props);
        width_features.append(&mut widths);
    }

    io::write_properties(out_prefix, &properties)?;
    io::write_segments(out_prefix, &segments)?;
    io::write_width_features(out_prefix, &width_features)?;
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
    let decision = classify::decide_array(arr, pers, &widths, cfg);

    // Copy decision fields into properties.
    props_m35.class = decision.class;
    props_m35.base_width_bp = decision.base_width_bp.or(props_m35.base_width_bp);
    props_m35.hor_k = decision.hor_k;
    props_m35.hor_length_bp = decision.hor_length_bp;
    props_m35.n_complete_copies = decision.n_complete_copies;
    props_m35.column_conservation = decision.column_conservation;
    props_m35.phase_separation = decision.phase_separation;
    props_m35.inter_monomer_identity = decision.inter_monomer_identity;
    props_m35.reason = decision.reason;

    // If classification picked a different base_width than M3.5's
    // heuristic, rerun Pass-A/B at the new base_width to keep the
    // phase-shift positions/offsets consistent.
    if let Some(target_w) = decision.base_width_bp {
        let already_done = props_m35.base_width_bp == Some(target_w);
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

/// M3 per-array work: M2 + edge field + Pass-A shift signal.
/// Fills `vertical_edge_rate`, `column_edge_autocorr_*`,
/// `mean_shift_bp`, `wobble_amplitude_bp`, `n_phase_shifts`.
fn run_array_m3(
    arr: &ArrayRecord,
    pers: &[PeriodCandidate],
    cfg: &DetectorConfig,
) -> (Properties, Vec<WidthFeatures>) {
    let bg = wrap::Background::compute(&arr.seq);
    let widths = widths::expand(pers, cfg, arr.length);
    let mut out = Vec::with_capacity(widths.len());
    for w in widths {
        let stats = wrap::wrap_and_ic(&arr.seq, w, &bg, cfg);
        let (rows, ic, fc) = match &stats {
            Some(s) => (s.n_rows, Some(s.mean_column_ic), Some(s.fraction_conserved)),
            None => (0usize, None, None),
        };
        let (r_lag1, best_lag, best_lag_score) = if stats.is_some() {
            let embs = embed::embed_rows(&arr.seq, w, cfg);
            let summary = autocorr::compute(&embs, cfg.max_hor_k);
            (summary.r_lag1, summary.best_lag, summary.best_lag_score)
        } else {
            (None, None, None)
        };
        let edge = if rows >= 2 {
            edges::compute(&arr.seq, w, rows)
        } else {
            None
        };
        let shift = if rows >= 2 {
            shift::compute(&arr.seq, w, rows, cfg)
        } else {
            None
        };
        let (vertical_edge_rate, column_edge_autocorr_k, column_edge_autocorr_score) =
            match &edge {
                Some(e) => (
                    Some(e.vertical_edge_rate),
                    e.column_edge_autocorr_k,
                    e.column_edge_autocorr_score,
                ),
                None => (None, None, None),
            };
        let (mean_shift_bp, wobble_amplitude_bp, n_phase_shifts) = match &shift {
            Some(s) => (
                Some(s.mean_shift_bp),
                Some(s.wobble_amplitude_bp),
                s.breakpoints.len(),
            ),
            None => (None, None, 0),
        };
        out.push(WidthFeatures {
            array_id: arr.id.clone(),
            width_bp: w,
            rows,
            column_ic: ic,
            fraction_conserved_columns: fc,
            row_lag1_similarity: r_lag1,
            best_lag,
            best_lag_score,
            phase_separation: None, // M4
            vertical_edge_rate,
            column_edge_autocorr_k,
            column_edge_autocorr_score,
            mean_shift_bp,
            wobble_amplitude_bp,
            n_phase_shifts,
            irregularity_score: None, // M4
            class_hint: ClassHint::UnsupportedWidth, // M4
        });
    }
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
) -> Result<usize> {
    use rayon::prelude::*;
    std::fs::create_dir_all(out_dir)?;
    let pairs = discover_pairs(fasta_dir, periods_dir)?;
    let n = pairs.len();
    pairs.par_iter().try_for_each(|(fa, pe, stem)| -> Result<()> {
        let prefix = out_dir.join(stem);
        run_one(fa, pe, &prefix, cfg, viz_flags).map(|_| ())
    })?;
    Ok(n)
}

fn discover_pairs(fasta_dir: &Path, periods_dir: &Path) -> Result<Vec<(std::path::PathBuf, std::path::PathBuf, String)>> {
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
        let periods_path = periods_dir.join(format!("{stem}.periods.tsv"));
        if !periods_path.exists() {
            anyhow::bail!(
                "FASTA {:?} has no matching periods TSV at {:?}",
                p, periods_path
            );
        }
        out.push((p, periods_path, stem));
    }
    out.sort();
    Ok(out)
}

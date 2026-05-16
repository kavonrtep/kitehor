//! Shared data types for the line-width detector.
//!
//! Schema for `properties.tsv` is **frozen at M0** per the
//! implementation plan §10.1 — every later milestone fills more
//! columns with real values but no column is added or removed.

use serde::{Deserialize, Serialize};

/// One of the five summary class labels (`detect_spec.md` §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Class {
    #[serde(rename = "simple_TR")]
    SimpleTR,
    HOR,
    #[serde(rename = "irregular_HOR")]
    IrregularHOR,
    Mixed,
    Ambiguous,
}

impl Class {
    pub fn as_str(&self) -> &'static str {
        match self {
            Class::SimpleTR => "simple_TR",
            Class::HOR => "HOR",
            Class::IrregularHOR => "irregular_HOR",
            Class::Mixed => "mixed",
            Class::Ambiguous => "ambiguous",
        }
    }
}

/// Class hint emitted by a single tested width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassHint {
    SimpleTRBaseWidth,
    HORBaseWidth { k: usize },
    HORUnitWidth,
    UnsupportedWidth,
}

impl ClassHint {
    pub fn as_str(&self) -> &'static str {
        match self {
            ClassHint::SimpleTRBaseWidth => "simple_TR_base",
            ClassHint::HORBaseWidth { .. } => "HOR_base",
            ClassHint::HORUnitWidth => "HOR_unit",
            ClassHint::UnsupportedWidth => "unsupported",
        }
    }
}

/// Final per-array property vector. See `detect_impl_plan.md §10.1`
/// for the frozen TSV schema.
#[derive(Debug, Clone)]
pub struct Properties {
    pub array_id: String,
    pub length_bp: usize,
    pub class: Class,
    pub base_width_bp: Option<usize>,
    pub hor_k: Option<usize>,
    pub hor_length_bp: Option<usize>,
    pub n_complete_copies: Option<usize>,
    pub column_conservation: Option<f64>,
    pub phase_separation: Option<f64>,
    pub mean_shift_bp: Option<f64>,
    pub wobble_amplitude_bp: Option<f64>,
    pub wobble_periodicity_bp: Option<f64>,
    pub n_phase_shifts: usize,
    pub phase_shift_positions: Vec<usize>,
    pub phase_shift_offsets: Vec<i64>,
    pub irregularity_score: Option<f64>,
    /// **Approximation.** Review-2026-05-16 #5: in the current
    /// implementation this column carries `R(1)` at the chosen base
    /// width — k-mer-composition row similarity, not pairwise
    /// sequence identity between inferred slot consensuses. It is
    /// useful as a regime indicator (≥ 0.95 in regime B, lower in
    /// regime C) but should NOT be interpreted as a calibrated
    /// biological identity. Schema is frozen at M0 so the field
    /// name stays; a future major version will either rename it
    /// or fill it with a real mean-pairwise-identity calculation.
    pub inter_monomer_identity: Option<f64>,
    /// Heuristic confidence score, not a calibrated probability.
    /// Review-2026-05-16 #6.
    pub confidence: Option<f64>,
    pub n_segments: usize,
    pub reason: String,
}

impl Properties {
    /// Build a "no detection yet" placeholder for M0. Every property
    /// is empty / NA except `array_id`, `length_bp`, and the
    /// `Ambiguous` class default.
    pub fn placeholder(array_id: &str, length_bp: usize) -> Self {
        Self {
            array_id: array_id.to_string(),
            length_bp,
            class: Class::Ambiguous,
            base_width_bp: None,
            hor_k: None,
            hor_length_bp: None,
            n_complete_copies: None,
            column_conservation: None,
            phase_separation: None,
            mean_shift_bp: None,
            wobble_amplitude_bp: None,
            wobble_periodicity_bp: None,
            n_phase_shifts: 0,
            phase_shift_positions: Vec::new(),
            phase_shift_offsets: Vec::new(),
            irregularity_score: None,
            inter_monomer_identity: None,
            confidence: None,
            n_segments: 1,
            reason: "M0 scaffolding — detection logic not yet wired".to_string(),
        }
    }
}

/// Frozen header for `properties.tsv` (20 columns, A4 — includes
/// `inter_monomer_identity`).
pub const PROPERTIES_HEADER: &str = "\
array_id\tlength_bp\tclass\tbase_width_bp\thor_k\thor_length_bp\
\tn_complete_copies\tcolumn_conservation\tphase_separation\tmean_shift_bp\
\twobble_amplitude_bp\twobble_periodicity_bp\tn_phase_shifts\
\tphase_shift_positions\tphase_shift_offsets\tirregularity_score\
\tinter_monomer_identity\tconfidence\tn_segments\treason";

/// Frozen header for `segments.tsv` (11 columns).
pub const SEGMENTS_HEADER: &str = "\
array_id\tsegment_id\tstart_bp\tend_bp\tclass\tbase_width_bp\thor_k\
\tcolumn_conservation\tphase_separation\twobble_amplitude_bp\
\tirregularity_score";

/// Frozen header for `width_features.tsv` (17 columns).
pub const WIDTH_FEATURES_HEADER: &str = "\
array_id\twidth_bp\trows\tcolumn_IC\tfraction_conserved_columns\
\trow_lag1_similarity\tbest_lag\tbest_lag_score\tphase_separation\
\tvertical_edge_rate\tcolumn_edge_autocorr_k\tcolumn_edge_autocorr_score\
\tmean_shift_bp\twobble_amplitude_bp\tn_phase_shifts\tirregularity_score\
\tclass_hint";

/// Per-segment row for `segments.tsv`. Populated when
/// `n_segments > 1` (M4).
#[derive(Debug, Clone)]
pub struct Segment {
    pub array_id: String,
    pub segment_id: usize,
    pub start_bp: usize,
    pub end_bp: usize,
    pub class: Class,
    pub base_width_bp: Option<usize>,
    pub hor_k: Option<usize>,
    pub column_conservation: Option<f64>,
    pub phase_separation: Option<f64>,
    pub wobble_amplitude_bp: Option<f64>,
    pub irregularity_score: Option<f64>,
}

/// Per-width diagnostic row for `width_features.tsv`. M0 emits one
/// placeholder per (array, tested-width) once width expansion lands
/// in M1; in M0 the file is header-only.
#[derive(Debug, Clone)]
pub struct WidthFeatures {
    pub array_id: String,
    pub width_bp: usize,
    pub rows: usize,
    pub column_ic: Option<f64>,
    pub fraction_conserved_columns: Option<f64>,
    pub row_lag1_similarity: Option<f64>,
    pub best_lag: Option<usize>,
    pub best_lag_score: Option<f64>,
    pub phase_separation: Option<f64>,
    pub vertical_edge_rate: Option<f64>,
    pub column_edge_autocorr_k: Option<usize>,
    pub column_edge_autocorr_score: Option<f64>,
    pub mean_shift_bp: Option<f64>,
    pub wobble_amplitude_bp: Option<f64>,
    pub n_phase_shifts: usize,
    pub irregularity_score: Option<f64>,
    pub class_hint: ClassHint,
}

/// A single period candidate, as emitted by `kitehor synth` and any
/// equivalent upstream period generator. Schema spelled out in
/// `detect_impl_plan.md §7` (A12).
#[derive(Debug, Clone)]
pub struct PeriodCandidate {
    pub array_id: String,
    pub period_bp: usize,
    pub period_score: f64,
    pub source: String,
}

//! `truth.tsv` writer.
//!
//! One row per array, mirroring the schema in `docs/new/taxonomy.md §5.2`.
//! Columns whose values are unknown until later milestones (wobble,
//! events) are emitted as NA in M4 and back-filled by callers.

use crate::synth::blocks::{FillerKind, SimState};
use crate::synth::config::{Block, Config, Template};
use crate::synth::events::{to_events_json, EventLog};
use crate::synth::fasta::with_ext;
use crate::synth::grammar::to_expression;
use crate::synth::noise::NoiseLog;
use crate::synth::wobble::WobbleLog;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::Path;

/// Header of `truth.tsv`. Exposed so tests can assert against it.
pub const TRUTH_HEADER: &str =
    "array_id\tlength_bp\ttruth_class\tbase_width_bp\thor_k\thor_length_bp\
\tn_complete_copies\twobble_amplitude_bp\twobble_periodicity_bp\
\tn_phase_shifts\tphase_shift_positions\tphase_shift_offsets\
\tn_segments\treason\tstructural_expression\tschema_version\tevents_json";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TruthClass {
    SimpleTR,
    HOR,
    Mixed,
    Random,
}

impl TruthClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            TruthClass::SimpleTR => "simple_TR",
            TruthClass::HOR => "HOR",
            TruthClass::Mixed => "mixed",
            TruthClass::Random => "random",
        }
    }
}

#[derive(Debug, Clone)]
pub struct TruthRow {
    pub array_id: String,
    pub length_bp: usize,
    pub truth_class: TruthClass,
    pub base_width_bp: usize,
    pub hor_k: Option<usize>,
    pub hor_length_bp: Option<usize>,
    pub n_complete_copies: usize,
    pub wobble_amplitude_bp: f64,
    pub wobble_periodicity_bp: Option<f64>,
    pub n_phase_shifts: usize,
    pub phase_shift_positions: Vec<usize>,
    pub phase_shift_offsets: Vec<i64>,
    pub n_segments: usize,
    pub reason: String,
    pub structural_expression: String,
    pub schema_version: u32,
    pub events_json: String,
}

impl TruthRow {
    pub fn to_tsv(&self) -> String {
        let na = "NA";
        let fmt_opt_usize =
            |o: &Option<usize>| o.map(|x| x.to_string()).unwrap_or_else(|| na.to_string());
        let fmt_opt_f64 = |o: &Option<f64>| {
            o.map(|x| format!("{:.4}", x))
                .unwrap_or_else(|| na.to_string())
        };
        let fmt_list_usize = |v: &[usize]| {
            if v.is_empty() {
                na.to_string()
            } else {
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
        };
        let fmt_list_i64 = |v: &[i64]| {
            if v.is_empty() {
                na.to_string()
            } else {
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.array_id,
            self.length_bp,
            self.truth_class.as_str(),
            self.base_width_bp,
            fmt_opt_usize(&self.hor_k),
            fmt_opt_usize(&self.hor_length_bp),
            self.n_complete_copies,
            self.wobble_amplitude_bp,
            fmt_opt_f64(&self.wobble_periodicity_bp),
            self.n_phase_shifts,
            fmt_list_usize(&self.phase_shift_positions),
            fmt_list_i64(&self.phase_shift_offsets),
            self.n_segments,
            self.reason,
            self.structural_expression,
            self.schema_version,
            self.events_json,
        )
    }
}

pub fn build_truth(
    cfg: &Config,
    state: &SimState,
    array_id: &str,
    _noise: &NoiseLog,
    wobble: &WobbleLog,
    events: &[EventLog],
) -> TruthRow {
    let (truth_class, base_width_bp, hor_k, n_complete_copies, reason) = classify(cfg);
    let hor_length_bp = hor_k.map(|k| base_width_bp * k);

    let (phase_shift_positions, phase_shift_offsets) = collect_phase_shifts(state);
    let n_phase_shifts = phase_shift_offsets.len();

    TruthRow {
        array_id: array_id.to_string(),
        length_bp: state.sequence.len(),
        truth_class,
        base_width_bp,
        hor_k,
        hor_length_bp,
        n_complete_copies,
        wobble_amplitude_bp: wobble.realised_amplitude_bp,
        wobble_periodicity_bp: wobble.realised_periodicity_bp,
        n_phase_shifts,
        phase_shift_positions,
        phase_shift_offsets,
        n_segments: 1 + n_phase_shifts,
        reason,
        structural_expression: to_expression(cfg),
        schema_version: cfg.schema_version,
        events_json: to_events_json(events),
    }
}

/// Walk `state.filler_spans` and pull out the SHIFT entries.
fn collect_phase_shifts(state: &SimState) -> (Vec<usize>, Vec<i64>) {
    let mut positions = Vec::new();
    let mut offsets = Vec::new();
    for fs in &state.filler_spans {
        if let FillerKind::Shift { offset_bp } = fs.kind {
            positions.push(fs.realised_start_bp);
            offsets.push(offset_bp);
        }
    }
    (positions, offsets)
}

/// Decide the array's truth class from the structure.
///
/// MVP rules:
/// - Exactly one HOR block (with k≥2, div>0) → HOR.
/// - Exactly one SIMPLE_TR block (or HOR block with div==0) → simple_TR.
/// - Multiple repeat blocks, all with the same base_width and k → HOR
///   or simple_TR (phase shifts / insertions between identical blocks
///   are properties, not class changes — see plan A5).
/// - Multiple repeat blocks with different base_width or k → mixed.
/// - Random insertion only → random (degenerate; mostly for T17).
fn classify(cfg: &Config) -> (TruthClass, usize, Option<usize>, usize, String) {
    let repeats: Vec<&Block> = cfg
        .structure
        .iter()
        .filter(|b| matches!(b, Block::HOR { .. } | Block::SIMPLE_TR { .. }))
        .collect();

    if repeats.is_empty() {
        return (TruthClass::Random, 0, None, 0, "no repeat blocks".into());
    }

    let descriptions: Vec<(usize, Option<usize>, usize)> =
        repeats.iter().map(|b| describe_repeat(cfg, b)).collect();
    let first = descriptions[0];
    let mismatch = descriptions
        .iter()
        .any(|d| d.0 != first.0 || d.1 != first.1);

    if mismatch {
        let total_copies: usize = descriptions.iter().map(|d| d.2).sum();
        return (
            TruthClass::Mixed,
            first.0,
            first.1,
            total_copies,
            "multiple repeat blocks with different base_width or k".into(),
        );
    }

    let total_copies: usize = descriptions.iter().map(|d| d.2).sum();
    let class = if first.1.is_some() {
        TruthClass::HOR
    } else {
        TruthClass::SimpleTR
    };
    let reason = match repeats.len() {
        1 => "single repeat block".into(),
        n => format!("{} identical-architecture blocks (treated as one)", n),
    };
    (class, first.0, first.1, total_copies, reason)
}

/// Returns (base_width_bp, Option<hor_k>, n_copies).
fn describe_repeat(cfg: &Config, b: &Block) -> (usize, Option<usize>, usize) {
    match b {
        Block::HOR { template, n_copies } => {
            let (l, k, div) = match cfg.templates.get(template) {
                Some(Template::HOR_slots {
                    monomer_length_bp,
                    k,
                    inter_slot_divergence,
                    ..
                }) => (*monomer_length_bp, *k, *inter_slot_divergence),
                _ => (0, 1, 0.0),
            };
            // div==0 collapses HOR to simple_TR at base_width=L.
            if div == 0.0 {
                (l, None, *n_copies * k)
            } else {
                (l, Some(k), *n_copies)
            }
        }
        Block::SIMPLE_TR { template, n_copies } => {
            let l = match cfg.templates.get(template) {
                Some(Template::HOR_slots {
                    monomer_length_bp, ..
                })
                | Some(Template::monomer {
                    monomer_length_bp, ..
                }) => *monomer_length_bp,
                None => 0,
            };
            (l, None, *n_copies)
        }
        _ => (0, None, 0),
    }
}

pub fn write(prefix: &Path, row: &TruthRow) -> Result<()> {
    let path = with_ext(prefix, "truth.tsv");
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "{}", TRUTH_HEADER)?;
    writeln!(f, "{}", row.to_tsv())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::{blocks::expand, rng::Streams, templates::instantiate};
    use std::io::Write;

    fn parse(yaml: &str) -> Config {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        crate::synth::config::load_and_validate(f.path()).unwrap()
    }

    fn full_state(cfg: &Config) -> SimState {
        let s = Streams::new(cfg.seed);
        let mut rt = s.templates();
        let inst = instantiate(&cfg.templates, &mut rt);
        let mut rs = s.structure();
        expand(&cfg.structure, &inst, &mut rs).unwrap()
    }

    #[test]
    fn single_hor_classification() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 100
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr1",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        assert_eq!(row.truth_class, TruthClass::HOR);
        assert_eq!(row.base_width_bp, 171);
        assert_eq!(row.hor_k, Some(12));
        assert_eq!(row.hor_length_bp, Some(171 * 12));
        assert_eq!(row.n_complete_copies, 100);
        assert_eq!(row.n_phase_shifts, 0);
        assert_eq!(row.n_segments, 1);
    }

    #[test]
    fn simple_tr_classification() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 170
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 1000
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr1",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        assert_eq!(row.truth_class, TruthClass::SimpleTR);
        assert_eq!(row.base_width_bp, 170);
        assert!(row.hor_k.is_none());
        assert_eq!(row.n_complete_copies, 1000);
    }

    #[test]
    fn hor_with_zero_divergence_is_simple_tr() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.0
structure:
  - type: HOR
    template: t
    n_copies: 50
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr1",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        assert_eq!(row.truth_class, TruthClass::SimpleTR);
        assert!(row.hor_k.is_none());
        // Total copies = 50 * k = 200 monomers.
        assert_eq!(row.n_complete_copies, 200);
    }

    #[test]
    fn phase_shifted_hor_is_hor_with_one_shift() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.10
structure:
  - type: HOR
    template: t
    n_copies: 50
  - type: SHIFT
    offset_bp: 25
  - type: HOR
    template: t
    n_copies: 50
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr1",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        assert_eq!(row.truth_class, TruthClass::HOR);
        assert_eq!(row.n_phase_shifts, 1);
        assert_eq!(row.phase_shift_offsets, vec![25]);
        assert_eq!(row.phase_shift_positions, vec![100 * 4 * 50]);
        assert_eq!(row.n_segments, 2);
    }

    #[test]
    fn mixed_blocks_classified_as_mixed() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  a:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.15
  b:
    type: HOR_slots
    monomer_length_bp: 200
    k: 8
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: a
    n_copies: 50
  - type: HOR
    template: b
    n_copies: 50
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr1",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        assert_eq!(row.truth_class, TruthClass::Mixed);
    }

    #[test]
    fn tsv_row_roundtrips() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.10
structure:
  - type: HOR
    template: t
    n_copies: 5
"#,
        );
        let st = full_state(&cfg);
        let row = build_truth(
            &cfg,
            &st,
            "arr_test",
            &NoiseLog::default(),
            &WobbleLog::default(),
            &[],
        );
        let line = row.to_tsv();
        let fields: Vec<&str> = line.split('\t').collect();
        let header_fields: Vec<&str> = TRUTH_HEADER.split('\t').collect();
        assert_eq!(fields.len(), header_fields.len(), "field count mismatch");
        assert_eq!(fields[0], "arr_test");
        // length_bp = 5 * 4 * 100 = 2000
        assert_eq!(fields[1], "2000");
        assert_eq!(fields[2], "HOR");
    }
}

//! YAML config types + loader + validation.
//!
//! Validation layers (in order):
//!
//! 1. **JSON Schema** — `jsonschema` validates the parsed YAML
//!    against the canonical embedded schema (`docs/new/simulator_schema.json`).
//!    Catches type errors, range/pattern/minimum constraints, and
//!    `oneOf` discriminated-union mismatches. This is the same schema
//!    `synth-schema --print` emits, so `synth-validate` and the
//!    contract are guaranteed to agree.
//! 2. **Structural** — serde with `deny_unknown_fields` everywhere
//!    deserialises into the Rust types and catches anything the schema
//!    missed.
//! 3. **MVP business rules** — `validate_mvp_invariants` enforces the
//!    contract items the schema can't express:
//!    - **A1**: every `post_generation` event names a `block` index
//!      that points to a HOR/SIMPLE_TR block and whose copy range
//!      fits within that block's `n_copies`.
//!    - **A3**: `global.output` is silently ignored (non-fatal warn).
//!    - **Q5**: a negative `SHIFT` block must follow a HOR or
//!      SIMPLE_TR and `|offset_bp|` may not exceed
//!      `monomer_length / 2`.
//!    - **Q8**: `source: file` is rejected as not-implemented-in-MVP.
//!    - **F4**: when `source: sequence`, the inline string length must
//!      equal `monomer_length_bp` (and contain only ACGT).

use serde::Deserialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("reading {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("parsing {path}: {source}")]
    Yaml {
        path: PathBuf,
        #[source]
        source: serde_yaml::Error,
    },
    #[error("{0}")]
    Invariant(String),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub schema_version: u32,
    #[serde(default)]
    pub seed: u64,
    #[serde(default)]
    pub global: Global,
    #[serde(default)]
    pub templates: HashMap<String, Template>,
    pub structure: Vec<Block>,
    #[serde(default)]
    pub modifiers: Vec<Modifier>,
    #[serde(default)]
    pub post_generation: Vec<Event>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Global {
    #[serde(default)]
    pub mutation_rate: f64,
    #[serde(default)]
    pub indel_rate: f64,
    /// A3: tolerated in YAML for backwards-compat with the upstream
    /// `simulator_plan.md` but ignored at runtime. The CLI `-o/--out`
    /// flag is the single source of truth for output paths.
    #[serde(default)]
    pub output: Option<String>,
    #[serde(default)]
    pub array_id: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
#[allow(non_camel_case_types)]
pub enum Template {
    HOR_slots {
        monomer_length_bp: usize,
        k: usize,
        #[serde(default = "default_source")]
        source: Source,
        #[serde(default)]
        sequence: Option<String>,
        #[serde(default)]
        file: Option<PathBuf>,
        #[serde(default = "default_gc")]
        gc_content: f64,
        #[serde(default)]
        inter_slot_divergence: f64,
    },
    monomer {
        monomer_length_bp: usize,
        #[serde(default = "default_source")]
        source: Source,
        #[serde(default)]
        sequence: Option<String>,
        #[serde(default)]
        file: Option<PathBuf>,
        #[serde(default = "default_gc")]
        gc_content: f64,
    },
}

#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Source {
    Random,
    Sequence,
    File,
}

fn default_source() -> Source {
    Source::Random
}
fn default_gc() -> f64 {
    0.5
}
fn default_split() -> f64 {
    0.5
}
fn default_target() -> String {
    "all".to_string()
}
fn default_wobble_model() -> WobbleModel {
    WobbleModel::RandomWalk
}
fn default_realisation() -> WobbleRealisation {
    WobbleRealisation::IntegerEdits
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
#[allow(non_camel_case_types)]
pub enum Block {
    HOR {
        template: String,
        n_copies: usize,
    },
    SIMPLE_TR {
        template: String,
        n_copies: usize,
    },
    SHIFT {
        offset_bp: i64,
    },
    INSERTION {
        length_bp: usize,
        kind: InsertionKind,
    },
}

#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
pub enum InsertionKind {
    #[serde(rename = "random")]
    Random,
    #[serde(rename = "AT_rich")]
    AtRich,
    #[serde(rename = "GC_rich")]
    GcRich,
    #[serde(rename = "retro_like")]
    RetroLike,
    #[serde(rename = "segdup_like")]
    SegdupLike,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Modifier {
    #[serde(default = "default_target")]
    pub target: String,
    pub wobble: WobbleSpec,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WobbleSpec {
    pub amplitude_bp: f64,
    #[serde(default)]
    pub period_rows: u32,
    #[serde(default = "default_wobble_model")]
    pub model: WobbleModel,
    #[serde(default = "default_realisation")]
    pub realisation: WobbleRealisation,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum WobbleModel {
    Sinusoidal,
    RandomWalk,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum WobbleRealisation {
    IntegerEdits,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
#[allow(non_camel_case_types)]
pub enum Event {
    HYBRID {
        /// A1: 0-indexed into `structure`; target must be HOR or SIMPLE_TR.
        block: usize,
        /// 1-indexed within the targeted block.
        at_copy: usize,
        slot: usize,
        source_slots: [usize; 2],
        #[serde(default = "default_split")]
        split_fraction: f64,
    },
    INVERSION {
        block: usize,
        start_copy: usize,
        length_copies: usize,
    },
    DUPLICATION {
        block: usize,
        start_copy: usize,
        length_copies: usize,
    },
    DELETION {
        block: usize,
        start_copy: usize,
        length_copies: usize,
    },
}

// ---------- public API ----------

pub fn load_and_validate(path: &Path) -> Result<Config, ConfigError> {
    let text = std::fs::read_to_string(path).map_err(|e| ConfigError::Io {
        path: path.to_owned(),
        source: e,
    })?;

    // 1. JSON Schema validation. Parse YAML into a `serde_yaml::Value`
    //    first so we can also surface any parse error early; then
    //    cross-walk to `serde_json::Value` for the schema engine.
    let yaml_value: serde_yaml::Value =
        serde_yaml::from_str(&text).map_err(|e| ConfigError::Yaml {
            path: path.to_owned(),
            source: e,
        })?;
    let json_value: serde_json::Value = serde_json::to_value(&yaml_value)
        .map_err(|e| ConfigError::Invariant(format!("YAML→JSON conversion failed: {e}")))?;
    let schema: serde_json::Value = serde_json::from_str(crate::synth::CANONICAL_SCHEMA)
        .expect("embedded simulator.schema.json must be valid JSON");
    let validator = jsonschema::validator_for(&schema)
        .map_err(|e| ConfigError::Invariant(format!("schema compile failed: {e}")))?;
    let errors: Vec<String> = validator
        .iter_errors(&json_value)
        .map(|e| format!("  at `{}`: {}", e.instance_path, e))
        .collect();
    if !errors.is_empty() {
        return Err(ConfigError::Invariant(format!(
            "JSON Schema validation failed ({} error{}):\n{}",
            errors.len(),
            if errors.len() == 1 { "" } else { "s" },
            errors.join("\n")
        )));
    }

    // 2. Typed deserialisation.
    let cfg: Config = serde_yaml::from_value(yaml_value).map_err(|e| ConfigError::Yaml {
        path: path.to_owned(),
        source: e,
    })?;
    if cfg.schema_version != 1 {
        return Err(ConfigError::Invariant(format!(
            "schema_version must be 1 (got {})",
            cfg.schema_version
        )));
    }

    // 3. MVP business rules (A1, A3, Q5, Q8, F4).
    validate_mvp_invariants(&cfg)?;
    Ok(cfg)
}

fn validate_mvp_invariants(cfg: &Config) -> Result<(), ConfigError> {
    // Q8: source:file is not implemented in MVP. Reject explicitly
    // rather than silently letting generation fail later.
    // F4: source:sequence requires an inline string of length equal to
    //     monomer_length_bp (after uppercasing, A/C/G/T only — pattern
    //     enforced by the schema).
    for (name, tpl) in &cfg.templates {
        let (src, seq, monomer_length) = match tpl {
            Template::HOR_slots {
                source,
                sequence,
                monomer_length_bp,
                ..
            } => (*source, sequence.as_deref(), *monomer_length_bp),
            Template::monomer {
                source,
                sequence,
                monomer_length_bp,
                ..
            } => (*source, sequence.as_deref(), *monomer_length_bp),
        };
        if src == Source::File {
            return Err(ConfigError::Invariant(format!(
                "template `{}`: source: file is not implemented in MVP — use 'random' or 'sequence'",
                name
            )));
        }
        if src == Source::Sequence {
            let s = seq.ok_or_else(|| {
                ConfigError::Invariant(format!(
                    "template `{}`: source: sequence requires a `sequence` field",
                    name
                ))
            })?;
            if s.len() != monomer_length {
                return Err(ConfigError::Invariant(format!(
                    "template `{}`: source: sequence length {} must equal monomer_length_bp {}",
                    name,
                    s.len(),
                    monomer_length
                )));
            }
        }
    }

    // A3: warn (non-fatal) if global.output is set.
    if cfg.global.output.is_some() {
        log::warn!("global.output is ignored in MVP; use the CLI -o/--out flag");
    }

    // A1: validate every post_generation event names a real HOR/SIMPLE_TR
    // block and that its copy range fits.
    for (i, ev) in cfg.post_generation.iter().enumerate() {
        let (block_idx, copy_start, copy_end) = match ev {
            Event::HYBRID { block, at_copy, .. } => (*block, *at_copy, *at_copy),
            Event::INVERSION {
                block,
                start_copy,
                length_copies,
            }
            | Event::DUPLICATION {
                block,
                start_copy,
                length_copies,
            }
            | Event::DELETION {
                block,
                start_copy,
                length_copies,
            } => {
                if *length_copies == 0 {
                    return Err(ConfigError::Invariant(format!(
                        "post_generation[{}]: length_copies must be >= 1",
                        i
                    )));
                }
                (*block, *start_copy, *start_copy + *length_copies - 1)
            }
        };
        if block_idx >= cfg.structure.len() {
            return Err(ConfigError::Invariant(format!(
                "post_generation[{}]: block index {} out of range (structure has {} blocks)",
                i,
                block_idx,
                cfg.structure.len()
            )));
        }
        let (n_copies, kind) = match &cfg.structure[block_idx] {
            Block::HOR { n_copies, .. } => (*n_copies, "HOR"),
            Block::SIMPLE_TR { n_copies, .. } => (*n_copies, "SIMPLE_TR"),
            Block::SHIFT { .. } => {
                return Err(ConfigError::Invariant(format!(
                    "post_generation[{}]: target block {} is a SHIFT, not HOR/SIMPLE_TR",
                    i, block_idx
                )));
            }
            Block::INSERTION { .. } => {
                return Err(ConfigError::Invariant(format!(
                    "post_generation[{}]: target block {} is an INSERTION, not HOR/SIMPLE_TR",
                    i, block_idx
                )));
            }
        };
        if copy_start < 1 || copy_end > n_copies {
            return Err(ConfigError::Invariant(format!(
                "post_generation[{}]: copy range [{}..{}] does not fit in {} block {} (which has {} copies)",
                i, copy_start, copy_end, kind, block_idx, n_copies
            )));
        }

        // HYBRID also constrains slot.
        if let Event::HYBRID {
            slot, source_slots, ..
        } = ev
        {
            let template_k = template_k_of(cfg, block_idx);
            if let Some(k) = template_k {
                if *slot < 1 || *slot > k {
                    return Err(ConfigError::Invariant(format!(
                        "post_generation[{}]: slot {} out of range for HOR block (k={})",
                        i, slot, k
                    )));
                }
                for s in source_slots {
                    if *s < 1 || *s > k {
                        return Err(ConfigError::Invariant(format!(
                            "post_generation[{}]: source_slot {} out of range (k={})",
                            i, s, k
                        )));
                    }
                }
            }
        }
    }

    // Q5: negative SHIFT bounds — must follow a HOR/SIMPLE_TR; |offset|
    // ≤ monomer_length / 2.
    for (i, b) in cfg.structure.iter().enumerate() {
        if let Block::SHIFT { offset_bp } = b {
            if *offset_bp < 0 {
                let prev_tpl = (0..i).rev().find_map(|j| match &cfg.structure[j] {
                    Block::HOR { template, .. } | Block::SIMPLE_TR { template, .. } => {
                        Some(template.clone())
                    }
                    _ => None,
                });
                let tpl_name = prev_tpl.ok_or_else(|| {
                    ConfigError::Invariant(format!(
                        "structure[{}]: negative SHIFT must follow a HOR or SIMPLE_TR block",
                        i
                    ))
                })?;
                let monomer_len = match cfg.templates.get(&tpl_name) {
                    Some(Template::HOR_slots {
                        monomer_length_bp, ..
                    })
                    | Some(Template::monomer {
                        monomer_length_bp, ..
                    }) => *monomer_length_bp,
                    None => {
                        return Err(ConfigError::Invariant(format!(
                            "structure[{}]: preceding block references unknown template `{}`",
                            i, tpl_name
                        )));
                    }
                };
                let limit = (monomer_len as i64) / 2;
                if offset_bp.abs() > limit {
                    return Err(ConfigError::Invariant(format!(
                        "structure[{}]: negative SHIFT |offset_bp|={} exceeds monomer_length/2={} of preceding block (template `{}`)",
                        i, offset_bp.abs(), limit, tpl_name
                    )));
                }
            }
        }
    }

    // Block templates must exist.
    for (i, b) in cfg.structure.iter().enumerate() {
        if let Block::HOR { template, .. } | Block::SIMPLE_TR { template, .. } = b {
            if !cfg.templates.contains_key(template) {
                return Err(ConfigError::Invariant(format!(
                    "structure[{}]: references unknown template `{}`",
                    i, template
                )));
            }
            // A HOR block must point at an HOR_slots template; SIMPLE_TR can use either.
            if let Block::HOR { template, .. } = b {
                if let Some(Template::monomer { .. }) = cfg.templates.get(template) {
                    return Err(ConfigError::Invariant(format!(
                        "structure[{}]: HOR block uses template `{}` of type `monomer` (need `HOR_slots`)",
                        i, template
                    )));
                }
            }
        }
    }

    Ok(())
}

/// Returns Some(k) if the block at index `i` is HOR and resolves to an
/// HOR_slots template; None otherwise.
fn template_k_of(cfg: &Config, i: usize) -> Option<usize> {
    let name = match cfg.structure.get(i)? {
        Block::HOR { template, .. } | Block::SIMPLE_TR { template, .. } => template,
        _ => return None,
    };
    match cfg.templates.get(name)? {
        Template::HOR_slots { k, .. } => Some(*k),
        Template::monomer { .. } => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(yaml: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        f
    }

    const MINIMAL_HOR: &str = r#"
schema_version: 1
seed: 42
templates:
  alpha_A:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha_A
    n_copies: 100
"#;

    #[test]
    fn minimal_hor_validates() {
        let f = write_tmp(MINIMAL_HOR);
        let cfg = load_and_validate(f.path()).unwrap();
        assert_eq!(cfg.schema_version, 1);
        assert_eq!(cfg.seed, 42);
        assert_eq!(cfg.structure.len(), 1);
        assert!(cfg.templates.contains_key("alpha_A"));
    }

    #[test]
    fn schema_version_must_be_1() {
        let bad = MINIMAL_HOR.replace("schema_version: 1", "schema_version: 2");
        let f = write_tmp(&bad);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("schema_version"));
    }

    #[test]
    fn unknown_field_rejected() {
        let bad = format!("{}\nnonsense_field: 99\n", MINIMAL_HOR);
        let f = write_tmp(&bad);
        assert!(load_and_validate(f.path()).is_err());
    }

    #[test]
    fn source_file_rejected_with_mvp_message() {
        let yaml = r#"
schema_version: 1
templates:
  bad:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    source: file
    file: /nonexistent.fa
structure:
  - type: HOR
    template: bad
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("source: file is not implemented in MVP"),
            "got: {msg}"
        );
    }

    #[test]
    fn event_block_oor_rejected() {
        let yaml = format!(
            "{}\npost_generation:\n  - type: INVERSION\n    block: 5\n    start_copy: 1\n    length_copies: 1\n",
            MINIMAL_HOR
        );
        let f = write_tmp(&yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("out of range"));
    }

    #[test]
    fn event_copy_range_oor_rejected() {
        let yaml = format!(
            "{}\npost_generation:\n  - type: INVERSION\n    block: 0\n    start_copy: 95\n    length_copies: 20\n",
            MINIMAL_HOR
        );
        let f = write_tmp(&yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("does not fit"));
    }

    #[test]
    fn event_targeting_shift_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 50
  - type: SHIFT
    offset_bp: 25
post_generation:
  - type: INVERSION
    block: 1
    start_copy: 1
    length_copies: 1
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("SHIFT"), "got: {err}");
    }

    #[test]
    fn hybrid_slot_oor_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 50
post_generation:
  - type: HYBRID
    block: 0
    at_copy: 10
    slot: 9
    source_slots: [1, 2]
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("slot"), "got: {err}");
    }

    #[test]
    fn negative_shift_within_bound_ok() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 10
  - type: SHIFT
    offset_bp: -40
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        load_and_validate(f.path()).unwrap();
    }

    #[test]
    fn negative_shift_too_large_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 10
  - type: SHIFT
    offset_bp: -75
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("exceeds monomer_length/2"));
    }

    #[test]
    fn negative_shift_no_preceding_block_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: SHIFT
    offset_bp: -10
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("must follow"));
    }

    #[test]
    fn hor_with_monomer_template_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: monomer
    monomer_length_bp: 100
structure:
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("HOR_slots"));
    }

    #[test]
    fn unknown_template_rejected() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: nonsuch
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        assert!(format!("{err}").contains("unknown template"));
    }

    #[test]
    fn json_schema_rejects_mutation_rate_above_1() {
        // F1: schema says `global.mutation_rate` is in [0, 1]; runtime
        // must reject values above the bound.
        let yaml = r#"
schema_version: 1
global:
  mutation_rate: 1.5
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("Schema validation failed") || msg.contains("mutation_rate"),
            "expected schema error mentioning mutation_rate; got: {msg}"
        );
    }

    #[test]
    fn json_schema_rejects_k_below_2() {
        // F1: HOR_slots requires k >= 2.
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 1
structure:
  - type: HOR
    template: t
    n_copies: 10
"#;
        let f = write_tmp(yaml);
        assert!(load_and_validate(f.path()).is_err());
    }

    #[test]
    fn json_schema_rejects_inline_sequence_with_non_acgt() {
        // F1: schema's `sequence` pattern is ^[ACGTacgt]+$.
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: monomer
    monomer_length_bp: 10
    source: sequence
    sequence: "ACGTXNNNNN"
structure:
  - type: SIMPLE_TR
    template: t
    n_copies: 5
"#;
        let f = write_tmp(yaml);
        assert!(load_and_validate(f.path()).is_err());
    }

    #[test]
    fn json_schema_rejects_empty_structure() {
        // F1: schema's structure has minItems=1.
        let yaml = r#"
schema_version: 1
templates: {}
structure: []
"#;
        let f = write_tmp(yaml);
        assert!(load_and_validate(f.path()).is_err());
    }

    #[test]
    fn inline_sequence_wrong_length_rejected() {
        // F4: sequence length must equal monomer_length_bp.
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: monomer
    monomer_length_bp: 10
    source: sequence
    sequence: "ACGTACGT"
structure:
  - type: SIMPLE_TR
    template: t
    n_copies: 5
"#;
        let f = write_tmp(yaml);
        let err = load_and_validate(f.path()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("must equal monomer_length_bp"), "got: {msg}");
    }

    #[test]
    fn inline_sequence_correct_length_accepted() {
        let yaml = r#"
schema_version: 1
templates:
  t:
    type: monomer
    monomer_length_bp: 10
    source: sequence
    sequence: "ACGTACGTAC"
structure:
  - type: SIMPLE_TR
    template: t
    n_copies: 5
"#;
        let f = write_tmp(yaml);
        load_and_validate(f.path()).unwrap();
    }

    #[test]
    fn insertion_kinds_parse() {
        for kind in ["random", "AT_rich", "GC_rich", "retro_like", "segdup_like"] {
            let yaml = format!(
                r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
structure:
  - type: HOR
    template: t
    n_copies: 10
  - type: INSERTION
    length_bp: 500
    kind: {kind}
  - type: HOR
    template: t
    n_copies: 10
"#
            );
            let f = write_tmp(&yaml);
            load_and_validate(f.path()).unwrap_or_else(|e| panic!("kind={kind} failed: {e}"));
        }
    }
}

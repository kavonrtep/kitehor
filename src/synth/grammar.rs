//! Structural-expression emission (`docs/new/taxonomy.md §2` grammar).
//!
//! Walks the `structure` list of a [`crate::synth::Config`] and emits a
//! single-line string like:
//!
//! ```text
//! H([M_1..M_12],100,div=0.15)|shift(85)|H([M_1..M_12],100,div=0.15)
//! ```
//!
//! Used in `truth.tsv::structural_expression`. Post-generation events
//! and ambient mutation/indel rates are **not** included here — they
//! live in `events_json` and the noise truth columns respectively.

use crate::synth::config::{Block, Config, InsertionKind, Template};

pub fn to_expression(cfg: &Config) -> String {
    let mut out = String::new();
    let blocks = &cfg.structure;

    // Helper: render one repeat / insertion block. Returns Some(string)
    // for HOR/SIMPLE_TR/INSERTION and None for SHIFT (handled inline by
    // the separator logic).
    let mut renders: Vec<(usize, String)> = Vec::new();
    for (i, b) in blocks.iter().enumerate() {
        if let Some(s) = render_block(cfg, b) {
            renders.push((i, s));
        }
    }

    for w in renders.windows(2).enumerate() {
        let (idx, pair) = w;
        let (li, l_str) = &pair[0];
        let (ri, _) = &pair[1];
        if idx == 0 {
            out.push_str(l_str);
        }
        // Anything between li and ri in `structure` — should be SHIFTs (or
        // possibly more INSERTIONs which we already rendered). If exactly
        // one SHIFT lies between, emit a |shift(δ)| separator; otherwise
        // fall back to '+'.
        let between = &blocks[*li + 1..*ri];
        let mut joined = false;
        if between.len() == 1 {
            if let Block::SHIFT { offset_bp } = &between[0] {
                out.push_str(&format!("|shift({})|", offset_bp));
                joined = true;
            }
        }
        if !joined {
            out.push('+');
        }
        out.push_str(&pair[1].1);
    }

    // Single-block edge case.
    if renders.len() == 1 {
        out.push_str(&renders[0].1);
    }
    if renders.is_empty() {
        // All-SHIFT structure: nothing to render. Should never happen
        // because the validator rejects this, but produce something
        // safe.
        out.push_str("(empty)");
    }

    out
}

fn render_block(cfg: &Config, b: &Block) -> Option<String> {
    match b {
        Block::HOR {
            template,
            n_copies,
        } => {
            let (k, div) = template_k_div(cfg, template);
            let div_part = if div > 0.0 {
                format!(",div={}", trim_float(div))
            } else {
                String::new()
            };
            Some(format!("H([M_1..M_{}],{}{})", k, n_copies, div_part))
        }
        Block::SIMPLE_TR {
            template,
            n_copies,
        } => {
            let l = template_monomer_len(cfg, template);
            Some(format!("T(M({}),{})", l, n_copies))
        }
        Block::SHIFT { .. } => None,
        Block::INSERTION { length_bp, kind } => {
            Some(format!("INS({},{})", length_bp, kind_str(*kind)))
        }
    }
}

fn template_k_div(cfg: &Config, name: &str) -> (usize, f64) {
    match cfg.templates.get(name) {
        Some(Template::HOR_slots {
            k,
            inter_slot_divergence,
            ..
        }) => (*k, *inter_slot_divergence),
        _ => (1, 0.0),
    }
}

fn template_monomer_len(cfg: &Config, name: &str) -> usize {
    match cfg.templates.get(name) {
        Some(Template::HOR_slots {
            monomer_length_bp, ..
        })
        | Some(Template::monomer {
            monomer_length_bp, ..
        }) => *monomer_length_bp,
        None => 0,
    }
}

fn kind_str(k: InsertionKind) -> &'static str {
    match k {
        InsertionKind::Random => "random",
        InsertionKind::AtRich => "AT_rich",
        InsertionKind::GcRich => "GC_rich",
        InsertionKind::RetroLike => "retro_like",
        InsertionKind::SegdupLike => "segdup_like",
    }
}

fn trim_float(x: f64) -> String {
    // Round to 3 decimal places, strip trailing zeros / trailing '.'.
    let s = format!("{:.3}", x);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s.is_empty() {
        "0".to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn parse(yaml: &str) -> Config {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        crate::synth::config::load_and_validate(f.path()).unwrap()
    }

    #[test]
    fn single_hor() {
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
        assert_eq!(to_expression(&cfg), "H([M_1..M_12],100,div=0.15)");
    }

    #[test]
    fn phase_shifted_hor() {
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
  - type: SHIFT
    offset_bp: 85
  - type: HOR
    template: t
    n_copies: 100
"#,
        );
        assert_eq!(
            to_expression(&cfg),
            "H([M_1..M_12],100,div=0.15)|shift(85)|H([M_1..M_12],100,div=0.15)"
        );
    }

    #[test]
    fn hor_with_retro_insertion() {
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
    n_copies: 50
  - type: INSERTION
    length_bp: 5000
    kind: retro_like
  - type: HOR
    template: t
    n_copies: 50
"#,
        );
        assert_eq!(
            to_expression(&cfg),
            "H([M_1..M_12],50,div=0.15)+INS(5000,retro_like)+H([M_1..M_12],50,div=0.15)"
        );
    }

    #[test]
    fn simple_tr() {
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
        assert_eq!(to_expression(&cfg), "T(M(170),1000)");
    }

    #[test]
    fn zero_divergence_omitted() {
        let cfg = parse(
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
    n_copies: 50
"#,
        );
        // div=0 → div part should be absent.
        assert_eq!(to_expression(&cfg), "H([M_1..M_4],50)");
    }
}

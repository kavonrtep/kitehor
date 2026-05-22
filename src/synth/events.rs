//! Post-generation events.
//!
//! Each event resolves its logical (block, copy[, slot]) target via
//! `coord_map`, mutates `state.sequence` in place, and updates
//! `coord_map` + `filler_spans` to keep the realised coordinate model
//! consistent for downstream stages (noise) and the truth/diagnostics
//! pipeline.
//!
//! Conventions:
//! - **HYBRID**: replaces one slot's bytes with a chimera composed of
//!   `source_slots[0]`'s 5′ half + `source_slots[1]`'s 3′ half (split
//!   at `split_fraction × monomer_len`, rounded). Sequence length
//!   unchanged; coord_map unchanged.
//! - **INVERSION**: reverse-complements the bytes spanning the
//!   contiguous copy range. Sequence length unchanged. coord_map
//!   positions inside the range are mirrored so future lookups still
//!   point to *where the original bytes physically live* (now
//!   reverse-complemented).
//! - **DUPLICATION**: inserts a verbatim copy of the byte range
//!   immediately after the original. Sequence length += range_len.
//!   coord_map shifts via `apply_indels`; the duplicated bytes are
//!   uncovered (they're recorded only in `events_json` /
//!   `EventLog`).
//! - **DELETION**: removes the bytes of the copy range. Sequence
//!   length -= range_len. coord_map entries inside the range
//!   collapse to length 0 via `apply_indels`.

use crate::synth::blocks::SimState;
use crate::synth::config::{Block, Config, Event};
use crate::synth::coords::{apply_indels_to_span, shift_span_after};
use crate::synth::templates::InstantiatedTemplate;
use anyhow::{anyhow, bail, Result};
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub enum EventLog {
    Hybrid {
        block: usize,
        at_bp: usize,
        copy: usize,
        slot: usize,
        source_slots: [usize; 2],
    },
    Inversion {
        block: usize,
        start_bp: usize,
        length_bp: usize,
    },
    Duplication {
        block: usize,
        start_bp: usize,
        length_bp: usize,
    },
    Deletion {
        block: usize,
        start_bp: usize,
        length_bp: usize,
    },
}

pub fn apply(
    state: &mut SimState,
    events: &[Event],
    cfg: &Config,
    templates: &HashMap<String, InstantiatedTemplate>,
    _rng: &mut ChaCha20Rng,
) -> Result<Vec<EventLog>> {
    let mut out = Vec::with_capacity(events.len());
    for (i, ev) in events.iter().enumerate() {
        let log = match ev {
            Event::HYBRID {
                block,
                at_copy,
                slot,
                source_slots,
                split_fraction,
            } => apply_hybrid(
                state,
                *block,
                *at_copy,
                *slot,
                *source_slots,
                *split_fraction,
                cfg,
                templates,
            )
            .with_context_suffix(i)?,
            Event::INVERSION {
                block,
                start_copy,
                length_copies,
            } => apply_inversion(state, *block, *start_copy, *length_copies)
                .with_context_suffix(i)?,
            Event::DUPLICATION {
                block,
                start_copy,
                length_copies,
            } => apply_duplication(state, *block, *start_copy, *length_copies)
                .with_context_suffix(i)?,
            Event::DELETION {
                block,
                start_copy,
                length_copies,
            } => {
                apply_deletion(state, *block, *start_copy, *length_copies).with_context_suffix(i)?
            }
        };
        out.push(log);
    }
    Ok(out)
}

trait ContextSuffix<T> {
    fn with_context_suffix(self, i: usize) -> Result<T>;
}
impl<T> ContextSuffix<T> for Result<T> {
    fn with_context_suffix(self, i: usize) -> Result<T> {
        self.map_err(|e| anyhow!("post_generation[{}]: {}", i, e))
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_hybrid(
    state: &mut SimState,
    block: usize,
    at_copy: usize,
    slot: usize,
    source_slots: [usize; 2],
    split_fraction: f64,
    cfg: &Config,
    templates: &HashMap<String, InstantiatedTemplate>,
) -> Result<EventLog> {
    let entry = *state
        .coord_map
        .find(block, at_copy, slot)
        .ok_or_else(|| anyhow!("no coord entry for block={block} copy={at_copy} slot={slot}"))?;
    let template_name = repeat_template_name(cfg, block)?;
    let inst = templates
        .get(&template_name)
        .ok_or_else(|| anyhow!("template `{template_name}` not instantiated"))?;
    if inst.slots.len() < 2 {
        bail!("HYBRID needs an HOR_slots template (k>=2)");
    }
    let s_a = &inst.slots[source_slots[0] - 1];
    let s_b = &inst.slots[source_slots[1] - 1];
    if s_a.len() != entry.realised_len_bp || s_b.len() != entry.realised_len_bp {
        bail!(
            "HYBRID source slot lengths {}/{} do not match entry length {}",
            s_a.len(),
            s_b.len(),
            entry.realised_len_bp
        );
    }
    let split = ((entry.realised_len_bp as f64) * split_fraction).round() as usize;
    let split = split.min(entry.realised_len_bp);
    let mut chimera = Vec::with_capacity(entry.realised_len_bp);
    chimera.extend_from_slice(&s_a[..split]);
    chimera.extend_from_slice(&s_b[split..]);
    state.sequence[entry.realised_start_bp..entry.realised_start_bp + entry.realised_len_bp]
        .copy_from_slice(&chimera);
    Ok(EventLog::Hybrid {
        block,
        at_bp: entry.realised_start_bp,
        copy: at_copy,
        slot,
        source_slots,
    })
}

fn apply_inversion(
    state: &mut SimState,
    block: usize,
    start_copy: usize,
    length_copies: usize,
) -> Result<EventLog> {
    let (start_bp, end_bp) = copy_range_bp(state, block, start_copy, length_copies)?;
    let len = end_bp - start_bp;
    let slice = &mut state.sequence[start_bp..end_bp];
    slice.reverse();
    for b in slice.iter_mut() {
        *b = complement(*b);
    }
    // Mirror coord_map entries fully contained in [start_bp, end_bp).
    for e in &mut state.coord_map.entries {
        let s = e.realised_start_bp;
        let t = s + e.realised_len_bp;
        if s >= start_bp && t <= end_bp {
            // new_start = start_bp + (end_bp - t)
            e.realised_start_bp = start_bp + (end_bp - t);
        }
    }
    Ok(EventLog::Inversion {
        block,
        start_bp,
        length_bp: len,
    })
}

fn apply_duplication(
    state: &mut SimState,
    block: usize,
    start_copy: usize,
    length_copies: usize,
) -> Result<EventLog> {
    let (start_bp, end_bp) = copy_range_bp(state, block, start_copy, length_copies)?;
    let len = end_bp - start_bp;
    let dup: Vec<u8> = state.sequence[start_bp..end_bp].to_vec();
    state.sequence.splice(end_bp..end_bp, dup);
    // F2: the duplicated bytes are an *uncovered* structural filler.
    // Entries at or past `end_bp` must shift right by `len`; entries
    // wholly inside `[start_bp, end_bp)` (the original copies) are
    // untouched. Using `apply_indels` here would incorrectly absorb
    // the duplicate into whichever entry starts at `end_bp`.
    state.coord_map.shift_after(end_bp, len as i64);
    shift_fillers_after(state, end_bp, len as i64);
    Ok(EventLog::Duplication {
        block,
        start_bp,
        length_bp: len,
    })
}

fn apply_deletion(
    state: &mut SimState,
    block: usize,
    start_copy: usize,
    length_copies: usize,
) -> Result<EventLog> {
    let (start_bp, end_bp) = copy_range_bp(state, block, start_copy, length_copies)?;
    let len = end_bp - start_bp;
    state.sequence.drain(start_bp..end_bp);
    let indels: Vec<(usize, i32)> = (start_bp..end_bp).map(|p| (p, -1)).collect();
    state.coord_map.apply_indels(&indels);
    apply_to_fillers(state, &indels);
    Ok(EventLog::Deletion {
        block,
        start_bp,
        length_bp: len,
    })
}

fn apply_to_fillers(state: &mut SimState, indels: &[(usize, i32)]) {
    for fs in state.filler_spans.iter_mut() {
        let (s, l) = apply_indels_to_span(fs.realised_start_bp, fs.realised_len_bp, indels);
        fs.realised_start_bp = s;
        fs.realised_len_bp = l;
    }
}

fn shift_fillers_after(state: &mut SimState, pos: usize, delta: i64) {
    for fs in state.filler_spans.iter_mut() {
        let (s, l) = shift_span_after(fs.realised_start_bp, fs.realised_len_bp, pos, delta);
        fs.realised_start_bp = s;
        fs.realised_len_bp = l;
    }
}

/// Compute the contiguous bp span of copies `[start_copy, start_copy +
/// length_copies)` in `block`. Validates copies are present and
/// contiguous.
fn copy_range_bp(
    state: &SimState,
    block: usize,
    start_copy: usize,
    length_copies: usize,
) -> Result<(usize, usize)> {
    let end_copy_excl = start_copy + length_copies;
    // Find min start and max end among entries matching the criteria.
    let mut min_s = usize::MAX;
    let mut max_e = 0usize;
    let mut count = 0;
    for e in &state.coord_map.entries {
        if e.block_idx == block && e.copy_idx >= start_copy && e.copy_idx < end_copy_excl {
            min_s = min_s.min(e.realised_start_bp);
            max_e = max_e.max(e.realised_start_bp + e.realised_len_bp);
            count += 1;
        }
    }
    if count == 0 {
        bail!(
            "no coord entries match block={block} copies={start_copy}..{} ",
            end_copy_excl - 1
        );
    }
    Ok((min_s, max_e))
}

fn repeat_template_name(cfg: &Config, block: usize) -> Result<String> {
    match cfg.structure.get(block) {
        Some(Block::HOR { template, .. }) | Some(Block::SIMPLE_TR { template, .. }) => {
            Ok(template.clone())
        }
        _ => bail!("block {block} is not HOR/SIMPLE_TR"),
    }
}

#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        x => x,
    }
}

pub fn to_events_json(logs: &[EventLog]) -> String {
    if logs.is_empty() {
        return "[]".to_string();
    }
    let parts: Vec<String> = logs
        .iter()
        .map(|e| match e {
            EventLog::Hybrid {
                block,
                at_bp,
                copy,
                slot,
                source_slots,
            } => format!(
                r#"{{"type":"HYBRID","block":{block},"at_bp":{at_bp},"copy":{copy},"slot":{slot},"source_slots":[{},{}]}}"#,
                source_slots[0], source_slots[1]
            ),
            EventLog::Inversion {
                block,
                start_bp,
                length_bp,
            } => format!(
                r#"{{"type":"INVERSION","block":{block},"start_bp":{start_bp},"length_bp":{length_bp}}}"#
            ),
            EventLog::Duplication {
                block,
                start_bp,
                length_bp,
            } => format!(
                r#"{{"type":"DUPLICATION","block":{block},"start_bp":{start_bp},"length_bp":{length_bp}}}"#
            ),
            EventLog::Deletion {
                block,
                start_bp,
                length_bp,
            } => format!(
                r#"{{"type":"DELETION","block":{block},"start_bp":{start_bp},"length_bp":{length_bp}}}"#
            ),
        })
        .collect();
    format!("[{}]", parts.join(","))
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

    fn build_state(cfg: &Config) -> (SimState, HashMap<String, InstantiatedTemplate>) {
        let s = Streams::new(cfg.seed);
        let mut rt = s.templates();
        let inst = instantiate(&cfg.templates, &mut rt);
        let mut rs = s.structure();
        let state = expand(&cfg.structure, &inst, &mut rs).unwrap();
        (state, inst)
    }

    #[test]
    fn rc_palindrome() {
        // ACGT reverse-complemented is ACGT.
        let mut bs = b"ACGT".to_vec();
        bs.reverse();
        for b in bs.iter_mut() {
            *b = complement(*b);
        }
        assert_eq!(bs, b"ACGT".to_vec());
    }

    #[test]
    fn rc_aaaa_to_tttt() {
        let mut bs = b"AAAA".to_vec();
        bs.reverse();
        for b in bs.iter_mut() {
            *b = complement(*b);
        }
        assert_eq!(bs, b"TTTT".to_vec());
    }

    #[test]
    fn hybrid_replaces_slot_bytes() {
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.20
structure:
  - type: HOR
    template: t
    n_copies: 50
post_generation:
  - type: HYBRID
    block: 0
    at_copy: 27
    slot: 3
    source_slots: [3, 4]
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        let entry_before = *state.coord_map.find(0, 27, 3).unwrap();
        let pre_len = state.sequence.len();
        let mut rng = Streams::new(1).events();
        let logs = apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        assert_eq!(state.sequence.len(), pre_len, "HYBRID must preserve length");
        // The bytes in the slot should now be the chimera: first 50 of
        // slot 3, next 50 of slot 4.
        let s_3 = &inst["t"].slots[2];
        let s_4 = &inst["t"].slots[3];
        let actual =
            &state.sequence[entry_before.realised_start_bp..entry_before.realised_start_bp + 100];
        assert_eq!(&actual[..50], &s_3[..50]);
        assert_eq!(&actual[50..], &s_4[50..]);
        // Event log
        match &logs[0] {
            EventLog::Hybrid {
                block,
                at_bp,
                copy,
                slot,
                source_slots,
            } => {
                assert_eq!(*block, 0);
                assert_eq!(*at_bp, entry_before.realised_start_bp);
                assert_eq!(*copy, 27);
                assert_eq!(*slot, 3);
                assert_eq!(*source_slots, [3, 4]);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn inversion_preserves_length() {
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
post_generation:
  - type: INVERSION
    block: 0
    start_copy: 11
    length_copies: 10
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        let pre_len = state.sequence.len();
        let mut rng = Streams::new(1).events();
        let logs = apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        assert_eq!(state.sequence.len(), pre_len);
        match &logs[0] {
            EventLog::Inversion {
                block,
                start_bp,
                length_bp,
            } => {
                assert_eq!(*block, 0);
                // copy 11..21, each 4 slots × 100 bp = 4000 bp.
                assert_eq!(*length_bp, 4000);
                assert_eq!(*start_bp, 10 * 4 * 100);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn inversion_actually_reverse_complements_bytes() {
        // Use a deterministic template so we can read bytes back.
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
    n_copies: 5
post_generation:
  - type: INVERSION
    block: 0
    start_copy: 2
    length_copies: 2
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        // Snapshot the bytes that will be inverted.
        let pre_seq = state.sequence.clone();
        let start_bp = 4 * 100;
        let end_bp = 3 * 4 * 100;
        let pre_range: Vec<u8> = pre_seq[start_bp..end_bp].to_vec();
        let mut rng = Streams::new(1).events();
        apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        let post_range: Vec<u8> = state.sequence[start_bp..end_bp].to_vec();
        // Reverse-complement of pre_range should equal post_range.
        let mut expected: Vec<u8> = pre_range.iter().rev().map(|b| complement(*b)).collect();
        assert_eq!(post_range, expected);
        expected.clear();
    }

    #[test]
    fn duplication_extends_sequence_by_range_len() {
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
post_generation:
  - type: DUPLICATION
    block: 0
    start_copy: 20
    length_copies: 5
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        let pre_len = state.sequence.len();
        let mut rng = Streams::new(1).events();
        let logs = apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        // 5 copies × 4 slots × 100 bp = 2000 bp added.
        assert_eq!(state.sequence.len(), pre_len + 2000);
        match &logs[0] {
            EventLog::Duplication { length_bp, .. } => assert_eq!(*length_bp, 2000),
            _ => panic!(),
        }
    }

    #[test]
    fn deletion_shrinks_sequence_by_range_len() {
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
post_generation:
  - type: DELETION
    block: 0
    start_copy: 20
    length_copies: 5
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        let pre_len = state.sequence.len();
        let mut rng = Streams::new(1).events();
        let logs = apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        assert_eq!(state.sequence.len(), pre_len - 2000);
        match &logs[0] {
            EventLog::Deletion { length_bp, .. } => assert_eq!(*length_bp, 2000),
            _ => panic!(),
        }
    }

    #[test]
    fn events_json_round_trip() {
        let logs = vec![
            EventLog::Hybrid {
                block: 0,
                at_bp: 1234,
                copy: 27,
                slot: 4,
                source_slots: [4, 5],
            },
            EventLog::Inversion {
                block: 1,
                start_bp: 10_000,
                length_bp: 4000,
            },
            EventLog::Duplication {
                block: 0,
                start_bp: 50_000,
                length_bp: 2000,
            },
            EventLog::Deletion {
                block: 2,
                start_bp: 60_000,
                length_bp: 1000,
            },
        ];
        let j = to_events_json(&logs);
        let v: serde_json::Value = serde_json::from_str(&j).expect("emit valid JSON");
        assert_eq!(v.as_array().unwrap().len(), 4);
        assert_eq!(v[0]["type"], "HYBRID");
        assert_eq!(v[0]["block"], 0);
        assert_eq!(v[0]["copy"], 27);
        assert_eq!(v[1]["type"], "INVERSION");
        assert_eq!(v[1]["block"], 1);
        assert_eq!(v[2]["type"], "DUPLICATION");
        assert_eq!(v[2]["block"], 0);
        assert_eq!(v[3]["type"], "DELETION");
        assert_eq!(v[3]["block"], 2);
    }

    #[test]
    fn duplication_downstream_entries_shift_not_extend() {
        // F2 regression: reviewer's explicit ask — duplicate copies 20..24,
        // assert copy 25 starts at old_end_bp + duplicated_len.
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 50
post_generation:
  - type: DUPLICATION
    block: 0
    start_copy: 20
    length_copies: 5
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        // Capture pre-event coords.
        let copy25_pre = *state.coord_map.find(0, 25, 1).unwrap();
        let dup_end_pre = state.coord_map.find(0, 24, 4).unwrap().end_bp(); // == copy25_pre.start
        assert_eq!(copy25_pre.realised_start_bp, dup_end_pre);

        let mut rng = Streams::new(1).events();
        apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();

        // 5 copies × 4 slots × 100 bp = 2000 bp duplicated.
        let dup_len = 5 * 4 * 100;
        let copy25_post = state.coord_map.find(0, 25, 1).unwrap();
        assert_eq!(
            copy25_post.realised_start_bp,
            dup_end_pre + dup_len,
            "copy 25 must shift past the duplicated region, not absorb it"
        );
        assert_eq!(
            copy25_post.realised_len_bp, 100,
            "copy 25 length must NOT grow"
        );
        // And copy 20 (first of the duplicated range) stays at its
        // original position.
        let copy20_post = state.coord_map.find(0, 20, 1).unwrap();
        assert_eq!(copy20_post.realised_start_bp, 19 * 4 * 100);
    }

    #[test]
    fn event_chain_dup_then_inv_then_del_keeps_coords_consistent() {
        // Reviewer's "additional improvements" ask: chain events with
        // coordinates downstream of each prior event and verify the
        // running coord_map still resolves logical (block, copy, slot)
        // to valid byte ranges.
        let cfg = parse(
            r#"
schema_version: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 100
post_generation:
  - type: DUPLICATION
    block: 0
    start_copy: 20
    length_copies: 5
  - type: INVERSION
    block: 0
    start_copy: 40
    length_copies: 3
  - type: DELETION
    block: 0
    start_copy: 80
    length_copies: 2
"#,
        );
        let (mut state, inst) = build_state(&cfg);
        let pre_len = state.sequence.len();
        let mut rng = Streams::new(1).events();
        let logs = apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut rng).unwrap();
        // Final length = pre + dup - del = pre + 2000 - 800.
        assert_eq!(state.sequence.len(), pre_len + 2000 - 800);
        assert_eq!(logs.len(), 3);
        // Every remaining (block, copy, slot) lookup must yield a span
        // that lies entirely within the post-event sequence.
        for entry in &state.coord_map.entries {
            let e = entry.realised_start_bp + entry.realised_len_bp;
            assert!(
                e <= state.sequence.len(),
                "coord entry {:?} extends past sequence (len={})",
                entry,
                state.sequence.len()
            );
        }
    }

    #[test]
    fn empty_events_json_is_bracket_pair() {
        assert_eq!(to_events_json(&[]), "[]");
    }
}

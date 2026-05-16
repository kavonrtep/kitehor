//! Block expansion.
//!
//! Walks the YAML `structure` list, appends bytes to the output
//! sequence, and records one `CoordEntry` per emitted slot of every
//! `HOR` or `SIMPLE_TR` block. `SHIFT` and `INSERTION` blocks emit
//! bytes that are intentionally **not** indexed in the coord_map —
//! they're structural fillers, not slot bytes. The truth/periods
//! pipeline (M4) walks the original `structure` list separately to
//! emit shift coordinates and period candidates.
//!
//! Realised positions of SHIFT and INSERTION blocks are remembered in
//! `SimState::filler_spans` so the truth file can report
//! `phase_shift_positions` and `events_json` insertion start positions
//! against the **post-expansion** sequence (subsequent passes —
//! wobble, events, noise — will adjust these via `apply_indels`).

use crate::synth::config::{Block, InsertionKind};
use crate::synth::coords::{CoordEntry, CoordMap};
use crate::synth::templates::InstantiatedTemplate;
use anyhow::{anyhow, bail, Result};
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;

/// Span emitted by a non-repeat block, so the truth/diagnostics layer
/// can report it without rescanning the sequence. Realised positions
/// are post-expansion; subsequent stages (wobble, events, noise)
/// shift them via `CoordMap::apply_indels`.
#[derive(Debug, Clone, Copy)]
pub struct FillerSpan {
    pub block_idx: usize,
    pub kind: FillerKind,
    pub realised_start_bp: usize,
    pub realised_len_bp: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillerKind {
    Shift { offset_bp: i64 },
    Insertion(InsertionKind),
}

#[derive(Debug, Default)]
pub struct SimState {
    pub sequence: Vec<u8>,
    pub coord_map: CoordMap,
    pub filler_spans: Vec<FillerSpan>,
}

pub fn expand(
    structure: &[Block],
    templates: &HashMap<String, InstantiatedTemplate>,
    rng: &mut ChaCha20Rng,
) -> Result<SimState> {
    let mut st = SimState::default();
    for (bi, block) in structure.iter().enumerate() {
        match block {
            Block::HOR {
                template,
                n_copies,
            } => emit_hor(bi, template, *n_copies, templates, &mut st)?,
            Block::SIMPLE_TR {
                template,
                n_copies,
            } => emit_simple_tr(bi, template, *n_copies, templates, &mut st)?,
            Block::SHIFT { offset_bp } => emit_shift(bi, *offset_bp, &mut st, rng)?,
            Block::INSERTION { length_bp, kind } => {
                emit_insertion(bi, *length_bp, *kind, &mut st, rng)?
            }
        }
    }
    Ok(st)
}

fn emit_hor(
    bi: usize,
    template: &str,
    n_copies: usize,
    templates: &HashMap<String, InstantiatedTemplate>,
    st: &mut SimState,
) -> Result<()> {
    let inst = templates
        .get(template)
        .ok_or_else(|| anyhow!("block {bi}: template `{template}` not instantiated"))?;
    let k = inst.slots.len();
    if k < 2 {
        bail!("block {bi}: HOR template `{template}` has k={k}; HOR requires k>=2");
    }
    for copy_idx in 1..=n_copies {
        for slot_idx in 1..=k {
            let slot = &inst.slots[slot_idx - 1];
            let start = st.sequence.len();
            st.sequence.extend_from_slice(slot);
            st.coord_map.push(CoordEntry {
                block_idx: bi,
                copy_idx,
                slot_idx,
                realised_start_bp: start,
                realised_len_bp: slot.len(),
            });
        }
    }
    Ok(())
}

fn emit_simple_tr(
    bi: usize,
    template: &str,
    n_copies: usize,
    templates: &HashMap<String, InstantiatedTemplate>,
    st: &mut SimState,
) -> Result<()> {
    let inst = templates
        .get(template)
        .ok_or_else(|| anyhow!("block {bi}: template `{template}` not instantiated"))?;
    let slot = &inst.slots[0];
    for copy_idx in 1..=n_copies {
        let start = st.sequence.len();
        st.sequence.extend_from_slice(slot);
        st.coord_map.push(CoordEntry {
            block_idx: bi,
            copy_idx,
            slot_idx: 1,
            realised_start_bp: start,
            realised_len_bp: slot.len(),
        });
    }
    Ok(())
}

fn emit_shift(
    bi: usize,
    offset_bp: i64,
    st: &mut SimState,
    rng: &mut ChaCha20Rng,
) -> Result<()> {
    if offset_bp > 0 {
        let n = offset_bp as usize;
        let start = st.sequence.len();
        let bytes = draw_local_composition(&st.sequence, n, rng);
        st.sequence.extend_from_slice(&bytes);
        st.filler_spans.push(FillerSpan {
            block_idx: bi,
            kind: FillerKind::Shift { offset_bp },
            realised_start_bp: start,
            realised_len_bp: n,
        });
    } else if offset_bp < 0 {
        let n = offset_bp.unsigned_abs() as usize;
        if n > st.sequence.len() {
            bail!(
                "block {bi}: negative SHIFT |{}| > current sequence length {}",
                n,
                st.sequence.len()
            );
        }
        let new_len = st.sequence.len() - n;
        // Update coord_map: any entry whose end extends past new_len
        // shrinks; any entry entirely past new_len has its length zeroed
        // (validator should already preclude this).
        for e in &mut st.coord_map.entries {
            let end = e.realised_start_bp + e.realised_len_bp;
            if e.realised_start_bp >= new_len {
                e.realised_len_bp = 0;
            } else if end > new_len {
                e.realised_len_bp = new_len - e.realised_start_bp;
            }
        }
        st.sequence.truncate(new_len);
        st.filler_spans.push(FillerSpan {
            block_idx: bi,
            kind: FillerKind::Shift { offset_bp },
            realised_start_bp: new_len,
            realised_len_bp: 0,
        });
    }
    Ok(())
}

fn emit_insertion(
    bi: usize,
    length_bp: usize,
    kind: InsertionKind,
    st: &mut SimState,
    rng: &mut ChaCha20Rng,
) -> Result<()> {
    let start = st.sequence.len();
    let bytes = match kind {
        InsertionKind::Random => random_dna_bernoulli(length_bp, 0.5, rng),
        InsertionKind::AtRich => random_dna_bernoulli(length_bp, 0.2, rng),
        InsertionKind::GcRich => random_dna_bernoulli(length_bp, 0.8, rng),
        InsertionKind::RetroLike => retro_like(length_bp, rng),
        InsertionKind::SegdupLike => segdup_like(length_bp, &st.sequence, rng),
    };
    st.sequence.extend_from_slice(&bytes);
    st.filler_spans.push(FillerSpan {
        block_idx: bi,
        kind: FillerKind::Insertion(kind),
        realised_start_bp: start,
        realised_len_bp: length_bp,
    });
    Ok(())
}

fn draw_local_composition(out: &[u8], n: usize, rng: &mut ChaCha20Rng) -> Vec<u8> {
    let window: usize = 50;
    let source: &[u8] = if out.len() >= window {
        &out[out.len() - window..]
    } else if !out.is_empty() {
        out
    } else {
        // Empty prefix → uniform random ACGT
        return (0..n)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
    };
    (0..n)
        .map(|_| source[rng.random_range(0..source.len())])
        .collect()
}

fn random_dna_bernoulli(n: usize, gc: f64, rng: &mut ChaCha20Rng) -> Vec<u8> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        let is_gc = rng.random::<f64>() < gc;
        let coin = rng.random::<f64>() < 0.5;
        v.push(match (is_gc, coin) {
            (true, true) => b'G',
            (true, false) => b'C',
            (false, true) => b'A',
            (false, false) => b'T',
        });
    }
    v
}

/// Synthetic LTR-internal-LTR. Two identical ~200 bp LTRs flank a
/// random internal region. If `length_bp < 400`, degenerate to a flat
/// random insertion.
fn retro_like(length_bp: usize, rng: &mut ChaCha20Rng) -> Vec<u8> {
    if length_bp < 400 {
        return random_dna_bernoulli(length_bp, 0.5, rng);
    }
    let ltr_len = 200.min(length_bp / 3);
    let ltr = random_dna_bernoulli(ltr_len, 0.5, rng);
    let internal_len = length_bp - 2 * ltr_len;
    let internal = random_dna_bernoulli(internal_len, 0.5, rng);
    let mut v = Vec::with_capacity(length_bp);
    v.extend_from_slice(&ltr);
    v.extend_from_slice(&internal);
    v.extend_from_slice(&ltr);
    v
}

/// Copy a random subsequence of the requested length from the existing
/// output. Falls back to random if no source range of that length is
/// available.
fn segdup_like(length_bp: usize, src: &[u8], rng: &mut ChaCha20Rng) -> Vec<u8> {
    if src.len() < length_bp {
        return random_dna_bernoulli(length_bp, 0.5, rng);
    }
    let max_start = src.len() - length_bp;
    let start = rng.random_range(0..=max_start);
    src[start..start + length_bp].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::config::{Source, Template};
    use crate::synth::rng::Streams;
    use crate::synth::templates::instantiate;

    fn templates_with(name: &str, tpl: Template) -> HashMap<String, Template> {
        let mut h = HashMap::new();
        h.insert(name.to_string(), tpl);
        h
    }

    fn hor_slots(monomer_length_bp: usize, k: usize, d: f64) -> Template {
        Template::HOR_slots {
            monomer_length_bp,
            k,
            source: Source::Random,
            sequence: None,
            file: None,
            gc_content: 0.5,
            inter_slot_divergence: d,
        }
    }

    fn monomer_tpl(monomer_length_bp: usize) -> Template {
        Template::monomer {
            monomer_length_bp,
            source: Source::Random,
            sequence: None,
            file: None,
            gc_content: 0.5,
        }
    }

    #[test]
    fn hor_block_emits_n_copies_times_k_bases() {
        let cfgs = templates_with("t", hor_slots(100, 4, 0.10));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);

        let blocks = vec![Block::HOR {
            template: "t".into(),
            n_copies: 10,
        }];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        assert_eq!(st.sequence.len(), 100 * 4 * 10);
        assert_eq!(st.coord_map.len(), 4 * 10);
    }

    #[test]
    fn simple_tr_block_emits_n_copies_times_monomer_len() {
        let cfgs = templates_with("m", monomer_tpl(171));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);

        let blocks = vec![Block::SIMPLE_TR {
            template: "m".into(),
            n_copies: 20,
        }];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        assert_eq!(st.sequence.len(), 171 * 20);
        assert_eq!(st.coord_map.len(), 20);
    }

    #[test]
    fn coord_map_lookup_matches_realised_positions() {
        let cfgs = templates_with("t", hor_slots(100, 4, 0.10));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);

        let blocks = vec![Block::HOR {
            template: "t".into(),
            n_copies: 5,
        }];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        // copy 3, slot 2 → realised_start = 2 copies × 4 slots × 100 + 1 slot × 100 = 900
        let e = st.coord_map.find(0, 3, 2).unwrap();
        assert_eq!(e.realised_start_bp, 900);
        assert_eq!(e.realised_len_bp, 100);
        // and the bytes at that position match slot 2.
        let actual = &st.sequence[e.realised_start_bp..e.end_bp()];
        assert_eq!(actual, inst["t"].slots[1].as_slice());
    }

    #[test]
    fn two_hor_blocks_concatenate() {
        let mut cfgs = HashMap::new();
        cfgs.insert("a".into(), hor_slots(50, 3, 0.10));
        cfgs.insert("b".into(), hor_slots(40, 5, 0.15));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);

        let blocks = vec![
            Block::HOR {
                template: "a".into(),
                n_copies: 10,
            },
            Block::HOR {
                template: "b".into(),
                n_copies: 7,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        assert_eq!(st.sequence.len(), 50 * 3 * 10 + 40 * 5 * 7);
        assert_eq!(st.coord_map.len(), 3 * 10 + 5 * 7);
        // The first entry of block 1 should start where block 0 ends.
        let b0_last = st.coord_map.find(0, 10, 3).unwrap();
        let b1_first = st.coord_map.find(1, 1, 1).unwrap();
        assert_eq!(b1_first.realised_start_bp, b0_last.end_bp());
    }

    #[test]
    fn same_template_referenced_twice_shares_slots() {
        let cfgs = templates_with("alpha", hor_slots(80, 6, 0.15));
        let s = Streams::new(42);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);

        let blocks = vec![
            Block::HOR {
                template: "alpha".into(),
                n_copies: 5,
            },
            Block::HOR {
                template: "alpha".into(),
                n_copies: 5,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        // Slot 3 of copy 1 in each block should be the byte-identical
        // slot consensus from the cached InstantiatedTemplate.
        let e0 = st.coord_map.find(0, 1, 3).unwrap();
        let e1 = st.coord_map.find(1, 1, 3).unwrap();
        let b0 = &st.sequence[e0.realised_start_bp..e0.end_bp()];
        let b1 = &st.sequence[e1.realised_start_bp..e1.end_bp()];
        assert_eq!(b0, b1, "shared template must produce byte-identical slots");
    }

    #[test]
    fn positive_shift_appends_bytes_with_no_coord_entry() {
        let cfgs = templates_with("t", hor_slots(100, 4, 0.1));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![
            Block::HOR {
                template: "t".into(),
                n_copies: 10,
            },
            Block::SHIFT { offset_bp: 85 },
            Block::HOR {
                template: "t".into(),
                n_copies: 10,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        let cm_total: usize = st.coord_map.entries.iter().map(|e| e.realised_len_bp).sum();
        // 2 HOR blocks × 10 copies × 4 slots × 100 bp = 8000; plus 85 bp shift = 8085.
        assert_eq!(st.sequence.len(), 8085);
        assert_eq!(cm_total, 8000);
        // The second HOR block must start at position 4100 (4000 of first block + 100 of shift… wait, 85).
        let b2_first = st.coord_map.find(2, 1, 1).unwrap();
        assert_eq!(b2_first.realised_start_bp, 4000 + 85);
        // FillerSpan present
        assert_eq!(st.filler_spans.len(), 1);
        assert_eq!(st.filler_spans[0].realised_start_bp, 4000);
        assert_eq!(st.filler_spans[0].realised_len_bp, 85);
    }

    #[test]
    fn negative_shift_truncates_preceding_block() {
        let cfgs = templates_with("t", hor_slots(100, 4, 0.1));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![
            Block::HOR {
                template: "t".into(),
                n_copies: 10,
            },
            Block::SHIFT { offset_bp: -30 },
            Block::HOR {
                template: "t".into(),
                n_copies: 10,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        // 4000 - 30 + 4000 = 7970 total bp.
        assert_eq!(st.sequence.len(), 7970);
        // The last HOR copy of block 0 (copy 10, slot 4) should be shorter by 30.
        let last = st.coord_map.find(0, 10, 4).unwrap();
        assert_eq!(last.realised_len_bp, 70);
        // Block 2's first slot starts where block 0 ended (after truncation).
        let b2_first = st.coord_map.find(2, 1, 1).unwrap();
        assert_eq!(b2_first.realised_start_bp, 3970);
    }

    #[test]
    fn insertion_block_appends_bytes_with_no_coord_entry() {
        use crate::synth::config::InsertionKind;
        let cfgs = templates_with("t", hor_slots(100, 4, 0.1));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![
            Block::HOR {
                template: "t".into(),
                n_copies: 50,
            },
            Block::INSERTION {
                length_bp: 5000,
                kind: InsertionKind::RetroLike,
            },
            Block::HOR {
                template: "t".into(),
                n_copies: 50,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        assert_eq!(st.sequence.len(), 20_000 + 5_000 + 20_000);
        let cm_total: usize = st.coord_map.entries.iter().map(|e| e.realised_len_bp).sum();
        assert_eq!(cm_total, 40_000);
        // Second HOR block starts at 25_000.
        let b2_first = st.coord_map.find(2, 1, 1).unwrap();
        assert_eq!(b2_first.realised_start_bp, 25_000);
        // retro_like: first 200 bytes should equal the last 200 bytes of the insertion.
        let ins_start = 20_000;
        let ins_end = 25_000;
        let head = &st.sequence[ins_start..ins_start + 200];
        let tail = &st.sequence[ins_end - 200..ins_end];
        assert_eq!(head, tail, "retro_like LTRs must match");
    }

    #[test]
    fn insertion_gc_skew_in_range() {
        use crate::synth::config::InsertionKind;
        let cfgs = templates_with("t", hor_slots(100, 4, 0.1));
        let s = Streams::new(1);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        // GC-rich block large enough for a useful sample.
        let blocks = vec![
            Block::HOR {
                template: "t".into(),
                n_copies: 1,
            },
            Block::INSERTION {
                length_bp: 5000,
                kind: InsertionKind::GcRich,
            },
        ];
        let mut rs = s.structure();
        let st = expand(&blocks, &inst, &mut rs).unwrap();
        // Insertion lives at positions 400..5400.
        let ins = &st.sequence[400..400 + 5000];
        let gc = ins.iter().filter(|b| matches!(**b, b'G' | b'C')).count() as f64 / 5000.0;
        assert!(gc > 0.70, "expected GC-rich >0.7, got {gc}");
    }

    #[test]
    fn hor_block_with_k1_template_rejected() {
        // HOR_slots with k=1 is rejected at the schema level (minimum: 2),
        // but emit_hor must defend against constructed inputs too.
        // Note: serde would never produce this from YAML; we build by hand.
        let mut cfgs = HashMap::new();
        cfgs.insert("bad".into(), monomer_tpl(50));
        let s = Streams::new(0);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![Block::HOR {
            template: "bad".into(),
            n_copies: 1,
        }];
        let mut rs = s.structure();
        let err = expand(&blocks, &inst, &mut rs).unwrap_err();
        assert!(format!("{err}").contains("k=1"));
    }
}

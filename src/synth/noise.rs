//! Final noise pass: per-base substitutions + indels.
//!
//! Walks `state.sequence` once. For each base, with probability
//! `mutation_rate` substitute with a different base; with probability
//! `indel_rate / 2` insert a random base immediately before;
//! with probability `indel_rate / 2` delete this base.
//!
//! Returns a [`NoiseLog`] carrying counts and the indel position list
//! (pre-noise coordinates) so the truth/diagnostics pipeline can
//! record what happened and `coord_map` can shift to match the new
//! sequence length.

use crate::synth::blocks::SimState;
use crate::synth::config::Global;
use crate::synth::coords::apply_indels_to_span;
use rand::Rng;
use rand_chacha::ChaCha20Rng;

#[derive(Debug, Clone, Default)]
pub struct NoiseLog {
    pub n_substitutions: usize,
    pub n_insertions: usize,
    pub n_deletions: usize,
    /// `(position_in_pre_noise_sequence, delta)`. `delta` is +1 for an
    /// insertion at that position or -1 for a deletion of that
    /// position. Positions are sorted ascending (we visit in order).
    pub indels: Vec<(usize, i32)>,
}

pub fn apply(state: &mut SimState, global: &Global, rng: &mut ChaCha20Rng) -> NoiseLog {
    let mut log = NoiseLog::default();
    if global.mutation_rate <= 0.0 && global.indel_rate <= 0.0 {
        return log;
    }
    let mut new_seq = Vec::with_capacity(state.sequence.len());
    let mut_rate = global.mutation_rate;
    let ins_rate = global.indel_rate * 0.5;
    let del_rate = global.indel_rate * 0.5;

    for i in 0..state.sequence.len() {
        let mut b = state.sequence[i];
        if mut_rate > 0.0 && rng.random::<f64>() < mut_rate {
            let alts: &[u8] = match b {
                b'A' => b"CGT",
                b'C' => b"AGT",
                b'G' => b"ACT",
                b'T' => b"ACG",
                _ => b"ACGT",
            };
            b = alts[rng.random_range(0..alts.len())];
            log.n_substitutions += 1;
        }
        if ins_rate > 0.0 && rng.random::<f64>() < ins_rate {
            new_seq.push(random_base(rng));
            log.indels.push((i, 1));
            log.n_insertions += 1;
        }
        if del_rate > 0.0 && rng.random::<f64>() < del_rate {
            log.indels.push((i, -1));
            log.n_deletions += 1;
            continue;
        }
        new_seq.push(b);
    }
    state.sequence = new_seq;
    state.coord_map.apply_indels(&log.indels);
    // F3: also update filler_spans so phase_shift_positions and
    // insertion spans reported in truth/diagnostics reflect the
    // post-noise FASTA, not the pre-noise sequence.
    for fs in state.filler_spans.iter_mut() {
        let (s, l) =
            apply_indels_to_span(fs.realised_start_bp, fs.realised_len_bp, &log.indels);
        fs.realised_start_bp = s;
        fs.realised_len_bp = l;
    }
    log
}

fn random_base(rng: &mut ChaCha20Rng) -> u8 {
    b"ACGT"[rng.random_range(0..4)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::config::{Source, Template};
    use crate::synth::rng::Streams;

    fn flat_seq(byte: u8, len: usize) -> SimState {
        SimState {
            sequence: vec![byte; len],
            coord_map: Default::default(),
            filler_spans: Vec::new(),
        }
    }

    #[test]
    fn zero_rates_no_op() {
        let mut s = flat_seq(b'A', 1000);
        let g = Global::default();
        let mut r = Streams::new(0).noise();
        let log = apply(&mut s, &g, &mut r);
        assert_eq!(log.n_substitutions, 0);
        assert_eq!(log.n_insertions, 0);
        assert_eq!(log.n_deletions, 0);
        assert_eq!(s.sequence, vec![b'A'; 1000]);
    }

    #[test]
    fn deterministic_same_seed() {
        let mut a = flat_seq(b'A', 5000);
        let mut b = flat_seq(b'A', 5000);
        let g = Global {
            mutation_rate: 0.05,
            indel_rate: 0.02,
            ..Default::default()
        };
        let mut r1 = Streams::new(123).noise();
        let mut r2 = Streams::new(123).noise();
        let _ = apply(&mut a, &g, &mut r1);
        let _ = apply(&mut b, &g, &mut r2);
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn mutation_rate_within_tolerance() {
        let n = 200_000;
        let mut s = flat_seq(b'A', n);
        let g = Global {
            mutation_rate: 0.05,
            indel_rate: 0.0,
            ..Default::default()
        };
        let mut r = Streams::new(7).noise();
        let log = apply(&mut s, &g, &mut r);
        let exp = n as f64 * 0.05;
        let rel_err = (log.n_substitutions as f64 - exp).abs() / exp;
        assert!(
            rel_err < 0.05,
            "expected ~{exp} subs (±5%), got {} (rel err {:.3})",
            log.n_substitutions,
            rel_err
        );
        assert_eq!(s.sequence.len(), n);
    }

    #[test]
    fn indel_rate_within_tolerance() {
        let n = 200_000;
        let mut s = flat_seq(b'A', n);
        let g = Global {
            mutation_rate: 0.0,
            indel_rate: 0.02, // 1% insertions + 1% deletions
            ..Default::default()
        };
        let mut r = Streams::new(11).noise();
        let log = apply(&mut s, &g, &mut r);
        let exp = n as f64 * 0.01;
        let ins_err = (log.n_insertions as f64 - exp).abs() / exp;
        let del_err = (log.n_deletions as f64 - exp).abs() / exp;
        assert!(ins_err < 0.05, "insertions: exp ~{exp}, got {}", log.n_insertions);
        assert!(del_err < 0.05, "deletions: exp ~{exp}, got {}", log.n_deletions);
        // Length: n + insertions - deletions
        let expected_len = n + log.n_insertions - log.n_deletions;
        assert_eq!(s.sequence.len(), expected_len);
    }

    #[test]
    fn mutations_use_only_other_three_bases() {
        let mut s = flat_seq(b'A', 10_000);
        let g = Global {
            mutation_rate: 1.0,
            indel_rate: 0.0,
            ..Default::default()
        };
        let mut r = Streams::new(1).noise();
        let log = apply(&mut s, &g, &mut r);
        assert_eq!(log.n_substitutions, 10_000);
        for b in &s.sequence {
            assert_ne!(*b, b'A', "every base must have been substituted away from 'A'");
            assert!(matches!(*b, b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn noise_updates_filler_spans() {
        // F3 regression: after noise inserts/deletes bytes BEFORE a
        // SHIFT filler, the filler's realised_start_bp must move with
        // the sequence. Otherwise truth's phase_shift_positions point
        // at stale pre-noise coordinates.
        use crate::synth::blocks::{expand, FillerKind};
        use crate::synth::config::{Block, Source, Template};
        use crate::synth::templates::instantiate;
        use std::collections::HashMap;

        // Build a HOR with a SHIFT halfway through. Then run noise with
        // a guaranteed insertion at the very start (mutation_rate=0,
        // indel_rate forced via a high rate).
        let mut cfgs: HashMap<String, Template> = HashMap::new();
        cfgs.insert(
            "t".into(),
            Template::HOR_slots {
                monomer_length_bp: 100,
                k: 4,
                source: Source::Random,
                sequence: None,
                file: None,
                gc_content: 0.5,
                inter_slot_divergence: 0.10,
            },
        );
        let s = crate::synth::rng::Streams::new(42);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![
            Block::HOR { template: "t".into(), n_copies: 50 },
            Block::SHIFT { offset_bp: 25 },
            Block::HOR { template: "t".into(), n_copies: 50 },
        ];
        let mut rs = s.structure();
        let mut state = expand(&blocks, &inst, &mut rs).unwrap();
        let shift_pre = state.filler_spans.iter().find(|f| matches!(f.kind, FillerKind::Shift {..})).unwrap().realised_start_bp;

        // Now apply a noise pass with non-trivial indel_rate.
        let g = Global { mutation_rate: 0.0, indel_rate: 0.05, ..Default::default() };
        let mut rn = s.noise();
        let log = apply(&mut state, &g, &mut rn);

        let net_indels_before_shift: i64 = log
            .indels
            .iter()
            .filter(|(p, _)| *p < shift_pre)
            .map(|(_, d)| *d as i64)
            .sum();
        // Must have at least one indel for the test to be meaningful.
        // With seed=42 and 20_000 bp pre-shift at 5% indel rate, expect ~1000.
        assert!(net_indels_before_shift.abs() > 0, "test setup: noise produced no indels before shift");

        let shift_post = state.filler_spans.iter().find(|f| matches!(f.kind, FillerKind::Shift {..})).unwrap().realised_start_bp;
        let expected = (shift_pre as i64 + net_indels_before_shift) as usize;
        assert_eq!(
            shift_post, expected,
            "SHIFT filler position must shift with net indels before it (pre={shift_pre}, post={shift_post}, expected={expected})"
        );
    }

    #[test]
    fn coord_map_lengths_consistent_after_noise() {
        // Build a small HOR via the M2 path, run noise, then verify
        // that sum of coord_map lens + leading-flank shift == sequence
        // length.
        use crate::synth::blocks::expand;
        use crate::synth::config::Block;
        use crate::synth::templates::instantiate;
        use std::collections::HashMap;

        let mut cfgs: HashMap<String, Template> = HashMap::new();
        cfgs.insert(
            "t".into(),
            Template::HOR_slots {
                monomer_length_bp: 200,
                k: 4,
                source: Source::Random,
                sequence: None,
                file: None,
                gc_content: 0.5,
                inter_slot_divergence: 0.10,
            },
        );
        let s = Streams::new(42);
        let mut rt = s.templates();
        let inst = instantiate(&cfgs, &mut rt);
        let blocks = vec![Block::HOR {
            template: "t".into(),
            n_copies: 50,
        }];
        let mut rs = s.structure();
        let mut state = expand(&blocks, &inst, &mut rs).unwrap();
        let pre_len = state.sequence.len();
        assert_eq!(pre_len, 200 * 4 * 50);

        let g = Global {
            mutation_rate: 0.03,
            indel_rate: 0.02,
            ..Default::default()
        };
        let mut rn = s.noise();
        let log = apply(&mut state, &g, &mut rn);

        // Sum of coord_map lens should equal the post-noise sequence
        // length, since every base of every entry is covered (blocks
        // are contiguous) and no out-of-block bytes exist yet.
        let cm_total: usize = state.coord_map.entries.iter().map(|e| e.realised_len_bp).sum();
        assert_eq!(
            cm_total,
            state.sequence.len(),
            "coord_map lens ({cm_total}) must match post-noise sequence length ({}); \
             pre={pre_len}, +ins={}, -del={}",
            state.sequence.len(),
            log.n_insertions,
            log.n_deletions
        );

        // And expected length matches.
        assert_eq!(state.sequence.len(), pre_len + log.n_insertions - log.n_deletions);
    }
}

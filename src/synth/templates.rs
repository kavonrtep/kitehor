//! Template instantiation.
//!
//! For `HOR_slots`:
//! 1. Draw slot 1 (random with requested GC, or take the user-supplied
//!    sequence).
//! 2. Derive slots 2..k as independent per-base mutations of slot 1 at
//!    rate `inter_slot_divergence / 2`. The realised pairwise
//!    divergence approximates the requested target (small under-shoot
//!    at short monomer lengths — recorded in `realised_inter_slot_divergence`).
//!
//! Templates are instantiated **once** in deterministic
//! sorted-by-name order, so two structure blocks referencing the same
//! template name share slot consensuses byte-for-byte. That equality
//! is load-bearing: it's the structural distinction between *one
//! phase-shifted HOR* and *two unrelated HORs with the same k*.

use crate::synth::config::{Source, Template};
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct InstantiatedTemplate {
    /// One Vec per slot. `slots[0]` is slot 1 in the YAML's 1-indexing.
    pub slots: Vec<Vec<u8>>,
    /// Mean pairwise divergence across all `C(k, 2)` unordered slot
    /// pairs — i.e. the same metric the upstream `inter_slot_divergence`
    /// targets. Approximately equal to `inter_slot_divergence` in the
    /// large-k limit; somewhat lower for small k because the derived
    /// slots all share slot 1 as a common ancestor (see F5 in the
    /// implementation review).
    pub realised_inter_slot_divergence: f64,
}

pub fn instantiate(
    cfg_templates: &HashMap<String, Template>,
    rng: &mut ChaCha20Rng,
) -> HashMap<String, InstantiatedTemplate> {
    // Sorted-by-name order so the RNG sees a deterministic sequence of
    // draws regardless of the YAML's map iteration order.
    let mut names: Vec<&String> = cfg_templates.keys().collect();
    names.sort();
    let mut out = HashMap::with_capacity(names.len());
    for name in names {
        let tpl = &cfg_templates[name];
        let inst = match tpl {
            Template::HOR_slots {
                monomer_length_bp,
                k,
                source,
                sequence,
                gc_content,
                inter_slot_divergence,
                ..
            } => instantiate_hor_slots(
                *monomer_length_bp,
                *k,
                *source,
                sequence.as_deref(),
                *gc_content,
                *inter_slot_divergence,
                rng,
            ),
            Template::monomer {
                monomer_length_bp,
                source,
                sequence,
                gc_content,
                ..
            } => instantiate_monomer(
                *monomer_length_bp,
                *source,
                sequence.as_deref(),
                *gc_content,
                rng,
            ),
        };
        out.insert(name.clone(), inst);
    }
    out
}

fn random_dna(len: usize, gc: f64, rng: &mut ChaCha20Rng) -> Vec<u8> {
    let mut v = Vec::with_capacity(len);
    for _ in 0..len {
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

fn mutate(seed: &[u8], rate: f64, rng: &mut ChaCha20Rng) -> Vec<u8> {
    if rate <= 0.0 {
        return seed.to_vec();
    }
    let mut out = seed.to_vec();
    for b in out.iter_mut() {
        if rng.random::<f64>() < rate {
            let alts: &[u8] = match *b {
                b'A' => b"CGT",
                b'C' => b"AGT",
                b'G' => b"ACT",
                b'T' => b"ACG",
                _ => b"ACGT",
            };
            *b = alts[rng.random_range(0..alts.len())];
        }
    }
    out
}

fn upper_dna(s: &str) -> Vec<u8> {
    s.as_bytes()
        .iter()
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

fn instantiate_hor_slots(
    monomer_length_bp: usize,
    k: usize,
    source: Source,
    sequence: Option<&str>,
    gc_content: f64,
    inter_slot_divergence: f64,
    rng: &mut ChaCha20Rng,
) -> InstantiatedTemplate {
    let slot1: Vec<u8> = match source {
        Source::Random => random_dna(monomer_length_bp, gc_content, rng),
        Source::Sequence => upper_dna(sequence.expect("validated by load_and_validate")),
        Source::File => unreachable!("source: file rejected in MVP validator"),
    };
    let rate = inter_slot_divergence / 2.0;
    let mut slots = Vec::with_capacity(k);
    slots.push(slot1.clone());
    for _ in 1..k {
        slots.push(mutate(&slot1, rate, rng));
    }
    let realised = if k >= 2 {
        // F5: mean over all unordered slot pairs, matching the upstream
        // definition of inter_slot_divergence.
        let mut sum = 0.0;
        let mut pair_count = 0usize;
        let slot_len = slots[0].len() as f64;
        for i in 0..slots.len() {
            for j in (i + 1)..slots.len() {
                let diffs = slots[i]
                    .iter()
                    .zip(slots[j].iter())
                    .filter(|(a, b)| a != b)
                    .count();
                sum += diffs as f64 / slot_len;
                pair_count += 1;
            }
        }
        sum / pair_count as f64
    } else {
        0.0
    };
    InstantiatedTemplate {
        slots,
        realised_inter_slot_divergence: realised,
    }
}

fn instantiate_monomer(
    monomer_length_bp: usize,
    source: Source,
    sequence: Option<&str>,
    gc_content: f64,
    rng: &mut ChaCha20Rng,
) -> InstantiatedTemplate {
    let s = match source {
        Source::Random => random_dna(monomer_length_bp, gc_content, rng),
        Source::Sequence => upper_dna(sequence.expect("validated")),
        Source::File => unreachable!("source: file rejected in MVP"),
    };
    InstantiatedTemplate {
        slots: vec![s],
        realised_inter_slot_divergence: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::rng::Streams;

    fn one_hor_slots(
        monomer_length_bp: usize,
        k: usize,
        gc: f64,
        d: f64,
        seed: u64,
    ) -> InstantiatedTemplate {
        let mut rng = Streams::new(seed).templates();
        instantiate_hor_slots(monomer_length_bp, k, Source::Random, None, gc, d, &mut rng)
    }

    #[test]
    fn slot1_length_correct() {
        let t = one_hor_slots(171, 12, 0.5, 0.15, 42);
        assert_eq!(t.slots.len(), 12);
        for s in &t.slots {
            assert_eq!(s.len(), 171);
        }
    }

    #[test]
    fn deterministic_same_seed() {
        let a = one_hor_slots(200, 6, 0.45, 0.20, 7);
        let b = one_hor_slots(200, 6, 0.45, 0.20, 7);
        assert_eq!(a.slots, b.slots);
    }

    #[test]
    fn different_seed_differs() {
        let a = one_hor_slots(200, 4, 0.5, 0.15, 1);
        let b = one_hor_slots(200, 4, 0.5, 0.15, 2);
        assert_ne!(a.slots, b.slots);
    }

    #[test]
    fn gc_content_within_tolerance() {
        let t = one_hor_slots(2000, 1, 0.5, 0.0, 42);
        // 1-slot HOR_slots template is technically a degenerate case
        // but exercises the random_dna() path with a large sample.
        let n = t.slots[0].len() as f64;
        let gc = t.slots[0]
            .iter()
            .filter(|b| matches!(**b, b'G' | b'C'))
            .count() as f64
            / n;
        assert!(
            (gc - 0.5).abs() < 0.03,
            "expected GC ~ 0.5 +- 0.03; got {gc}"
        );
    }

    #[test]
    fn realised_divergence_in_band() {
        // F5: realised metric is mean pairwise across C(k, 2) pairs.
        //
        // Mutating slot 1 at rate d/2 to derive slots 2..k gives:
        //   slot 1 ↔ slot_i           : rate d/2  ((k-1) pairs)
        //   slot_i ↔ slot_j (i, j ≠ 1) : rate d   (C(k-1, 2) pairs)
        //
        // Mean over all C(k, 2) = k(k-1)/2 pairs:
        //   E[d_realised] = ((k-1)(d/2) + (k-1)(k-2)/2 · d) / (k(k-1)/2)
        //                 = d · (1/k + (k-2)/k)
        //                 = d · (k-1)/k
        //
        // For k=8, d=0.20 → E ≈ 0.175.
        let t = one_hor_slots(1000, 8, 0.5, 0.20, 99);
        let d = t.realised_inter_slot_divergence;
        let expected = 0.20 * 7.0 / 8.0;
        assert!(
            (d - expected).abs() < 0.03,
            "expected realised mean pairwise ≈ {expected} ± 0.03; got {d}"
        );
    }

    #[test]
    fn realised_divergence_approaches_d_for_large_k() {
        // Sanity check: at k=16, realised ≈ d · 15/16 ≈ 0.94 d.
        let t = one_hor_slots(2000, 16, 0.5, 0.10, 5);
        let d = t.realised_inter_slot_divergence;
        let expected = 0.10 * 15.0 / 16.0;
        assert!(
            (d - expected).abs() < 0.02,
            "expected realised ≈ {expected} ± 0.02; got {d}"
        );
    }

    #[test]
    fn zero_divergence_collapses_to_identical_slots() {
        let t = one_hor_slots(100, 4, 0.5, 0.0, 11);
        assert_eq!(t.realised_inter_slot_divergence, 0.0);
        for s in &t.slots[1..] {
            assert_eq!(s, &t.slots[0]);
        }
    }

    #[test]
    fn sequence_source_used_verbatim() {
        let mut rng = Streams::new(0).templates();
        let user_seq = "acgtacgtac";
        let t = instantiate_hor_slots(10, 3, Source::Sequence, Some(user_seq), 0.5, 0.0, &mut rng);
        assert_eq!(t.slots[0], b"ACGTACGTAC".to_vec());
    }

    #[test]
    fn instantiate_caches_by_name_deterministically() {
        use std::collections::HashMap;
        let mut cfgs: HashMap<String, Template> = HashMap::new();
        cfgs.insert(
            "alpha".to_string(),
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
        cfgs.insert(
            "beta".to_string(),
            Template::monomer {
                monomer_length_bp: 200,
                source: Source::Random,
                sequence: None,
                file: None,
                gc_content: 0.5,
            },
        );
        let mut r1 = Streams::new(7).templates();
        let mut r2 = Streams::new(7).templates();
        let a = instantiate(&cfgs, &mut r1);
        let b = instantiate(&cfgs, &mut r2);
        assert_eq!(a["alpha"].slots, b["alpha"].slots);
        assert_eq!(a["beta"].slots, b["beta"].slots);
    }

    #[test]
    fn monomer_template_has_k1() {
        let mut rng = Streams::new(1).templates();
        let t = instantiate_monomer(50, Source::Random, None, 0.5, &mut rng);
        assert_eq!(t.slots.len(), 1);
        assert_eq!(t.slots[0].len(), 50);
        assert_eq!(t.realised_inter_slot_divergence, 0.0);
    }
}

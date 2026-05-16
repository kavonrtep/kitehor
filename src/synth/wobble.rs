//! Wobble realisation via residual-accumulator integer edits.
//!
//! Two models:
//!
//! - **sinusoidal**: `δ(r) = amplitude_bp · sin(2π · r / period_rows)`.
//! - **random_walk**: cumulative-sum Gaussian noise smoothed with a
//!   moving average of width `period_rows / 4` (default 50 if
//!   `period_rows = 0`). σ is chosen so the smoothed series has std ≈
//!   `amplitude_bp`.
//!
//! Walks the sequence row by row, where one "row" = the monomer
//! length of the array's first repeat block. At each row boundary,
//! accumulates `δ(r) - δ(r-1)` into a residual; while the residual
//! exceeds ±1 bp, inserts or deletes a single base near the boundary.
//! This produces non-integer mean drift via integer base-level edits.
//!
//! Inserted bases come from local composition (last 50 bp). Deletions
//! remove the last byte of the row. `coord_map.apply_indels` updates
//! every entry to track the new positions.

use crate::synth::blocks::{FillerSpan, SimState};
use crate::synth::config::{Block, Config, Modifier, Template, WobbleModel};
use anyhow::{bail, Result};
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;

#[derive(Debug, Clone, Default)]
pub struct WobbleLog {
    /// Standard deviation of the realised δ(r) curve (the curve we
    /// actually emit, after integer-edit truncation).
    pub realised_amplitude_bp: f64,
    /// Recovered period in bp (= period_rows × monomer_length), or
    /// `None` for random-walk wobble.
    pub realised_periodicity_bp: Option<f64>,
    pub n_insertions: usize,
    pub n_deletions: usize,
}

pub fn apply(
    state: &mut SimState,
    modifiers: &[Modifier],
    cfg: &Config,
    rng: &mut ChaCha20Rng,
) -> Result<WobbleLog> {
    let mut log = WobbleLog::default();
    if modifiers.is_empty() {
        return Ok(log);
    }
    let row_size = first_repeat_monomer_len(cfg)?;
    for m in modifiers {
        if m.target != "all" {
            bail!(
                "wobble: target {:?} not supported in MVP (only 'all')",
                m.target
            );
        }
        apply_one(state, m, row_size, &mut log, rng)?;
    }
    Ok(log)
}

fn apply_one(
    state: &mut SimState,
    m: &Modifier,
    row_size: usize,
    log: &mut WobbleLog,
    rng: &mut ChaCha20Rng,
) -> Result<()> {
    let n_rows = state.sequence.len() / row_size;
    if n_rows == 0 || m.wobble.amplitude_bp <= 0.0 {
        return Ok(());
    }
    let amplitude = m.wobble.amplitude_bp;
    let period_rows = m.wobble.period_rows as f64;

    let delta = match m.wobble.model {
        WobbleModel::Sinusoidal => {
            let p = period_rows.max(1.0);
            (0..n_rows)
                .map(|r| amplitude * (2.0 * PI * r as f64 / p).sin())
                .collect::<Vec<f64>>()
        }
        WobbleModel::RandomWalk => {
            // AR(1) / Ornstein-Uhlenbeck-style process. y[r] = ρ·y[r-1] + ε,
            // ε ~ N(0, σ_step). Stationary std = σ_step / sqrt(1 - ρ²) ≈
            // σ_step · sqrt(W/2) for ρ = 1 − 1/W. Solving for σ_step:
            //     σ_step = amplitude · sqrt(2/W)
            // This bounds the wobble std at `amplitude_bp` regardless of
            // n_rows (a plain cumulative-sum walk would grow unbounded).
            let window = if period_rows > 0.0 {
                (period_rows / 4.0).round().max(2.0) as usize
            } else {
                50
            };
            let w = window as f64;
            let rho = 1.0 - 1.0 / w;
            let sigma_step = amplitude * (2.0 / w).sqrt();
            let normal = Normal::new(0.0, sigma_step).expect("valid normal");
            let mut y = vec![0.0; n_rows];
            let mut prev = 0.0;
            // A short burn-in lets the process reach stationary distribution
            // before we sample its values into the curve.
            for _ in 0..(3 * window) {
                prev = rho * prev + normal.sample(rng);
            }
            for r in 0..n_rows {
                prev = rho * prev + normal.sample(rng);
                y[r] = prev;
            }
            y
        }
    };

    // Convert δ into integer indels at row boundaries.
    let mut new_seq = Vec::with_capacity(state.sequence.len());
    let mut indels: Vec<(usize, i32)> = Vec::new();
    let mut residual = 0.0;
    for r in 0..n_rows {
        let row_start = r * row_size;
        let row_end = (r + 1) * row_size;
        new_seq.extend_from_slice(&state.sequence[row_start..row_end]);
        let target_change = if r == 0 { delta[0] } else { delta[r] - delta[r - 1] };
        residual += target_change;
        while residual >= 1.0 {
            // Insert a base at row_end (pre-wobble coord).
            new_seq.push(random_local_base(&new_seq, rng));
            indels.push((row_end, 1));
            residual -= 1.0;
            log.n_insertions += 1;
        }
        while residual <= -1.0 {
            // Delete the last base of the row.
            if new_seq.is_empty() {
                break;
            }
            new_seq.pop();
            indels.push((row_end - 1, -1));
            residual += 1.0;
            log.n_deletions += 1;
        }
    }
    // Carry over any tail bytes that didn't fit a whole row.
    if state.sequence.len() > n_rows * row_size {
        new_seq.extend_from_slice(&state.sequence[n_rows * row_size..]);
    }

    state.sequence = new_seq;
    state.coord_map.apply_indels(&indels);
    apply_to_fillers(&mut state.filler_spans, &indels);

    log.realised_amplitude_bp = std_dev(&delta);
    log.realised_periodicity_bp = match m.wobble.model {
        WobbleModel::Sinusoidal if period_rows > 0.0 => Some(period_rows * row_size as f64),
        _ => None,
    };
    Ok(())
}

fn apply_to_fillers(fillers: &mut [FillerSpan], indels: &[(usize, i32)]) {
    for fs in fillers.iter_mut() {
        let s = fs.realised_start_bp;
        let e = s + fs.realised_len_bp;
        let mut shift: i64 = 0;
        let mut len_delta: i64 = 0;
        for &(pos, delta) in indels {
            if pos < s {
                shift += delta as i64;
            } else if pos < e {
                len_delta += delta as i64;
            }
        }
        fs.realised_start_bp = (s as i64 + shift) as usize;
        fs.realised_len_bp = ((fs.realised_len_bp as i64) + len_delta).max(0) as usize;
    }
}

fn first_repeat_monomer_len(cfg: &Config) -> Result<usize> {
    for b in &cfg.structure {
        if let Block::HOR { template, .. } | Block::SIMPLE_TR { template, .. } = b {
            return match cfg.templates.get(template) {
                Some(Template::HOR_slots {
                    monomer_length_bp, ..
                })
                | Some(Template::monomer {
                    monomer_length_bp, ..
                }) => Ok(*monomer_length_bp),
                None => bail!("wobble: template `{template}` not found"),
            };
        }
    }
    bail!("wobble: no HOR/SIMPLE_TR block in structure (cannot determine row size)")
}

fn random_local_base(out: &[u8], rng: &mut ChaCha20Rng) -> u8 {
    let window: usize = 50;
    let src: &[u8] = if out.len() >= window {
        &out[out.len() - window..]
    } else if !out.is_empty() {
        out
    } else {
        return b"ACGT"[rng.random_range(0..4)];
    };
    src[rng.random_range(0..src.len())]
}

fn std_dev(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mean = xs.iter().sum::<f64>() / xs.len() as f64;
    let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / xs.len() as f64;
    var.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::{blocks::expand, config::Source, rng::Streams, templates::instantiate};
    use std::collections::HashMap;
    use std::io::Write;

    fn parse(yaml: &str) -> Config {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        crate::synth::config::load_and_validate(f.path()).unwrap()
    }

    fn build_state(cfg: &Config, seed: u64) -> SimState {
        let s = Streams::new(seed);
        let mut rt = s.templates();
        let inst = instantiate(&cfg.templates, &mut rt);
        let mut rs = s.structure();
        expand(&cfg.structure, &inst, &mut rs).unwrap()
    }

    fn long_hor(n_copies: usize) -> Config {
        let yaml = format!(
            r#"
schema_version: 1
seed: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.10
structure:
  - type: HOR
    template: t
    n_copies: {n_copies}
"#
        );
        parse(&yaml)
    }

    #[test]
    fn zero_amplitude_is_noop() {
        let cfg = long_hor(50);
        let mut state = build_state(&cfg, 1);
        let pre_len = state.sequence.len();
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 0.0,
                period_rows: 100,
                model: WobbleModel::Sinusoidal,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        let log = apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        assert_eq!(state.sequence.len(), pre_len);
        assert_eq!(log.n_insertions, 0);
        assert_eq!(log.n_deletions, 0);
    }

    #[test]
    fn sinusoidal_realised_amplitude_close_to_target() {
        let cfg = long_hor(5000); // long enough for stable stats
        let mut state = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 1.5,
                period_rows: 500,
                model: WobbleModel::Sinusoidal,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        let log = apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        // Realised amplitude is std of sin curve with peak 1.5 → std = 1.5/sqrt(2) ≈ 1.06.
        let expected = 1.5 / std::f64::consts::SQRT_2;
        let err = (log.realised_amplitude_bp - expected).abs() / expected;
        assert!(
            err < 0.15,
            "expected sinusoidal std ≈ {expected:.3} ± 15%, got {:.3}",
            log.realised_amplitude_bp
        );
    }

    #[test]
    fn random_walk_realised_amplitude_close_to_target() {
        let cfg = long_hor(5000);
        let mut state = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 2.0,
                period_rows: 0,
                model: WobbleModel::RandomWalk,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        let log = apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        assert!(
            log.realised_amplitude_bp > 0.5 && log.realised_amplitude_bp < 4.0,
            "random_walk realised std {:.3} should be near 2.0",
            log.realised_amplitude_bp
        );
    }

    #[test]
    fn determinism_same_seed() {
        let cfg = long_hor(2000);
        let mut a = build_state(&cfg, 1);
        let mut b = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 1.0,
                period_rows: 300,
                model: WobbleModel::Sinusoidal,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut r1 = Streams::new(7).wobble();
        let mut r2 = Streams::new(7).wobble();
        apply(&mut a, &[m.clone_for_test()], &cfg, &mut r1).unwrap();
        apply(&mut b, &[m], &cfg, &mut r2).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }

    // Helper for test cloning since Modifier doesn't derive Clone.
    impl Modifier {
        fn clone_for_test(&self) -> Modifier {
            Modifier {
                target: self.target.clone(),
                wobble: crate::synth::config::WobbleSpec {
                    amplitude_bp: self.wobble.amplitude_bp,
                    period_rows: self.wobble.period_rows,
                    model: self.wobble.model,
                    realisation: self.wobble.realisation,
                },
            }
        }
    }

    #[test]
    fn coord_map_remains_consistent_after_wobble() {
        // Single HOR block, wobble applied. Sum of coord_map lens
        // should still equal sequence length (no fillers exist).
        let cfg = long_hor(500);
        let mut state = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 1.5,
                period_rows: 200,
                model: WobbleModel::Sinusoidal,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        let cm_total: usize = state.coord_map.entries.iter().map(|e| e.realised_len_bp).sum();
        assert_eq!(cm_total, state.sequence.len());
    }

    #[test]
    fn sinusoidal_periodicity_in_bp() {
        // period_rows × monomer_length should be reported as
        // realised_periodicity_bp.
        let cfg = long_hor(2000);
        let mut state = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 1.0,
                period_rows: 500,
                model: WobbleModel::Sinusoidal,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        let log = apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        // monomer_length = 100 → 500 rows × 100 bp = 50_000 bp.
        assert_eq!(log.realised_periodicity_bp, Some(50_000.0));
    }

    #[test]
    fn random_walk_periodicity_is_none() {
        let cfg = long_hor(1000);
        let mut state = build_state(&cfg, 1);
        let m = Modifier {
            target: "all".into(),
            wobble: crate::synth::config::WobbleSpec {
                amplitude_bp: 1.0,
                period_rows: 0,
                model: WobbleModel::RandomWalk,
                realisation: crate::synth::config::WobbleRealisation::IntegerEdits,
            },
        };
        let mut rng = Streams::new(1).wobble();
        let log = apply(&mut state, &[m], &cfg, &mut rng).unwrap();
        assert!(log.realised_periodicity_bp.is_none());
    }

    #[test]
    fn requires_repeat_block_for_row_size() {
        // Empty structure isn't allowed by schema; the validator
        // forbids `structure: []` (minItems=1). Try the next-best
        // proxy: a structure with only INSERTIONs.
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
        let mut cfg = parse(
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
"#,
        );
        // Replace the structure with one that has no repeat block.
        cfg.structure = vec![Block::INSERTION {
            length_bp: 100,
            kind: crate::synth::config::InsertionKind::Random,
        }];
        let res = first_repeat_monomer_len(&cfg);
        assert!(res.is_err());
    }
}

//! Period candidate generator (`periods.tsv`).
//!
//! Mimics what an upstream period generator (TideHunter / TideCluster /
//! NTRprism) would feed the detector. For each unique repeat block in
//! the config, we emit:
//!
//! - the true base width (high score)
//! - the true HOR unit length if `k > 1` (high score)
//! - a near-miss `±[2..4]` bp (mid score)
//! - a `2× base_width` or `3× base_width` harmonic (mid score)
//! - one false positive in 100..5000 bp not matching any real period
//!   (low score)
//!
//! Score values are documentation only. Distinct values prevent the
//! corpus from being unrealistically clean.

use crate::synth::config::{Block, Config, Template};
use crate::synth::fasta::with_ext;
use anyhow::{Context, Result};
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct PeriodCandidate {
    pub period_bp: usize,
    pub period_score: f64,
    pub source: &'static str,
}

pub fn build(cfg: &Config, array_id: &str, rng: &mut ChaCha20Rng) -> Vec<PeriodCandidate> {
    let mut out = Vec::new();
    let mut real: BTreeSet<usize> = BTreeSet::new();

    for b in &cfg.structure {
        match b {
            Block::HOR {
                template,
                n_copies: _,
            } => {
                if let Some(Template::HOR_slots {
                    monomer_length_bp,
                    k,
                    inter_slot_divergence,
                    ..
                }) = cfg.templates.get(template)
                {
                    if real.insert(*monomer_length_bp) {
                        out.push(PeriodCandidate {
                            period_bp: *monomer_length_bp,
                            period_score: 0.94,
                            source: "true_base",
                        });
                    }
                    // HOR unit only when div > 0 — otherwise HOR is
                    // effectively a SIMPLE_TR.
                    if *inter_slot_divergence > 0.0 && *k >= 2 {
                        let unit = monomer_length_bp * k;
                        if real.insert(unit) {
                            out.push(PeriodCandidate {
                                period_bp: unit,
                                period_score: 0.88,
                                source: "true_hor_unit",
                            });
                        }
                    }
                }
            }
            Block::SIMPLE_TR { template, .. } => {
                let l = match cfg.templates.get(template) {
                    Some(Template::HOR_slots {
                        monomer_length_bp, ..
                    })
                    | Some(Template::monomer {
                        monomer_length_bp, ..
                    }) => *monomer_length_bp,
                    None => continue,
                };
                if real.insert(l) {
                    out.push(PeriodCandidate {
                        period_bp: l,
                        period_score: 0.94,
                        source: "true_base",
                    });
                }
            }
            _ => {}
        }
    }

    // Add a near-miss and harmonic for the **first** discovered base
    // width — that's the primary period for HOR detection purposes.
    if let Some(&base) = real.iter().next() {
        let nm_offset =
            rng.random_range(2..=4) as i64 * (if rng.random::<f64>() < 0.5 { 1 } else { -1 });
        let nm = (base as i64 + nm_offset).max(20) as usize;
        if !real.contains(&nm) {
            out.push(PeriodCandidate {
                period_bp: nm,
                period_score: 0.71,
                source: "near_miss",
            });
        }
        let harm = if rng.random::<f64>() < 0.5 {
            2 * base
        } else {
            3 * base
        };
        if !real.contains(&harm) {
            out.push(PeriodCandidate {
                period_bp: harm,
                period_score: 0.65,
                source: "harmonic",
            });
        }
        // One low-score false positive.
        let fp = loop {
            let p = rng.random_range(100..=5000);
            if !real.contains(&p) {
                break p;
            }
        };
        out.push(PeriodCandidate {
            period_bp: fp,
            period_score: 0.42,
            source: "false_positive",
        });
    }

    let _ = array_id;
    out
}

pub fn write(prefix: &Path, array_id: &str, candidates: &[PeriodCandidate]) -> Result<()> {
    let path = with_ext(prefix, "periods.tsv");
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {:?}", path))?;
    writeln!(f, "array_id\tperiod_bp\tperiod_score\tsource")?;
    for c in candidates {
        writeln!(
            f,
            "{}\t{}\t{:.4}\t{}",
            array_id, c.period_bp, c.period_score, c.source
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::rng::Streams;
    use std::io::Write;

    fn parse(yaml: &str) -> Config {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        crate::synth::config::load_and_validate(f.path()).unwrap()
    }

    #[test]
    fn hor_emits_base_and_unit() {
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
        let mut rng = Streams::new(0).structure();
        let ps = build(&cfg, "arr", &mut rng);
        let sources: Vec<&str> = ps.iter().map(|p| p.source).collect();
        assert!(sources.contains(&"true_base"));
        assert!(sources.contains(&"true_hor_unit"));
        assert!(sources.contains(&"near_miss"));
        assert!(sources.contains(&"harmonic"));
        assert!(sources.contains(&"false_positive"));
        let base = ps.iter().find(|p| p.source == "true_base").unwrap();
        assert_eq!(base.period_bp, 171);
        let unit = ps.iter().find(|p| p.source == "true_hor_unit").unwrap();
        assert_eq!(unit.period_bp, 171 * 12);
    }

    #[test]
    fn simple_tr_emits_no_hor_unit() {
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
        let mut rng = Streams::new(0).structure();
        let ps = build(&cfg, "arr", &mut rng);
        let sources: Vec<&str> = ps.iter().map(|p| p.source).collect();
        assert!(sources.contains(&"true_base"));
        assert!(!sources.contains(&"true_hor_unit"));
    }

    #[test]
    fn coexisting_periods_both_listed() {
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
        let mut rng = Streams::new(0).structure();
        let ps = build(&cfg, "arr", &mut rng);
        let bases: Vec<usize> = ps
            .iter()
            .filter(|p| p.source == "true_base")
            .map(|p| p.period_bp)
            .collect();
        assert!(bases.contains(&171));
        assert!(bases.contains(&200));
    }
}

//! v2 YAML-driven synthetic tandem-repeat array generator.
//!
//! Design lives in `docs/new/simulator_impl_plan.md` (kitehor-specific
//! implementation contract) and its companions: `taxonomy.md`,
//! `simulator_plan.md`, `simulator_schema.json`, `detect_spec.md`.
//!
//! Built in milestones M1..M7 per `simulator_impl_plan.md §9`. This
//! ships M1: schema embed + config loader + serde-first validation
//! with MVP business rules (A1 event block targeting, Q5 negative
//! SHIFT bounds, Q8 source:file rejection).

pub mod blocks;
pub mod config;
pub mod coords;
pub mod diagnostics;
pub mod events;
pub mod fasta;
pub mod grammar;
pub mod noise;
pub mod periods;
pub mod rng;
pub mod templates;
pub mod truth;
pub mod wobble;

/// Canonical JSON Schema embedded at build time. The file under
/// `src/synth/simulator.schema.json` is an embedded copy of
/// `docs/new/simulator_schema.json`. Drift between the two is caught
/// by `tests/synth_schema_drift.rs`.
pub const CANONICAL_SCHEMA: &str = include_str!("simulator.schema.json");

pub use blocks::{expand, SimState};
pub use config::{load_and_validate, Config, ConfigError};
pub use coords::{CoordEntry, CoordMap};
pub use rng::Streams;
pub use templates::{instantiate, InstantiatedTemplate};

use anyhow::Result;
use rayon::prelude::*;
use std::path::{Path, PathBuf};

/// End-to-end pipeline for a single YAML config. Writes
/// `{prefix}.fa`, `{prefix}.truth.tsv`, `{prefix}.periods.tsv`, and
/// optionally `{prefix}.diagnostics.json`.
pub fn run_one(
    cfg_path: &Path,
    out_prefix: &Path,
    seed_override: Option<u64>,
    diagnostics: bool,
) -> Result<()> {
    let cfg = load_and_validate(cfg_path).map_err(|e| anyhow::anyhow!("{}", e))?;
    let top_seed = seed_override.unwrap_or(cfg.seed);
    let streams = Streams::new(top_seed);

    let mut rt = streams.templates();
    let inst = instantiate(&cfg.templates, &mut rt);

    let mut rs = streams.structure();
    let mut state = expand(&cfg.structure, &inst, &mut rs)?;

    let mut rw = streams.wobble();
    let wobble_log = wobble::apply(&mut state, &cfg.modifiers, &cfg, &mut rw)?;

    let mut re = streams.events();
    let event_logs = events::apply(&mut state, &cfg.post_generation, &cfg, &inst, &mut re)?;

    let mut rn = streams.noise();
    let noise_log = noise::apply(&mut state, &cfg.global, &mut rn);

    let array_id = cfg
        .global
        .array_id
        .clone()
        .unwrap_or_else(|| array_id_from_prefix(out_prefix));

    fasta::write(out_prefix, &array_id, &state.sequence)?;
    let truth_row =
        truth::build_truth(&cfg, &state, &array_id, &noise_log, &wobble_log, &event_logs);
    truth::write(out_prefix, &truth_row)?;

    let mut rp = streams.structure(); // period candidates piggy-back on the structure stream
    let pers = periods::build(&cfg, &array_id, &mut rp);
    periods::write(out_prefix, &array_id, &pers)?;

    if diagnostics {
        diagnostics::write(
            out_prefix,
            cfg_path,
            top_seed,
            &state,
            &inst,
            &wobble_log,
            &event_logs,
            &noise_log,
            &cfg,
        )?;
    }

    Ok(())
}

/// Iterate every `*.yaml` under `config_dir` (skipping `.deferred.yaml`),
/// run `run_one` for each in parallel via rayon, and write outputs to
/// `out_dir/<stem>.*`. Per-config seed is derived as `seed_offset XOR
/// fnv1a(filename)` so different files produce different sequences
/// while a fixed `seed_offset` keeps the corpus byte-reproducible.
///
/// Returns the number of configs successfully run.
pub fn run_batch(
    config_dir: &Path,
    out_dir: &Path,
    seed_offset: u64,
    diagnostics: bool,
) -> Result<usize> {
    std::fs::create_dir_all(out_dir)?;
    let configs = discover_configs(config_dir)?;
    configs.par_iter().try_for_each(|p| -> Result<()> {
        let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("array");
        let prefix = out_dir.join(stem);
        let per_seed = rng::derive(seed_offset, stem);
        run_one(p, &prefix, Some(per_seed), diagnostics)
    })?;
    Ok(configs.len())
}

fn discover_configs(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        let name = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let is_yaml = p.extension().and_then(|x| x.to_str()) == Some("yaml");
        let is_deferred = name.ends_with(".deferred.yaml");
        if is_yaml && !is_deferred {
            out.push(p);
        }
    }
    out.sort();
    Ok(out)
}

fn array_id_from_prefix(p: &Path) -> String {
    p.file_name()
        .and_then(|n| n.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| "array".to_string())
}

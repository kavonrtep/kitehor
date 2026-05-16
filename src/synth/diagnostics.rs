//! Per-array `diagnostics.json` writer.
//!
//! Emitted only when `synth --diagnostics` is set. Carries information
//! that the truth file deliberately omits (RNG seeds, realised
//! template slot bytes, per-block coordinates, noise/wobble counts).
//! Schema follows `docs/new/simulator_plan.md §5.4`.

use crate::synth::blocks::{FillerKind, SimState};
use crate::synth::config::Config;
use crate::synth::events::EventLog;
use crate::synth::fasta::with_ext;
use crate::synth::noise::NoiseLog;
use crate::synth::rng::derive;
use crate::synth::templates::InstantiatedTemplate;
use crate::synth::wobble::WobbleLog;
use anyhow::Result;
use serde_json::json;
use std::collections::HashMap;
use std::path::Path;

pub fn write(
    prefix: &Path,
    cfg_path: &Path,
    top_seed: u64,
    state: &SimState,
    templates: &HashMap<String, InstantiatedTemplate>,
    wobble: &WobbleLog,
    events: &[EventLog],
    noise: &NoiseLog,
    cfg: &Config,
) -> Result<()> {
    let tmpl_json: serde_json::Value = templates
        .iter()
        .map(|(name, inst)| {
            (
                name.clone(),
                json!({
                    "slots": inst.slots.iter()
                        .map(|s| std::str::from_utf8(s).unwrap_or("").to_string())
                        .collect::<Vec<_>>(),
                    "realised_inter_slot_divergence": inst.realised_inter_slot_divergence,
                }),
            )
        })
        .collect::<serde_json::Map<_, _>>()
        .into();

    let blocks_json: Vec<serde_json::Value> = cfg
        .structure
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let kind = match b {
                crate::synth::config::Block::HOR { .. } => "HOR",
                crate::synth::config::Block::SIMPLE_TR { .. } => "SIMPLE_TR",
                crate::synth::config::Block::SHIFT { .. } => "SHIFT",
                crate::synth::config::Block::INSERTION { .. } => "INSERTION",
            };
            // Block extent in realised coords: union of all entries +
            // filler_spans whose block_idx == i.
            let mut s = usize::MAX;
            let mut e = 0usize;
            let mut covered = false;
            for entry in &state.coord_map.entries {
                if entry.block_idx == i {
                    s = s.min(entry.realised_start_bp);
                    e = e.max(entry.realised_start_bp + entry.realised_len_bp);
                    covered = true;
                }
            }
            for fs in &state.filler_spans {
                if fs.block_idx == i {
                    s = s.min(fs.realised_start_bp);
                    e = e.max(fs.realised_start_bp + fs.realised_len_bp);
                    covered = true;
                }
            }
            if !covered {
                s = 0;
                e = 0;
            }
            json!({ "index": i, "type": kind, "start_bp": s, "end_bp": e })
        })
        .collect();

    let filler_json: Vec<serde_json::Value> = state
        .filler_spans
        .iter()
        .map(|fs| {
            let kind = match fs.kind {
                FillerKind::Shift { offset_bp } => {
                    json!({"type": "SHIFT", "offset_bp": offset_bp})
                }
                FillerKind::Insertion(k) => json!({
                    "type": "INSERTION",
                    "kind": format!("{:?}", k),
                }),
            };
            json!({
                "block_idx": fs.block_idx,
                "kind": kind,
                "start_bp": fs.realised_start_bp,
                "length_bp": fs.realised_len_bp,
            })
        })
        .collect();

    let events_json: Vec<serde_json::Value> = events
        .iter()
        .map(|e| {
            let s = crate::synth::events::to_events_json(std::slice::from_ref(e));
            serde_json::from_str::<serde_json::Value>(&s)
                .map(|v| v[0].clone())
                .unwrap_or(serde_json::Value::Null)
        })
        .collect();

    let out = json!({
        "config_path": cfg_path.display().to_string(),
        "rng_seeds": {
            "top":       top_seed,
            "templates": derive(top_seed, "templates"),
            "structure": derive(top_seed, "structure"),
            "wobble":    derive(top_seed, "wobble"),
            "events":    derive(top_seed, "events"),
            "noise":     derive(top_seed, "noise"),
        },
        "templates": tmpl_json,
        "blocks": blocks_json,
        "filler_spans": filler_json,
        "wobble_realised": {
            "amplitude_bp_std": wobble.realised_amplitude_bp,
            "periodicity_bp":   wobble.realised_periodicity_bp,
            "n_insertions":     wobble.n_insertions,
            "n_deletions":      wobble.n_deletions,
        },
        "events": events_json,
        "noise": {
            "n_substitutions": noise.n_substitutions,
            "n_insertions":    noise.n_insertions,
            "n_deletions":     noise.n_deletions,
        },
        "sequence_length_bp": state.sequence.len(),
    });

    let path = with_ext(prefix, "diagnostics.json");
    let s = serde_json::to_string_pretty(&out)?;
    std::fs::write(path, s)?;
    Ok(())
}

//! `kitehor` — kite-first sequence-agnostic HOR detector.
//!
//! The library exposes the building blocks of the CLI:
//!
//! - [`kite`] — k-mer immediate-neighbour histogram and peak scoring,
//!   a Rust port of `TideCluster/tarean/kite.R`.
//! - [`rule_classify`] — HOR / simple_tr / unresolved classifier on
//!   kite peaks (port of `tools/rule_proto/rule_proto.py`).
//! - [`subrepeat`] — spatial alternation / nested-TR detector (port
//!   of `tools/rule_proto/subrepeat_scan.py`).
//! - [`ssr`] — TideCluster-style SSR scan with kite-driven consensus
//!   (port of `tools/rule_proto/ssr_scan.py`).
//! - [`hor_validate`] — within-tile + spatial-density HOR validator
//!   (port of `tools/rule_proto/hor_within_tile_check.py`).
//! - [`summary`] — 8-rule combined-class merger (port of
//!   `tools/rule_proto/summary.py`).
//! - [`analyze`] — end-to-end orchestrator running all five stages.
//! - [`monomer_model`] — block-mean homology probe (`probe_period`).
//! - [`simulate`] / [`simulate_grid`] / [`synth`] — synthetic
//!   HOR / tandem arrays for testing.

pub mod cli;
pub mod errors;
pub mod io;
pub mod sequence;

pub mod kite;
pub mod monomer_model;
pub mod rule_classify;
pub mod simulate;
pub mod simulate_grid;
pub mod ssr;
pub mod subrepeat;
pub mod summary;
pub mod synth;

pub mod analyze;
pub mod detect;
pub mod emit_periods;
pub mod hor_validate;

pub use errors::{HordetectError, Result};

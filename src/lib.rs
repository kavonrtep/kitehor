//! `kitehor` — kite-first sequence-agnostic HOR detector.
//!
//! The library exposes the building blocks of the CLI:
//!
//! - [`kite`] — k-mer immediate-neighbour histogram and peak scoring,
//!   a Rust port of `TideCluster/tarean/kite.R`.
//! - [`rule_classify`] — HOR / simple_tr / unresolved classifier on
//!   kite peaks (port of `tools/rule_proto/rule_proto.py`).
//! - [`tandem_validate`] — unified spatial-localization subrepeat
//!   detector (port of `tools/rule_proto/tandem_validate.py`, spec v5);
//!   replaces the older `hor_within_tile_check` + `subrepeat_scan`
//!   prototypes.
//! - [`ssr`] — TideCluster-style SSR scan with kite-driven consensus
//!   (port of `tools/rule_proto/ssr_scan.py`).
//! - [`summary`] — 7-rule combined-class merger (port of
//!   `tools/rule_proto/summary_unified.py`).
//! - [`irregularity`] — indel-event scan (distance-residual + phase-bin
//!   clustering; port of `tools/rule_proto/irregularity_v2.py`).
//! - [`analyze`] — end-to-end orchestrator running all five stages.
//! - [`monomer_model`] — block-mean homology probe (`probe_period`).
//! - [`simulate`] / [`simulate_grid`] / [`synth`] — synthetic
//!   HOR / tandem arrays for testing.

pub mod cli;
pub mod errors;
pub mod io;
pub mod sequence;

pub mod irregularity;
pub mod kite;
pub mod monomer_model;
pub mod periodogram;
pub mod report;
pub mod rule_classify;
pub mod simulate;
pub mod simulate_grid;
pub mod ssr;
pub mod summary;
pub mod synth;

pub mod analyze;
pub mod detect;
pub mod emit_periods;
pub mod tandem_validate;

pub use errors::{HordetectError, Result};

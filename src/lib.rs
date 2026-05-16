//! `kitehor` — kite-first probabilistic Higher-Order Repeat detector.
//!
//! The library exposes the building blocks of the CLI:
//!
//! - [`kite`] — k-mer immediate-neighbour histogram and peak scoring,
//!   a Rust port of `TideCluster/tarean/kite.R`.
//! - [`features`] — per-record feature extraction for the classifier.
//! - [`classifier`] — random-forest loader, Platt scaler, config.
//! - [`classify`] — verdict orchestrator: features → RF → Platt →
//!   4-category + family demotion + k-recovery.
//! - [`hor_call`] — rule-based HOR-vs-tandem layer (legacy fallback).
//! - [`monomer_model`] — block-mean homology probe (`probe_period`).
//! - [`simulate`] / [`simulate_grid`] — synthetic HOR / tandem arrays
//!   for testing and training-set construction.

pub mod cli;
pub mod errors;
pub mod io;
pub mod sequence;

pub mod classifier;
pub mod classify;
pub mod coverage;
pub mod features;
pub mod hor_call;
pub mod kite;
pub mod monomer_model;
pub mod rule;
pub mod simulate;
pub mod simulate_grid;
pub mod synth;

pub mod detect;

pub use errors::{HordetectError, Result};

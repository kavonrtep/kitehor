//! M7.2 same-width mixed override (`docs/new/detect_m7_plan.md` §6).
//!
//! Locked-in tests:
//!
//! - **Positive:** T20 — two HOR blocks at the same `(base_width, k)`
//!   with different monomer templates. The pre-M7 detector collapsed
//!   these to a single HOR call because no upstream signal
//!   distinguishes same-`k` mixed families; M7.2's per-block
//!   consensus-identity override correctly reports `mixed`.
//! - **Negative `T15_stratification`:** two `SIMPLE_TR` blocks with the
//!   same period but different monomers stay `simple_TR` — the M7.2
//!   override is gated to `HOR / IrregularHOR` only.
//! - **Negative `T05_hor_clean`:** clean HOR stays `HOR`.
//! - **Negative `T03_wobble_aperiodic`:** wobble stays as documented
//!   (`simple_TR` per `detect_expectations.tsv` — wobble at k=1 is
//!   a property, not a class change).
//! - **Negative `T10_phase_shift`:** phase-shift HOR stays `HOR`.
//! - **Schema:** segments.tsv emits the 13-column M7.2 header.

use std::path::{Path, PathBuf};
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("synth_configs")
}

fn synth_and_detect(stem: &str) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let synth = dir.path().join(stem);
    let cfg = corpus_dir().join(format!("{stem}.yaml"));
    let o = Command::new(kitehor_bin())
        .args(["synth"])
        .arg(&cfg)
        .arg("-o")
        .arg(&synth)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "synth failed for {stem}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let mut periods = synth.as_os_str().to_owned();
    periods.push(".periods.tsv");
    let det = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(synth.with_extension("fa"))
        .arg("--periods")
        .arg(PathBuf::from(periods))
        .arg("-o")
        .arg(&det)
        .arg("--allow-missing-periods")
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "detect failed for {stem}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    (dir, det)
}

fn class_of(det: &Path) -> String {
    let mut p = det.as_os_str().to_owned();
    p.push(".properties.tsv");
    let s = std::fs::read_to_string(PathBuf::from(p)).unwrap();
    let row = s.lines().nth(1).unwrap();
    row.split('\t').nth(2).unwrap().to_string()
}

fn segments_header(det: &Path) -> String {
    let mut p = det.as_os_str().to_owned();
    p.push(".segments.tsv");
    let s = std::fs::read_to_string(PathBuf::from(p)).unwrap();
    s.lines().next().unwrap().to_string()
}

// ---- positive ----

#[test]
fn t20_same_width_mixed_hor_fires_override() {
    let (_dir, det) = synth_and_detect("T20_same_width_mixed_hor");
    let class = class_of(&det);
    assert_eq!(
        class, "mixed",
        "T20 (two HOR blocks at same base_width+k, different monomers) \
         should fire the M7.2 mixed override; got `{class}`"
    );
}

// ---- negatives: M7.2 must NOT fire on these ----

#[test]
fn t15_stratification_stays_simple_tr() {
    let (_dir, det) = synth_and_detect("T15_stratification");
    let class = class_of(&det);
    assert_eq!(
        class, "simple_TR",
        "T15 (two simple_TR blocks, different monomers, same period) \
         must stay simple_TR — M7.2 override is HOR/IrregularHOR only; got `{class}`"
    );
}

#[test]
fn t05_clean_hor_stays_hor() {
    let (_dir, det) = synth_and_detect("T05_hor_clean");
    let class = class_of(&det);
    assert_eq!(class, "HOR", "T05 (clean HOR) must stay HOR; got `{class}`");
}

#[test]
fn t10_phase_shift_stays_hor() {
    let (_dir, det) = synth_and_detect("T10_phase_shift");
    let class = class_of(&det);
    assert_eq!(
        class, "HOR",
        "T10 (phase-shift HOR, single family) must stay HOR; got `{class}`"
    );
}

#[test]
fn t03_wobble_aperiodic_stays_simple_tr() {
    // T03 has expected_class=simple_TR in detect_expectations.tsv —
    // wobble at k=1 is a property, not a class change. M7.2 must
    // not perturb this.
    let (_dir, det) = synth_and_detect("T03_wobble_aperiodic");
    let class = class_of(&det);
    assert_eq!(
        class, "simple_TR",
        "T03 (aperiodic wobble at k=1) must stay simple_TR; got `{class}`"
    );
}

// ---- schema ----

#[test]
fn segments_tsv_has_m7_2_columns() {
    let (_dir, det) = synth_and_detect("T05_hor_clean");
    let header = segments_header(&det);
    let cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(
        cols.len(),
        13,
        "segments.tsv header should have 13 columns post-M7.2; got {} ({header})",
        cols.len()
    );
    assert!(
        cols.contains(&"consensus_identity_to_reference"),
        "expected `consensus_identity_to_reference` column; header was `{header}`"
    );
    assert!(
        cols.contains(&"consensus_identity_coverage"),
        "expected `consensus_identity_coverage` column; header was `{header}`"
    );
}

// ---- M7.3 ----

#[test]
fn diagnostics_json_schema_version_is_2() {
    let (_dir, det) = synth_and_detect("T05_hor_clean");
    let mut p = det.as_os_str().to_owned();
    p.push(".diagnostics.json");
    let s = std::fs::read_to_string(PathBuf::from(p)).unwrap();
    let doc: serde_json::Value = serde_json::from_str(&s).unwrap();
    assert_eq!(
        doc.get("schema_version").and_then(|v| v.as_u64()),
        Some(2),
        "expected diagnostics schema_version 2 after M7.3"
    );
}

#[test]
fn mixed_array_emits_per_segment_monomers_to_consensus_fa() {
    // T20 is the mixed positive control — class becomes Mixed via
    // the M7.2 override. M7.3 says we should emit per-block monomer
    // consensus records (`<array_id>_seg{N}_monomer`) instead of the
    // whole-array .monomer (which doesn't exist for mixed).
    let (_dir, det) = synth_and_detect("T20_same_width_mixed_hor");
    let class = class_of(&det);
    assert_eq!(
        class, "mixed",
        "T20 should be mixed for this assertion to apply"
    );

    let mut p = det.as_os_str().to_owned();
    p.push(".consensus.fa");
    let s = std::fs::read_to_string(PathBuf::from(p)).unwrap();
    // Per-segment monomer records present.
    assert!(
        s.contains("_seg1_monomer"),
        "expected per-segment monomer record in consensus.fa; got:\n{s}"
    );
    // Whole-array monomer / hor_unit must NOT be present for mixed.
    let array_dot_monomer = format!(">{}.monomer", "T20_same_width_mixed_hor");
    assert!(
        !s.contains(&array_dot_monomer),
        "mixed array must not emit whole-array .monomer record; got:\n{s}"
    );
}

//! M3.5 acceptance: Pass-B phase-shift offset recovery.
//!
//! Runs `kitehor detect` on the T10 phase-shift fixture (offset 85 bp
//! on a 171 bp monomer) and asserts the recovered
//! `phase_shift_offsets` is within ±5 bp of truth.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn synth_then_detect(yaml: &str) -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let cfg = dir.path().join("cfg.yaml");
    std::fs::File::create(&cfg)
        .unwrap()
        .write_all(yaml.as_bytes())
        .unwrap();
    let synth = dir.path().join("arr");
    let o = Command::new(kitehor_bin())
        .args(["synth"])
        .arg(&cfg)
        .arg("-o")
        .arg(&synth)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let mut periods = synth.as_os_str().to_owned();
    periods.push(".periods.tsv");
    let det = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(synth.with_extension("fa"))
        .arg("--periods")
        .arg(std::path::PathBuf::from(periods))
        .arg("-o")
        .arg(&det)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    (dir, det)
}

fn read_properties_row(det_prefix: &std::path::Path) -> Vec<String> {
    let mut p = det_prefix.as_os_str().to_owned();
    p.push(".properties.tsv");
    let s = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    s.lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect()
}

#[test]
fn t10_phase_shift_offset_within_tolerance() {
    let (_dir, det) = synth_then_detect(
        r#"
schema_version: 1
seed: 10
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: 100
  - type: SHIFT
    offset_bp: 85
  - type: HOR
    template: alpha
    n_copies: 100
"#,
    );
    let row = read_properties_row(&det);
    // properties.tsv columns:
    //   12 = n_phase_shifts
    //   13 = phase_shift_positions
    //   14 = phase_shift_offsets
    //   18 = n_segments
    let n: usize = row[12].parse().unwrap();
    assert!(
        n >= 1,
        "expected ≥1 phase shift in properties; row={:?}",
        row
    );
    let offsets_str = &row[14];
    assert_ne!(offsets_str, "NA", "phase_shift_offsets should not be NA");
    let offsets: Vec<i64> = offsets_str.split(',').map(|s| s.parse().unwrap()).collect();
    // Truth: shift offset = 85 bp. Pass B reports ±85 (sign depends on
    // which direction we declare "forward"). Allow either sign and
    // ±5 bp tolerance per plan §10 M3.5 acceptance.
    let any_close = offsets.iter().any(|&o| (o.abs() - 85).abs() <= 5);
    assert!(
        any_close,
        "expected some offset ≈ ±85 bp; got {:?}",
        offsets
    );
}

#[test]
fn no_shift_yields_empty_phase_shift_lists() {
    let (_dir, det) = synth_then_detect(
        r#"
schema_version: 1
seed: 5
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: 200
"#,
    );
    let row = read_properties_row(&det);
    let n: usize = row[12].parse().unwrap();
    assert_eq!(n, 0, "clean HOR should have zero phase shifts; got {n}");
    assert_eq!(row[13], "NA");
    assert_eq!(row[14], "NA");
}

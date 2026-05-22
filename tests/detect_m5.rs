//! M5 acceptance: consensus FASTA + per-array viz export.
//!
//! - HOR-unit consensus is built from the HOR-unit width directly,
//!   NOT by repeating the base-width consensus (A9).
//! - `--viz-dir` emits a per-array subdir with the always-on TSVs
//!   (column_ic, column_edge_rate, rk, shift).
//! - With `--export-raster`, the raster TSV and (under the `viz`
//!   feature) PNG are also written.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn synth_then(yaml: &str, det_flags: &[&str]) -> (tempfile::TempDir, std::path::PathBuf) {
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
    let mut cmd = Command::new(kitehor_bin());
    cmd.arg("detect")
        .arg(synth.with_extension("fa"))
        .arg("--periods")
        .arg(std::path::PathBuf::from(periods))
        .arg("-o")
        .arg(&det);
    for f in det_flags {
        cmd.arg(f);
    }
    let o = cmd.output().unwrap();
    assert!(
        o.status.success(),
        "detect failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    (dir, det)
}

#[test]
fn consensus_fasta_emitted_for_hor() {
    let (_dir, det) = synth_then(
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
        &[],
    );
    let mut p = det.as_os_str().to_owned();
    p.push(".consensus.fa");
    let text = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    assert!(
        text.contains(">arr.monomer  length=171"),
        "consensus FASTA missing monomer record: {}",
        text.lines().next().unwrap_or("")
    );
    assert!(
        text.contains(">arr.hor_unit  length=2052  k=12"),
        "consensus FASTA missing hor_unit record"
    );
}

#[test]
fn consensus_fasta_emitted_for_simple_tr_without_hor_unit() {
    let (_dir, det) = synth_then(
        r#"
schema_version: 1
seed: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 170
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 500
"#,
        &[],
    );
    let mut p = det.as_os_str().to_owned();
    p.push(".consensus.fa");
    let text = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    assert!(text.contains(">arr.monomer  length=170"));
    assert!(
        !text.contains("hor_unit"),
        "simple_TR shouldn't have hor_unit record"
    );
}

#[test]
fn viz_dir_emits_per_array_tsvs() {
    let dir = tempfile::tempdir().unwrap();
    let viz = dir.path().join("viz");
    let yaml = r#"
schema_version: 1
seed: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 170
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 200
"#;
    let cfg = dir.path().join("cfg.yaml");
    std::fs::File::create(&cfg)
        .unwrap()
        .write_all(yaml.as_bytes())
        .unwrap();
    let synth = dir.path().join("arr");
    Command::new(kitehor_bin())
        .args(["synth"])
        .arg(&cfg)
        .arg("-o")
        .arg(&synth)
        .output()
        .unwrap();
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
        .arg("--viz-dir")
        .arg(&viz)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let arr_dir = viz.join("arr");
    let ic = arr_dir.join("column_ic_w170.tsv");
    let rk = arr_dir.join("rk_w170.tsv");
    let shift = arr_dir.join("shift_w170.tsv");
    let edge = arr_dir.join("column_edge_rate_w170.tsv");
    assert!(ic.exists(), "column_ic TSV missing");
    assert!(rk.exists(), "rk TSV missing");
    assert!(shift.exists(), "shift TSV missing");
    assert!(edge.exists(), "column_edge_rate TSV missing");
    // raster shouldn't be there without --export-raster.
    assert!(!arr_dir.join("raster_w170.tsv").exists());
}

#[test]
fn export_raster_emits_tsv_and_png() {
    let (dir, _det) = synth_then(
        r#"
schema_version: 1
seed: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 100
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 50
"#,
        &[],
    );
    let viz = dir.path().join("viz");
    let synth = dir.path().join("arr");
    let mut periods = synth.as_os_str().to_owned();
    periods.push(".periods.tsv");
    let det = dir.path().join("det2");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(synth.with_extension("fa"))
        .arg("--periods")
        .arg(std::path::PathBuf::from(periods))
        .arg("-o")
        .arg(&det)
        .arg("--viz-dir")
        .arg(&viz)
        .arg("--export-raster")
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let arr_dir = viz.join("arr");
    let tsv = arr_dir.join("raster_w100.tsv");
    let png = arr_dir.join("raster_w100.png");
    assert!(tsv.exists(), "raster TSV missing");
    assert!(
        png.exists(),
        "raster PNG missing (default viz build should produce it)"
    );
    let bytes = std::fs::read(&png).unwrap();
    assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n", "not a PNG file");
}

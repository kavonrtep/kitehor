//! M2 acceptance: oriented 4-mer row embeddings + `R(k)` populate
//! `row_lag1_similarity`, `best_lag`, `best_lag_score` per width.
//! On a clean k=12 HOR at the true base width, `best_lag` must be 12.

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

fn read_widths(prefix: &std::path::Path) -> Vec<Vec<String>> {
    let mut p = prefix.as_os_str().to_owned();
    p.push(".width_features.tsv");
    let s = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    s.lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect()
}

#[test]
fn clean_hor_k12_best_lag_is_12_at_base_width() {
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
    let rows = read_widths(&det);
    // Find the width=171 row and verify best_lag is 12 with high score.
    let r = rows
        .iter()
        .find(|r| r[1] == "171")
        .expect("width=171 missing from width_features");
    let best_lag: usize = r[6].parse().expect("best_lag not an int");
    let best_score: f64 = r[7].parse().expect("best_lag_score not f64");
    assert_eq!(
        best_lag, 12,
        "best_lag at base width should be 12; row={r:?}"
    );
    assert!(
        best_score > 0.6,
        "best_lag_score should be high at the true HOR multiplicity; got {best_score}"
    );
    let r1: f64 = r[5].parse().unwrap();
    assert!(
        best_score > r1 + 0.2,
        "best_lag_score ({best_score}) should beat R(1) ({r1}) by a clear margin at HOR base width"
    );
}

#[test]
fn simple_tr_best_lag_score_close_to_r1() {
    let (_dir, det) = synth_then_detect(
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
    );
    let rows = read_widths(&det);
    let r = rows
        .iter()
        .find(|r| r[1] == "170")
        .expect("width=170 missing");
    let r1: f64 = r[5].parse().unwrap();
    let best: f64 = r[7].parse().unwrap();
    // On a simple TR, R(1) is already very high; R(best_k) doesn't
    // exceed it by much.
    assert!(
        r1 > 0.8,
        "simple TR base width R(1) should be high; got {r1}"
    );
    assert!(
        (best - r1).abs() < 0.20,
        "simple TR: best lag score should be close to R(1); r1={r1} best={best}"
    );
}

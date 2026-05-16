//! M3 acceptance: edge field + Pass-A shift signal.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn synth_then_detect(yaml: &str) -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let cfg = dir.path().join("cfg.yaml");
    std::fs::File::create(&cfg).unwrap().write_all(yaml.as_bytes()).unwrap();
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
    s.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect()
}

#[test]
fn simple_tr_has_low_wobble_at_true_width() {
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
    let r = rows.iter().find(|r| r[1] == "170").expect("width=170 missing");
    let wobble: f64 = r[13].parse().unwrap();
    assert!(
        wobble < 1.5,
        "simple TR at true width should have low wobble; got {wobble} bp"
    );
}

#[test]
fn phase_shift_fixture_reports_at_least_one_breakpoint() {
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
    let rows = read_widths(&det);
    let r = rows.iter().find(|r| r[1] == "171").expect("width=171 missing");
    let n_phase_shifts: usize = r[14].parse().unwrap();
    assert!(
        n_phase_shifts >= 1,
        "expected ≥1 breakpoint detected at width=171; got {n_phase_shifts}"
    );
}

#[test]
fn clean_hor_has_strong_vertical_edge_rate_at_base_width() {
    // At base width = monomer length, adjacent rows are slots i ↔ i+1
    // (or i ↔ i+1 mod k at copy boundaries) — slots are diverged so
    // vertical_edge_rate is appreciable (NOT zero like at HOR-unit
    // width, where adjacent rows are the same unit).
    let (_dir, det) = synth_then_detect(
        r#"
schema_version: 1
seed: 5
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    inter_slot_divergence: 0.20
structure:
  - type: HOR
    template: alpha
    n_copies: 200
"#,
    );
    let rows = read_widths(&det);
    let r_base = rows.iter().find(|r| r[1] == "171").expect("width=171 missing");
    let v_base: f64 = r_base[9].parse().unwrap();
    assert!(
        v_base > 0.10,
        "HOR base width should have non-trivial vertical_edge_rate; got {v_base}"
    );
    // At HOR-unit width 2052 the same metric should be much lower
    // because adjacent rows are the same logical HOR unit.
    if let Some(r_unit) = rows.iter().find(|r| r[1] == "2052") {
        let v_unit: f64 = r_unit[9].parse().unwrap_or(1.0);
        assert!(
            v_unit < v_base,
            "vertical_edge_rate at HOR-unit width ({v_unit}) should be < base width ({v_base})"
        );
    }
}

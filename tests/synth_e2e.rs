//! M4 end-to-end: `kitehor synth <yaml> -o PREFIX` produces a valid
//! FASTA + truth + periods bundle.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn run_synth(yaml: &str, prefix: &std::path::Path) -> std::process::Output {
    let mut cfg = tempfile::NamedTempFile::new().unwrap();
    cfg.write_all(yaml.as_bytes()).unwrap();
    Command::new(kitehor_bin())
        .arg("synth")
        .arg(cfg.path())
        .arg("-o")
        .arg(prefix)
        .output()
        .expect("run synth")
}

fn read_seq(prefix: &std::path::Path) -> Vec<u8> {
    let fa = std::fs::read_to_string(prefix.with_extension("fa")).unwrap();
    let mut seq = Vec::new();
    for line in fa.lines() {
        if line.starts_with('>') {
            continue;
        }
        seq.extend_from_slice(line.as_bytes());
    }
    seq
}

fn read_truth(prefix: &std::path::Path) -> (Vec<String>, Vec<String>) {
    let mut p = prefix.as_os_str().to_owned();
    p.push(".truth.tsv");
    let t = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    let mut lines = t.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(str::to_string).collect();
    let row: Vec<String> = lines.next().unwrap().split('\t').map(str::to_string).collect();
    (header, row)
}

fn read_periods(prefix: &std::path::Path) -> Vec<Vec<String>> {
    let mut p = prefix.as_os_str().to_owned();
    p.push(".periods.tsv");
    let t = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    let mut rows = Vec::new();
    for (i, line) in t.lines().enumerate() {
        if i == 0 {
            continue; // header
        }
        rows.push(line.split('\t').map(str::to_string).collect());
    }
    rows
}

#[test]
fn t05_clean_hor_runs() {
    let yaml = r#"
schema_version: 1
seed: 5
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: 50
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t05");
    let out = run_synth(yaml, &prefix);
    assert!(
        out.status.success(),
        "synth failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 50);

    let (header, row) = read_truth(&prefix);
    let idx_class = header.iter().position(|h| h == "truth_class").unwrap();
    let idx_base = header.iter().position(|h| h == "base_width_bp").unwrap();
    let idx_k = header.iter().position(|h| h == "hor_k").unwrap();
    let idx_n = header.iter().position(|h| h == "n_complete_copies").unwrap();
    let idx_se = header
        .iter()
        .position(|h| h == "structural_expression")
        .unwrap();
    assert_eq!(row[idx_class], "HOR");
    assert_eq!(row[idx_base], "100");
    assert_eq!(row[idx_k], "4");
    assert_eq!(row[idx_n], "50");
    assert_eq!(row[idx_se], "H([M_1..M_4],50,div=0.15)");

    let periods = read_periods(&prefix);
    let bases: Vec<&str> = periods.iter().map(|r| r[3].as_str()).collect();
    assert!(bases.contains(&"true_base"));
    assert!(bases.contains(&"true_hor_unit"));
}

#[test]
fn t01_simple_tr_runs() {
    let yaml = r#"
schema_version: 1
seed: 1
global:
  mutation_rate: 0.02
templates:
  m:
    type: monomer
    monomer_length_bp: 171
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 200
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t01");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    // 171 * 200 + a small number of indels (= 0 here since indel_rate=0)
    assert_eq!(seq.len(), 171 * 200);
    let (header, row) = read_truth(&prefix);
    let idx_class = header.iter().position(|h| h == "truth_class").unwrap();
    assert_eq!(row[idx_class], "simple_TR");
}

#[test]
fn t10_phase_shift_runs() {
    let yaml = r#"
schema_version: 1
seed: 10
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: 50
  - type: SHIFT
    offset_bp: 25
  - type: HOR
    template: alpha
    n_copies: 50
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t10");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 50 + 25 + 100 * 4 * 50);
    let (header, row) = read_truth(&prefix);
    let idx_class = header.iter().position(|h| h == "truth_class").unwrap();
    let idx_shifts = header.iter().position(|h| h == "n_phase_shifts").unwrap();
    let idx_pos = header
        .iter()
        .position(|h| h == "phase_shift_positions")
        .unwrap();
    let idx_off = header
        .iter()
        .position(|h| h == "phase_shift_offsets")
        .unwrap();
    let idx_seg = header.iter().position(|h| h == "n_segments").unwrap();
    assert_eq!(row[idx_class], "HOR");
    assert_eq!(row[idx_shifts], "1");
    assert_eq!(row[idx_pos], "20000");
    assert_eq!(row[idx_off], "25");
    assert_eq!(row[idx_seg], "2");
}

#[test]
fn t11_insertion_runs() {
    let yaml = r#"
schema_version: 1
seed: 11
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: 30
  - type: INSERTION
    length_bp: 500
    kind: retro_like
  - type: HOR
    template: alpha
    n_copies: 30
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t11");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 30 + 500 + 100 * 4 * 30);
}

#[test]
fn determinism_same_seed_byte_identical() {
    let yaml = r#"
schema_version: 1
seed: 42
global:
  mutation_rate: 0.03
  indel_rate: 0.01
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.10
structure:
  - type: HOR
    template: t
    n_copies: 50
"#;
    let dir = tempfile::tempdir().unwrap();
    let p1 = dir.path().join("a");
    let p2 = dir.path().join("b");
    assert!(run_synth(yaml, &p1).status.success());
    assert!(run_synth(yaml, &p2).status.success());
    let s1 = read_seq(&p1);
    let s2 = read_seq(&p2);
    assert_eq!(s1, s2, "same seed must produce byte-identical FASTA");
}

#[test]
fn t03_wobble_aperiodic_runs() {
    let yaml = r#"
schema_version: 1
seed: 3
templates:
  m:
    type: monomer
    monomer_length_bp: 170
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 1000
modifiers:
  - wobble:
      amplitude_bp: 2.0
      model: random_walk
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t03");
    let out = run_synth(yaml, &prefix);
    assert!(
        out.status.success(),
        "synth failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (header, row) = read_truth(&prefix);
    let idx_amp = header
        .iter()
        .position(|h| h == "wobble_amplitude_bp")
        .unwrap();
    let idx_per = header
        .iter()
        .position(|h| h == "wobble_periodicity_bp")
        .unwrap();
    let amp: f64 = row[idx_amp].parse().unwrap();
    assert!(
        amp > 0.5 && amp < 4.0,
        "expected realised wobble amplitude ~2 bp, got {amp}"
    );
    assert_eq!(row[idx_per], "NA");
}

#[test]
fn t04_wobble_periodic_runs() {
    let yaml = r#"
schema_version: 1
seed: 4
templates:
  m:
    type: monomer
    monomer_length_bp: 100
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 2000
modifiers:
  - wobble:
      amplitude_bp: 1.5
      period_rows: 500
      model: sinusoidal
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t04");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let (header, row) = read_truth(&prefix);
    let idx_per = header
        .iter()
        .position(|h| h == "wobble_periodicity_bp")
        .unwrap();
    // 500 rows × 100 bp = 50_000 bp.
    assert!(row[idx_per].starts_with("50000"), "got {}", row[idx_per]);
}

#[test]
fn t16_hybrid_runs_and_events_json_populated() {
    let yaml = r#"
schema_version: 1
seed: 16
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.20
structure:
  - type: HOR
    template: t
    n_copies: 100
post_generation:
  - type: HYBRID
    block: 0
    at_copy: 50
    slot: 3
    source_slots: [3, 4]
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t16");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 100);
    let (header, row) = read_truth(&prefix);
    let idx = header.iter().position(|h| h == "events_json").unwrap();
    let json: serde_json::Value = serde_json::from_str(&row[idx]).unwrap();
    let arr = json.as_array().unwrap();
    assert_eq!(arr.len(), 1);
    assert_eq!(arr[0]["type"], "HYBRID");
    assert_eq!(arr[0]["copy"], 50);
    assert_eq!(arr[0]["slot"], 3);
}

#[test]
fn t12_inversion_preserves_length_and_logs_event() {
    let yaml = r#"
schema_version: 1
seed: 12
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 210
post_generation:
  - type: INVERSION
    block: 0
    start_copy: 101
    length_copies: 10
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t12");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 210);
    let (header, row) = read_truth(&prefix);
    let idx = header.iter().position(|h| h == "events_json").unwrap();
    let json: serde_json::Value = serde_json::from_str(&row[idx]).unwrap();
    let arr = json.as_array().unwrap();
    assert_eq!(arr[0]["type"], "INVERSION");
    assert_eq!(arr[0]["length_bp"], 4000);
    assert_eq!(arr[0]["start_bp"], 100 * 4 * 100);
}

#[test]
fn duplication_extends_array() {
    let yaml = r#"
schema_version: 1
seed: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 50
post_generation:
  - type: DUPLICATION
    block: 0
    start_copy: 10
    length_copies: 5
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("dup");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 50 + 100 * 4 * 5);
}

#[test]
fn deletion_shrinks_array() {
    let yaml = r#"
schema_version: 1
seed: 1
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: 50
post_generation:
  - type: DELETION
    block: 0
    start_copy: 10
    length_copies: 5
"#;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("del");
    let out = run_synth(yaml, &prefix);
    assert!(out.status.success());
    let seq = read_seq(&prefix);
    assert_eq!(seq.len(), 100 * 4 * 50 - 100 * 4 * 5);
}

#[test]
fn seed_cli_override_changes_output() {
    let yaml = r#"
schema_version: 1
seed: 42
templates:
  t:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    inter_slot_divergence: 0.10
structure:
  - type: HOR
    template: t
    n_copies: 50
"#;
    let mut cfg = tempfile::NamedTempFile::new().unwrap();
    cfg.write_all(yaml.as_bytes()).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let p1 = dir.path().join("a");
    let p2 = dir.path().join("b");
    let _ = Command::new(kitehor_bin())
        .arg("synth")
        .arg(cfg.path())
        .arg("-o")
        .arg(&p1)
        .arg("--seed")
        .arg("42")
        .output()
        .unwrap();
    let _ = Command::new(kitehor_bin())
        .arg("synth")
        .arg(cfg.path())
        .arg("-o")
        .arg(&p2)
        .arg("--seed")
        .arg("123")
        .output()
        .unwrap();
    let s1 = read_seq(&p1);
    let s2 = read_seq(&p2);
    assert_ne!(s1, s2);
}

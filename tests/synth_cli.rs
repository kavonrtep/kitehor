//! M1 acceptance: `kitehor synth-validate` / `synth-schema` round-trip.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    // Cargo sets CARGO_BIN_EXE_<name> for integration tests.
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

#[test]
fn synth_schema_prints_canonical_schema_verbatim() {
    let out = Command::new(kitehor_bin())
        .arg("synth-schema")
        .output()
        .expect("run synth-schema");
    assert!(
        out.status.success(),
        "synth-schema failed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let canonical = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("docs")
            .join("new")
            .join("simulator_schema.json"),
    )
    .expect("read canonical schema");
    let printed = String::from_utf8_lossy(&out.stdout).into_owned();
    assert_eq!(printed, canonical, "synth-schema stdout drift");
}

#[test]
fn synth_validate_accepts_minimal_hor() {
    let yaml = r#"
schema_version: 1
seed: 7
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
"#;
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(yaml.as_bytes()).unwrap();
    let out = Command::new(kitehor_bin())
        .args(["synth-validate"])
        .arg(f.path())
        .output()
        .expect("run synth-validate");
    assert!(
        out.status.success(),
        "expected validate to succeed; stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn synth_validate_rejects_source_file() {
    let yaml = r#"
schema_version: 1
templates:
  bad:
    type: HOR_slots
    monomer_length_bp: 100
    k: 4
    source: file
    file: /nope.fa
structure:
  - type: HOR
    template: bad
    n_copies: 5
"#;
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(yaml.as_bytes()).unwrap();
    let out = Command::new(kitehor_bin())
        .args(["synth-validate"])
        .arg(f.path())
        .output()
        .expect("run synth-validate");
    assert!(!out.status.success(), "expected validate to fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("source: file is not implemented in MVP"),
        "missing MVP message; stderr={stderr}"
    );
}

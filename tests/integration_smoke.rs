//! End-to-end smoke test: invoke the built binary on the shipped
//! synthetic fixture under `test_data/smoke/` and verify the verdicts.
//!
//! The fixture has 3 records:
//!   - `tandem_pure` (60 bp × 300 copies)        → tandem
//!   - `hor_k3`      (100 bp founder, k=3)       → hor, founder=100, tile=300
//!   - `hor_k5`      (150 bp founder, k=5)       → hor, founder=150, tile=750
//!
//! If this test breaks, either the model artefacts in `models/` were
//! changed, or the kite/feature/classifier code drifted.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn binary_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("target");
    p.push(if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    });
    p.push("kitehor");
    p
}

fn parse_tsv(text: &str) -> Vec<HashMap<String, String>> {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            header
                .iter()
                .zip(cols.iter())
                .map(|(h, c)| ((*h).to_string(), (*c).to_string()))
                .collect()
        })
        .collect()
}

#[test]
#[ignore = "requires `cargo build --release` first; run with --ignored"]
fn smoke_classifier_verdicts_match_truth() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fasta = manifest.join("test_data/smoke/sequences.fasta");
    let tmp = std::env::temp_dir().join("kitehor_smoke_out.tsv");
    let _ = std::fs::remove_file(&tmp);

    let bin = binary_path();
    assert!(
        bin.exists(),
        "binary {:?} not found — run `cargo build --release` first",
        bin
    );

    let status = Command::new(&bin)
        .args([
            "kite-periodicity",
            fasta.to_str().unwrap(),
            "-o",
            tmp.to_str().unwrap(),
            "--classify",
            "--no-hor-call",
            "--threads",
            "2",
        ])
        .status()
        .expect("spawning kitehor failed");
    assert!(status.success(), "kitehor exited non-zero: {status}");

    let out = std::fs::read_to_string(&tmp).expect("read output");
    let rows = parse_tsv(&out);
    assert_eq!(rows.len(), 3, "expected 3 records, got {}", rows.len());

    let by_id: HashMap<String, &HashMap<String, String>> =
        rows.iter().map(|r| (r["case_id"].clone(), r)).collect();

    let tandem = by_id.get("tandem_pure").expect("tandem_pure missing");
    assert_eq!(tandem["verdict"], "tandem", "tandem_pure verdict mismatch");

    let hor3 = by_id.get("hor_k3").expect("hor_k3 missing");
    assert_eq!(hor3["verdict"], "hor", "hor_k3 verdict");
    assert_eq!(hor3["founder"], "100", "hor_k3 founder");
    assert_eq!(hor3["multiplicity"], "3", "hor_k3 k");
    assert_eq!(hor3["tile"], "300", "hor_k3 tile");

    let hor5 = by_id.get("hor_k5").expect("hor_k5 missing");
    assert_eq!(hor5["verdict"], "hor", "hor_k5 verdict");
    assert_eq!(hor5["founder"], "150", "hor_k5 founder");
    assert_eq!(hor5["multiplicity"], "5", "hor_k5 k");
    assert_eq!(hor5["tile"], "750", "hor_k5 tile");
}

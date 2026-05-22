//! End-to-end orchestrator smoke test. Runs `kitehor analyze` on the
//! 3-record smoke FASTA and asserts all 9 per-stage TSVs are emitted
//! and the combined_class for each record matches expectations.

use std::path::PathBuf;
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

#[test]
fn analyze_emits_all_per_stage_tsvs() {
    let fasta =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("test_data/smoke/sequences.fasta");
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("pipe");
    let ok = Command::new(kitehor_bin())
        .args(["analyze"])
        .arg(&fasta)
        .arg("-o")
        .arg(&prefix)
        .status()
        .unwrap();
    assert!(ok.success(), "analyze failed");

    let prefix_s = prefix.to_str().unwrap();
    for suffix in [
        "kite.tsv",
        "kite.peaks.tsv",
        "verdicts.tsv",
        "subrepeat.tsv",
        "windows.tsv",
        "ssr.tsv",
        "ssr.regions.tsv",
        "hor_within_tile.tsv",
        "summary.tsv",
    ] {
        let p = format!("{prefix_s}.{suffix}");
        assert!(
            std::path::Path::new(&p).exists(),
            "missing expected output {p}"
        );
    }

    // summary.tsv: expect 3 records, hor_k3/k5 → hor, tandem_pure → tr.
    let body = std::fs::read_to_string(format!("{prefix_s}.summary.tsv")).unwrap();
    let mut by_id: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    for line in body.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        by_id.insert(cells[0], *cells.last().unwrap());
    }
    assert_eq!(by_id.get("hor_k3"), Some(&"hor"));
    assert_eq!(by_id.get("hor_k5"), Some(&"hor"));
    assert_eq!(by_id.get("tandem_pure"), Some(&"tr"));
}

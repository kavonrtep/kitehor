//! CI corpus regression: run `kitehor analyze` against the curated
//! `test_data/ci_corpus/sequences.fasta` and check the per-record
//! `combined_class` against the documented baseline in
//! `test_data/ci_corpus/manifest.tsv`.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn parse_manifest() -> HashMap<String, String> {
    let path = root().join("test_data/ci_corpus/manifest.tsv");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading manifest {path:?}: {e}"));
    let mut out = HashMap::new();
    for line in text.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        out.insert(cells[0].to_string(), cells[2].to_string());
    }
    out
}

#[test]
fn analyze_ci_corpus_matches_manifest() {
    let fasta = root().join("test_data/ci_corpus/sequences.fasta");
    if !fasta.exists() {
        eprintln!("skipping: ci_corpus/sequences.fasta not present");
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("ci");
    let ok = Command::new(kitehor_bin())
        .args(["analyze"])
        .arg(&fasta)
        .arg("-o")
        .arg(&prefix)
        .status()
        .unwrap();
    assert!(ok.success(), "analyze failed");

    let mut prefix_s = prefix.into_os_string();
    prefix_s.push(".summary.tsv");
    let summary = PathBuf::from(prefix_s);
    let body = std::fs::read_to_string(&summary).unwrap();
    let mut got: HashMap<String, String> = HashMap::new();
    for line in body.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        got.insert(cells[0].to_string(), cells.last().unwrap().to_string());
    }

    let expected = parse_manifest();
    let mut diffs: Vec<String> = Vec::new();
    for (rec, want) in &expected {
        match got.get(rec) {
            Some(actual) if actual == want => {}
            Some(actual) => diffs.push(format!(
                "{rec}: combined_class = {actual} (expected {want})"
            )),
            None => diffs.push(format!("{rec}: missing from summary.tsv")),
        }
    }
    assert!(
        diffs.is_empty(),
        "ci_corpus regression:\n  {}",
        diffs.join("\n  ")
    );
}

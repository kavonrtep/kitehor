//! M7 acceptance: `kitehor synth-batch` over the staged corpus.
//!
//! Walks every `*.yaml` under `tests/synth_configs/` (skipping
//! `.deferred.yaml`), runs in parallel via the binary, and verifies:
//!
//! 1. All 22 active fixtures produce non-empty `.fa` + `.truth.tsv` +
//!    `.periods.tsv` outputs.
//! 2. The deferred placeholder is **not** generated.
//! 3. Re-running the same `--seed-offset` against the same corpus
//!    yields byte-identical `.fa` outputs (determinism).
//! 4. The full batch finishes in well under 30 s on the container.

use std::process::Command;
use std::time::Instant;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn corpus_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("synth_configs")
}

fn count_yaml(dir: &std::path::Path, deferred: bool) -> usize {
    std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name().into_string().unwrap_or_default();
            let yaml = name.ends_with(".yaml");
            let is_deferred = name.ends_with(".deferred.yaml");
            yaml && (is_deferred == deferred)
        })
        .count()
}

#[test]
fn corpus_has_22_active_and_1_deferred() {
    let active = count_yaml(&corpus_dir(), false);
    let deferred = count_yaml(&corpus_dir(), true);
    assert_eq!(
        active, 22,
        "expected 22 active fixtures, found {active}"
    );
    assert_eq!(
        deferred, 1,
        "expected 1 deferred placeholder, found {deferred}"
    );
}

#[test]
fn batch_produces_all_active_bundles() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(kitehor_bin())
        .arg("synth-batch")
        .arg("--config-dir")
        .arg(corpus_dir())
        .arg("--out-dir")
        .arg(dir.path())
        .output()
        .expect("run synth-batch");
    assert!(
        out.status.success(),
        "synth-batch failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    // Each active config should have produced three files.
    let expected = 22;
    let fa_count = count_ext(dir.path(), "fa");
    let tr_count = count_ext_chain(dir.path(), "truth.tsv");
    let pe_count = count_ext_chain(dir.path(), "periods.tsv");
    assert_eq!(fa_count, expected, ".fa count mismatch");
    assert_eq!(tr_count, expected, ".truth.tsv count mismatch");
    assert_eq!(pe_count, expected, ".periods.tsv count mismatch");

    // Deferred placeholder must not have produced any output.
    let deferred_fa = dir.path().join("T09_nested_hor.deferred.fa");
    assert!(
        !deferred_fa.exists(),
        ".deferred.yaml fixtures must be skipped"
    );
}

#[test]
fn batch_is_deterministic_across_runs() {
    let d1 = tempfile::tempdir().unwrap();
    let d2 = tempfile::tempdir().unwrap();
    for d in [&d1, &d2] {
        let out = Command::new(kitehor_bin())
            .arg("synth-batch")
            .arg("--config-dir")
            .arg(corpus_dir())
            .arg("--out-dir")
            .arg(d.path())
            .output()
            .unwrap();
        assert!(out.status.success());
    }
    // Compare every produced .fa file byte-for-byte.
    for entry in std::fs::read_dir(d1.path()).unwrap() {
        let e = entry.unwrap();
        let name = e.file_name();
        if !name.to_string_lossy().ends_with(".fa") {
            continue;
        }
        let p1 = e.path();
        let p2 = d2.path().join(&name);
        let s1 = std::fs::read(&p1).unwrap();
        let s2 = std::fs::read(&p2).unwrap();
        assert_eq!(s1, s2, "drift in {}", name.to_string_lossy());
    }
}

#[test]
fn batch_finishes_well_under_30s() {
    let dir = tempfile::tempdir().unwrap();
    let t = Instant::now();
    let out = Command::new(kitehor_bin())
        .arg("synth-batch")
        .arg("--config-dir")
        .arg(corpus_dir())
        .arg("--out-dir")
        .arg(dir.path())
        .output()
        .unwrap();
    assert!(out.status.success());
    let elapsed = t.elapsed();
    assert!(
        elapsed.as_secs() < 30,
        "synth-batch took {elapsed:?}, expected <30s"
    );
}

#[test]
fn diagnostics_emitted_when_flag_set() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(kitehor_bin())
        .arg("synth-batch")
        .arg("--config-dir")
        .arg(corpus_dir())
        .arg("--out-dir")
        .arg(dir.path())
        .arg("--diagnostics")
        .output()
        .unwrap();
    assert!(out.status.success());
    let n = count_ext_chain(dir.path(), "diagnostics.json");
    assert_eq!(n, 22, "diagnostics file count mismatch");

    // Spot-check structure of one diagnostics file.
    let any = std::fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .find(|e| e.file_name().to_string_lossy().ends_with(".diagnostics.json"))
        .unwrap()
        .path();
    let j: serde_json::Value = serde_json::from_slice(&std::fs::read(&any).unwrap()).unwrap();
    assert!(j["rng_seeds"]["top"].is_number());
    assert!(j["templates"].is_object());
    assert!(j["blocks"].is_array());
    assert!(j["noise"]["n_substitutions"].is_number());
    assert!(j["sequence_length_bp"].is_number());
}

// ---- Helpers ----

fn count_ext(dir: &std::path::Path, ext: &str) -> usize {
    std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .and_then(|x| x.to_str())
                .map(|x| x == ext)
                .unwrap_or(false)
        })
        .count()
}

/// Match files whose name ends with `.<chain>` (e.g. `.truth.tsv`).
/// `Path::extension` only captures the segment after the final dot, so
/// it doesn't work for `truth.tsv` etc.
fn count_ext_chain(dir: &std::path::Path, chain: &str) -> usize {
    let suffix = format!(".{}", chain);
    std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().ends_with(&suffix))
        .count()
}

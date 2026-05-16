//! End-to-end integration test for the kite → detect pipeline
//! (`docs/reviews/kite_emit_periods_integration_review_2026-05-16.md`
//! review finding #4).
//!
//! Verifies:
//!   1. `kite-periodicity --classify --emit-periods` succeeds and
//!      writes a v2-schema `periods.tsv` for every record kite
//!      analysed.
//!   2. `detect --periods` accepts that file without error and
//!      writes a schema-valid `properties.tsv` with one row per
//!      FASTA record.
//!   3. `--use-ml-classifier` + `--emit-periods` is rejected at
//!      CLI parse time (clap conflict).
//!
//! Uses the committed `test_data/smoke/sequences.fasta` fixture
//! (3 records: tandem_pure, hor_k3, hor_k5) so the test runs in
//! milliseconds and doesn't depend on regenerated corpora.

use std::path::{Path, PathBuf};
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn smoke_fasta() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_data")
        .join("smoke")
        .join("sequences.fasta")
}

#[test]
fn kite_emit_periods_then_detect_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let kite_tsv = dir.path().join("smoke.kite.tsv");
    let periods_tsv = dir.path().join("smoke.kite.periods.tsv");

    // Stage 1: kite-periodicity --classify --emit-periods.
    let o = Command::new(kitehor_bin())
        .arg("kite-periodicity")
        .arg(smoke_fasta())
        .arg("-o")
        .arg(&kite_tsv)
        .arg("--classify")
        .arg("--emit-periods")
        .arg(&periods_tsv)
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "kite-periodicity failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(periods_tsv.exists(), "--emit-periods did not write the file");

    let periods = std::fs::read_to_string(&periods_tsv).unwrap();
    let mut lines = periods.lines();
    let header = lines.next().expect("periods.tsv is empty");
    assert_eq!(
        header, "array_id\tperiod_bp\tperiod_score\tsource",
        "periods.tsv header doesn't match v2 schema"
    );
    let data_lines: Vec<&str> = lines.collect();
    assert!(
        !data_lines.is_empty(),
        "smoke fixture must produce at least one periods row"
    );

    // Each data row must have 4 tab-separated fields, period_score
    // parseable and within [0, 1], and an array_id present in the
    // FASTA. The detector loader validates strictly, so this mirrors
    // its parser.
    let known_ids: std::collections::HashSet<&str> =
        ["tandem_pure", "hor_k3", "hor_k5"].into_iter().collect();
    let mut high_score_rows = 0usize;
    for line in &data_lines {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 4, "row has wrong column count: {:?}", line);
        let array_id = f[0];
        assert!(
            known_ids.contains(array_id),
            "unknown array_id `{array_id}` in periods row"
        );
        let _period: usize = f[1].parse().expect("period_bp not integer");
        let score: f64 = f[2].parse().expect("period_score not float");
        assert!(
            (0.0..=1.0).contains(&score) && score.is_finite(),
            "period_score {score} out of [0, 1]"
        );
        if score >= 0.85 {
            high_score_rows += 1;
        }
    }
    // The HOR fixtures (hor_k3, hor_k5) should each produce at least
    // one high-score row (founder=0.95). Combined ≥ 2.
    assert!(
        high_score_rows >= 2,
        "expected ≥ 2 high-score (founder/tile) rows; got {high_score_rows}"
    );

    // Stage 2: detect --periods <emitted>.
    let det_prefix = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(smoke_fasta())
        .arg("--periods")
        .arg(&periods_tsv)
        .arg("-o")
        .arg(&det_prefix)
        .arg("--allow-missing-periods")
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "detect failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );

    let props_path = {
        let mut p = det_prefix.into_os_string();
        p.push(".properties.tsv");
        PathBuf::from(p)
    };
    let props = std::fs::read_to_string(&props_path).unwrap();
    let prop_lines: Vec<&str> = props.lines().collect();
    // Header + one row per FASTA record (3 records in the fixture).
    assert_eq!(
        prop_lines.len(),
        4,
        "expected 1 header + 3 data rows in properties.tsv; got {} lines",
        prop_lines.len()
    );
    // Header column count must match the frozen schema (20 columns).
    let header_cols: Vec<&str> = prop_lines[0].split('\t').collect();
    assert_eq!(header_cols.len(), 20, "properties.tsv header column count drift");
    // Every data row matches header column count.
    for row in &prop_lines[1..] {
        let f: Vec<&str> = row.split('\t').collect();
        assert_eq!(
            f.len(), header_cols.len(),
            "row column count mismatch: {row:?}"
        );
    }
    // At least one record should produce a non-ambiguous class.
    // (The smoke fixture contains valid HOR + tandem arrays.)
    let class_col = header_cols.iter().position(|h| *h == "class").unwrap();
    let n_resolved = prop_lines[1..]
        .iter()
        .filter(|row| {
            let f: Vec<&str> = row.split('\t').collect();
            matches!(f[class_col], "HOR" | "simple_TR" | "irregular_HOR")
        })
        .count();
    assert!(
        n_resolved >= 1,
        "expected ≥ 1 resolved class on smoke fixture; got {n_resolved}"
    );
}

#[test]
fn emit_periods_conflicts_with_use_ml_classifier() {
    let dir = tempfile::tempdir().unwrap();
    let kite_tsv = dir.path().join("k.tsv");
    let periods_tsv = dir.path().join("k.periods.tsv");

    // Clap should reject this combination at parse time (Review-#2).
    let o = Command::new(kitehor_bin())
        .arg("kite-periodicity")
        .arg(smoke_fasta())
        .arg("-o")
        .arg(&kite_tsv)
        .arg("--classify")
        .arg("--use-ml-classifier")
        .arg("--emit-periods")
        .arg(&periods_tsv)
        .output()
        .unwrap();
    assert!(
        !o.status.success(),
        "--emit-periods + --use-ml-classifier should be rejected by clap"
    );
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(
        stderr.contains("conflicts") || stderr.contains("cannot be used with"),
        "expected clap conflict message; got: {stderr}"
    );
}

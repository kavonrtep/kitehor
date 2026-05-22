//! Regression: the new `kitehor rule-classify` Rust port must match the
//! Python prototype's `tools/rule_proto/test_rule_proto.py` semantics on
//! the 6 hand-curated fixtures in `tools/rule_proto/fixtures/`.
//!
//! The harness is band-tolerant (per `tools/rule_proto/expected.tsv`):
//! - `verdict` must match exactly.
//! - `founder ∈ [founder_min, founder_max]`.
//! - For HOR: `int(multiplicity) == int(expected.multiplicity)` and
//!   `tile ∈ [tile_min, tile_max]`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    // Tests live in `<repo>/tests/`; binary in `<repo>/target/release/`.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

#[derive(Debug, Clone)]
struct ExpectedRow {
    case_id: String,
    verdict: String,
    founder_min: f64,
    founder_max: f64,
    multiplicity: Option<u32>,
    tile_min: Option<f64>,
    tile_max: Option<f64>,
}

fn parse_expected(path: &Path) -> Vec<ExpectedRow> {
    let text = std::fs::read_to_string(path).expect("reading expected.tsv");
    let mut lines = text.lines();
    let _header = lines.next().expect("expected.tsv missing header");
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        rows.push(ExpectedRow {
            case_id: f[0].into(),
            verdict: f[1].into(),
            founder_min: f[2].parse().unwrap(),
            founder_max: f[3].parse().unwrap(),
            multiplicity: f[4].parse().ok(),
            tile_min: f[5].parse().ok(),
            tile_max: f[6].parse().ok(),
        });
    }
    rows
}

#[derive(Debug, Clone)]
struct PredictedRow {
    verdict: String,
    founder: Option<f64>,
    multiplicity: Option<u32>,
    tile: Option<f64>,
}

fn parse_verdicts(path: &Path) -> Vec<PredictedRow> {
    let text = std::fs::read_to_string(path).expect("reading verdicts.tsv");
    let mut lines = text.lines();
    let _header = lines.next().expect("verdicts.tsv missing header");
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        rows.push(PredictedRow {
            verdict: f[1].into(),
            founder: parse_opt(f[2]),
            multiplicity: parse_opt_int(f[3]),
            tile: parse_opt(f[4]),
        });
    }
    rows
}

fn parse_opt(s: &str) -> Option<f64> {
    if s.is_empty() {
        None
    } else {
        s.parse().ok()
    }
}

fn parse_opt_int(s: &str) -> Option<u32> {
    if s.is_empty() {
        None
    } else {
        s.parse().ok()
    }
}

#[test]
fn fixtures_match_expected() {
    let root = repo_root();
    let fixtures_dir = root.join("tools/rule_proto/fixtures");
    let expected_path = root.join("tools/rule_proto/expected.tsv");
    let expected = parse_expected(&expected_path);
    assert!(!expected.is_empty(), "expected.tsv had no rows");

    let tmp = tempfile::tempdir().expect("creating tempdir");
    let mut failures: Vec<String> = Vec::new();

    for exp in &expected {
        let peaks = fixtures_dir.join(format!("{}.peaks.tsv", exp.case_id));
        assert!(
            peaks.exists(),
            "fixture not found: {peaks:?} — check tools/rule_proto/fixtures/"
        );
        let out_prefix = tmp.path().join(format!("{}.out", exp.case_id));
        let status = Command::new(kitehor_bin())
            .arg("rule-classify")
            .arg(&peaks)
            .arg("-o")
            .arg(&out_prefix)
            .status()
            .expect("running kitehor rule-classify");
        assert!(status.success(), "rule-classify failed for {}", exp.case_id);

        let verdicts_path = {
            let mut p = out_prefix.into_os_string();
            p.push(".verdicts.tsv");
            PathBuf::from(p)
        };
        let predicted = parse_verdicts(&verdicts_path);
        assert_eq!(
            predicted.len(),
            1,
            "{}: expected 1 verdict, got {}",
            exp.case_id,
            predicted.len()
        );
        let pred = &predicted[0];

        if pred.verdict != exp.verdict {
            failures.push(format!(
                "{}: verdict={} (expected {})",
                exp.case_id, pred.verdict, exp.verdict
            ));
            continue;
        }
        // For both hor and simple_tr the founder must be inside the band.
        let Some(f_val) = pred.founder else {
            failures.push(format!("{}: founder is empty", exp.case_id));
            continue;
        };
        if !(exp.founder_min <= f_val && f_val <= exp.founder_max) {
            failures.push(format!(
                "{}: founder={} outside [{}, {}]",
                exp.case_id, f_val, exp.founder_min, exp.founder_max
            ));
            continue;
        }
        if exp.verdict == "hor" {
            let Some(want_k) = exp.multiplicity else {
                panic!("hor row missing multiplicity in expected.tsv")
            };
            let Some(got_k) = pred.multiplicity else {
                failures.push(format!("{}: multiplicity empty", exp.case_id));
                continue;
            };
            if want_k != got_k {
                failures.push(format!(
                    "{}: multiplicity={} (expected {})",
                    exp.case_id, got_k, want_k
                ));
                continue;
            }
            let Some(t_val) = pred.tile else {
                failures.push(format!("{}: tile empty", exp.case_id));
                continue;
            };
            let lo = exp.tile_min.unwrap();
            let hi = exp.tile_max.unwrap();
            if !(lo <= t_val && t_val <= hi) {
                failures.push(format!(
                    "{}: tile={} outside [{}, {}]",
                    exp.case_id, t_val, lo, hi
                ));
            }
        }
        eprintln!("[PASS] {}", exp.case_id);
    }

    assert!(
        failures.is_empty(),
        "rule-classify fixture regressions:\n  {}",
        failures.join("\n  ")
    );
}

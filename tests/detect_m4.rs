//! M4 acceptance: detector classifies the **8 core CI fixtures**
//! against `tests/detect_expectations.tsv` per plan §10 A8.
//!
//! Core fixtures (must pass exactly):
//!   T01_simple_tr, T05_hor_clean, T06_regime_A, T07_regime_C,
//!   T10_phase_shift, T13_coexisting_periods, T17_random_negative,
//!   T18_at_rich.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("synth_configs")
}

#[derive(Debug, Clone)]
struct Expectation {
    array_id: String,
    expected_class: String,
    expected_base_width_bp: Option<usize>,
    expected_hor_k: Option<usize>,
    base_width_tol_bp: usize,
    expected_reason_contains: String,
}

fn load_expectations() -> HashMap<String, Expectation> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("detect_expectations.tsv");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut out = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let id = f[0].to_string();
        let parse_opt_usize = |s: &str| -> Option<usize> {
            if s == "NA" {
                None
            } else {
                s.parse().ok()
            }
        };
        out.insert(
            id.clone(),
            Expectation {
                array_id: id,
                expected_class: f[1].to_string(),
                expected_base_width_bp: parse_opt_usize(f[2]),
                expected_hor_k: parse_opt_usize(f[3]),
                base_width_tol_bp: f[4].parse().unwrap_or(0),
                expected_reason_contains: f.get(5).map(|s| s.to_string()).unwrap_or_default(),
            },
        );
    }
    out
}

fn synth_and_detect(corpus_yaml: &Path) -> (tempfile::TempDir, PathBuf) {
    let stem = corpus_yaml.file_stem().unwrap().to_string_lossy().to_string();
    let dir = tempfile::tempdir().unwrap();
    let synth = dir.path().join(&stem);
    let o = Command::new(kitehor_bin())
        .args(["synth"])
        .arg(corpus_yaml)
        .arg("-o")
        .arg(&synth)
        .output()
        .unwrap();
    assert!(o.status.success(), "synth failed for {stem}: {}", String::from_utf8_lossy(&o.stderr));
    let mut periods = synth.as_os_str().to_owned();
    periods.push(".periods.tsv");
    let det = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(synth.with_extension("fa"))
        .arg("--periods")
        .arg(PathBuf::from(periods))
        .arg("-o")
        .arg(&det)
        .output()
        .unwrap();
    assert!(o.status.success(), "detect failed for {stem}: {}", String::from_utf8_lossy(&o.stderr));
    (dir, det)
}

fn read_properties_row(det: &Path) -> Vec<String> {
    let mut p = det.as_os_str().to_owned();
    p.push(".properties.tsv");
    let s = std::fs::read_to_string(PathBuf::from(p)).unwrap();
    s.lines().nth(1).unwrap().split('\t').map(String::from).collect()
}

fn check_fixture(stem: &str, exps: &HashMap<String, Expectation>) {
    let exp = exps.get(stem).expect("expectation row missing");
    let cfg_yaml = corpus_dir().join(format!("{stem}.yaml"));
    assert!(
        cfg_yaml.exists(),
        "corpus YAML not found: {}",
        cfg_yaml.display()
    );
    let (_dir, det) = synth_and_detect(&cfg_yaml);
    let row = read_properties_row(&det);
    let class = &row[2];
    let base_str = &row[3];
    let hor_k_str = &row[4];
    let reason = &row[19];

    assert_eq!(
        class, &exp.expected_class,
        "{stem}: expected class={} got={}; reason={}",
        exp.expected_class, class, reason
    );

    if let Some(expected_bw) = exp.expected_base_width_bp {
        let got: usize = base_str.parse().expect("base_width_bp not int");
        assert!(
            got.abs_diff(expected_bw) <= exp.base_width_tol_bp,
            "{stem}: expected base_width={} (±{}) got {}; reason={}",
            expected_bw, exp.base_width_tol_bp, got, reason
        );
    }

    if let Some(expected_k) = exp.expected_hor_k {
        let got: usize = hor_k_str.parse().expect("hor_k not int");
        assert_eq!(
            got, expected_k,
            "{stem}: expected hor_k={} got {}; reason={}",
            expected_k, got, reason
        );
    } else {
        assert_eq!(
            hor_k_str, "NA",
            "{stem}: expected hor_k=NA, got {}; reason={}",
            hor_k_str, reason
        );
    }

    if !exp.expected_reason_contains.is_empty() {
        assert!(
            reason.contains(&exp.expected_reason_contains),
            "{stem}: expected reason to contain {:?}; got {:?}",
            exp.expected_reason_contains, reason
        );
    }
}

#[test]
fn core_t01_simple_tr() {
    check_fixture("T01_simple_tr", &load_expectations());
}

#[test]
fn core_t05_hor_clean() {
    check_fixture("T05_hor_clean", &load_expectations());
}

#[test]
fn core_t06_regime_a() {
    check_fixture("T06_regime_A", &load_expectations());
}

#[test]
fn core_t07_regime_c() {
    check_fixture("T07_regime_C", &load_expectations());
}

#[test]
fn core_t10_phase_shift() {
    check_fixture("T10_phase_shift", &load_expectations());
}

#[test]
fn core_t13_coexisting_periods() {
    check_fixture("T13_coexisting_periods", &load_expectations());
}

#[test]
fn core_t17_random_negative() {
    check_fixture("T17_random_negative", &load_expectations());
}

#[test]
fn core_t18_at_rich() {
    check_fixture("T18_at_rich", &load_expectations());
}

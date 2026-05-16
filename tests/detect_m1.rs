//! M1 acceptance: `width_features.tsv` is populated with
//! background-corrected column IC and fraction-conserved-columns at
//! every tested width. A clean simple-TR at L=170 must show much
//! higher IC at width=170 than at the random decoy.

use std::io::Write;
use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn synth_one(yaml: &str) -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let cfg = dir.path().join("cfg.yaml");
    std::fs::File::create(&cfg).unwrap().write_all(yaml.as_bytes()).unwrap();
    let prefix = dir.path().join("arr");
    let o = Command::new(kitehor_bin())
        .args(["synth"])
        .arg(&cfg)
        .arg("-o")
        .arg(&prefix)
        .output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    (dir, prefix)
}

fn parse_width_features(prefix: &std::path::Path) -> Vec<Vec<String>> {
    let mut p = prefix.as_os_str().to_owned();
    p.push(".width_features.tsv");
    let s = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    let mut out = Vec::new();
    for (i, line) in s.lines().enumerate() {
        if i == 0 {
            continue;
        }
        out.push(line.split('\t').map(str::to_string).collect());
    }
    out
}

#[test]
fn simple_tr_has_high_ic_at_true_width() {
    let (dir, synth_prefix) = synth_one(
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
    let fa = synth_prefix.with_extension("fa");
    let periods = {
        let mut p = synth_prefix.as_os_str().to_owned();
        p.push(".periods.tsv");
        std::path::PathBuf::from(p)
    };
    let det_prefix = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(&fa)
        .arg("--periods")
        .arg(&periods)
        .arg("-o")
        .arg(&det_prefix)
        .output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let rows = parse_width_features(&det_prefix);
    let ic_at = |target: usize| -> Option<f64> {
        rows.iter()
            .find(|r| r[1].parse::<usize>().ok() == Some(target))
            .and_then(|r| r[3].parse::<f64>().ok())
    };
    let ic170 = ic_at(170).expect("expected width=170 in width_features");
    assert!(
        ic170 > 1.5,
        "simple_TR base width should have high column IC; got {ic170}"
    );

    // A decoy width — find any width != 170 and != 340 (harmonic) in
    // the output; its IC should be much lower.
    let decoy = rows
        .iter()
        .find(|r| {
            let w: usize = r[1].parse().unwrap_or(0);
            w != 170 && w != 340 && w != 510 && w != 169 && w != 171
        })
        .map(|r| r[3].parse::<f64>().unwrap_or(f64::NAN));
    if let Some(d) = decoy {
        if !d.is_nan() {
            assert!(
                d < ic170 - 0.5,
                "decoy width IC ({d}) should be much lower than true-width IC ({ic170})"
            );
        }
    }
}

#[test]
fn unsupported_widths_emit_na() {
    // Force a very-low array (only enough for ~3 rows at width=100)
    // and check that any width producing fewer than min_rows rows
    // shows up as `class_hint=unsupported`.
    let (dir, synth_prefix) = synth_one(
        r#"
schema_version: 1
seed: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 100
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 30
"#,
    );
    let fa = synth_prefix.with_extension("fa");
    let periods = {
        let mut p = synth_prefix.as_os_str().to_owned();
        p.push(".periods.tsv");
        std::path::PathBuf::from(p)
    };
    let det_prefix = dir.path().join("det");
    let o = Command::new(kitehor_bin())
        .arg("detect")
        .arg(&fa)
        .arg("--periods")
        .arg(&periods)
        .arg("-o")
        .arg(&det_prefix)
        .output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    let rows = parse_width_features(&det_prefix);
    // Look for any width >= 500 (would yield ≤ 6 rows, below default 8).
    let big_widths: Vec<&Vec<String>> = rows
        .iter()
        .filter(|r| r[1].parse::<usize>().unwrap_or(0) >= 500)
        .collect();
    if !big_widths.is_empty() {
        for r in big_widths {
            assert_eq!(r[3], "NA", "wide widths must emit NA for column_IC");
            assert_eq!(r[16], "unsupported");
        }
    }
}

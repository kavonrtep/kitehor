//! M0 acceptance: `kitehor detect` reads FASTA + periods.tsv, validates
//! config, and writes schema-correct (header + per-array placeholder)
//! TSVs. No detection logic is exercised — just IO + scaffolding.

use std::process::Command;

fn kitehor_bin() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn synth_one(config_yaml: &str) -> tempfile::TempDir {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let cfg = dir.path().join("cfg.yaml");
    std::fs::File::create(&cfg).unwrap().write_all(config_yaml.as_bytes()).unwrap();
    let prefix = dir.path().join("arr");
    let out = Command::new(kitehor_bin())
        .args(["synth"])
        .arg(&cfg)
        .arg("-o")
        .arg(&prefix)
        .output()
        .expect("run synth");
    assert!(out.status.success(), "synth failed: {}", String::from_utf8_lossy(&out.stderr));
    dir
}

fn run_detect(fa: &std::path::Path, periods: &std::path::Path, out_prefix: &std::path::Path) -> std::process::Output {
    Command::new(kitehor_bin())
        .arg("detect")
        .arg(fa)
        .arg("--periods")
        .arg(periods)
        .arg("-o")
        .arg(out_prefix)
        .output()
        .expect("run detect")
}

#[test]
fn detect_writes_header_only_segments_for_single_array() {
    // Use the simulator to produce one valid input bundle, then run
    // detect over it. M0 doesn't classify; just verify the bundle.
    let dir = synth_one(
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
    n_copies: 100
"#,
    );
    let fa = dir.path().join("arr.fa");
    let periods = dir.path().join("arr.periods.tsv");
    let out_prefix = dir.path().join("det");

    let out = run_detect(&fa, &periods, &out_prefix);
    assert!(
        out.status.success(),
        "detect failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // properties.tsv: header + one row
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".properties.tsv");
    let props = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    let lines: Vec<&str> = props.lines().collect();
    assert_eq!(lines.len(), 2, "expected header + 1 row, got {}", lines.len());
    let header_cols = lines[0].split('\t').count();
    let row_cols = lines[1].split('\t').count();
    assert_eq!(header_cols, 20, "properties has 20 columns");
    assert_eq!(row_cols, 20);
    let fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(fields[0], "arr");
    // length_bp = 170 * 100 = 17000
    assert_eq!(fields[1], "17000");
    // `class` must be one of the documented values (M4+ replaces M0's
    // placeholder `ambiguous` with a real call).
    assert!(
        matches!(
            fields[2],
            "simple_TR" | "HOR" | "irregular_HOR" | "mixed" | "ambiguous"
        ),
        "unexpected class value {:?}",
        fields[2]
    );

    // segments.tsv: header only
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".segments.tsv");
    let segs = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    assert_eq!(segs.lines().count(), 1);
    assert_eq!(segs.lines().next().unwrap().split('\t').count(), 11);

    // width_features.tsv: header + one row per tested width (M1).
    let mut p = out_prefix.as_os_str().to_owned();
    p.push(".width_features.tsv");
    let wf = std::fs::read_to_string(std::path::PathBuf::from(p)).unwrap();
    let wf_lines: Vec<&str> = wf.lines().collect();
    assert!(wf_lines.len() >= 2, "width_features should have header + ≥1 data row");
    assert_eq!(wf_lines[0].split('\t').count(), 17);
    let cols: Vec<&str> = wf_lines[1].split('\t').collect();
    assert_eq!(cols.len(), 17);
    assert_eq!(cols[0], "arr");
    // M1: column_IC and fraction_conserved_columns are populated (not NA)
    // for widths that satisfy min_rows_per_width.
    let ic = cols[3];
    let fc = cols[4];
    assert!(
        ic != "NA" || fc != "NA",
        "expected M1 to populate column_IC/fraction_conserved for at least the first width"
    );
}

#[test]
fn detect_fails_clearly_on_missing_periods_file() {
    let dir = synth_one(
        r#"
schema_version: 1
templates:
  m:
    type: monomer
    monomer_length_bp: 100
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: 50
"#,
    );
    let fa = dir.path().join("arr.fa");
    let missing = dir.path().join("nope.periods.tsv");
    let out_prefix = dir.path().join("det");
    let out = run_detect(&fa, &missing, &out_prefix);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("periods") || stderr.contains("No such"),
        "expected periods error in stderr, got: {stderr}"
    );
}

#[test]
fn detect_batch_pairs_by_stem() {
    use std::io::Write;
    // Synth two tiny arrays into a dir, then batch detect.
    let dir = tempfile::tempdir().unwrap();
    let cfg_a = dir.path().join("a.yaml");
    let cfg_b = dir.path().join("b.yaml");
    std::fs::File::create(&cfg_a).unwrap().write_all(
        b"schema_version: 1\nseed: 1\ntemplates:\n  m:\n    type: monomer\n    monomer_length_bp: 100\nstructure:\n  - type: SIMPLE_TR\n    template: m\n    n_copies: 20\n",
    ).unwrap();
    std::fs::File::create(&cfg_b).unwrap().write_all(
        b"schema_version: 1\nseed: 2\ntemplates:\n  m:\n    type: monomer\n    monomer_length_bp: 150\nstructure:\n  - type: SIMPLE_TR\n    template: m\n    n_copies: 20\n",
    ).unwrap();
    let synth_out = dir.path().join("synth");
    std::fs::create_dir(&synth_out).unwrap();
    // synth-batch the configs
    let out = Command::new(kitehor_bin())
        .args(["synth-batch", "--config-dir"])
        .arg(dir.path())
        .arg("--out-dir")
        .arg(&synth_out)
        .output().unwrap();
    assert!(out.status.success(), "synth-batch failed: {}", String::from_utf8_lossy(&out.stderr));

    // Now split fasta and periods into separate sibling dirs.
    let fa_dir = dir.path().join("fa");
    let pe_dir = dir.path().join("pe");
    std::fs::create_dir(&fa_dir).unwrap();
    std::fs::create_dir(&pe_dir).unwrap();
    for entry in std::fs::read_dir(&synth_out).unwrap() {
        let p = entry.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if name.ends_with(".fa") {
            std::fs::copy(&p, fa_dir.join(&name)).unwrap();
        } else if name.ends_with(".periods.tsv") {
            std::fs::copy(&p, pe_dir.join(&name)).unwrap();
        }
    }

    let det_out = dir.path().join("det");
    std::fs::create_dir(&det_out).unwrap();
    let out = Command::new(kitehor_bin())
        .arg("detect-batch")
        .arg("--fasta-dir")
        .arg(&fa_dir)
        .arg("--periods-dir")
        .arg(&pe_dir)
        .arg("--out-dir")
        .arg(&det_out)
        .output().unwrap();
    assert!(out.status.success(), "detect-batch failed: {}", String::from_utf8_lossy(&out.stderr));

    // Each stem should produce three TSVs.
    for stem in ["a", "b"] {
        for ext in ["properties.tsv", "segments.tsv", "width_features.tsv"] {
            let p = det_out.join(format!("{stem}.{ext}"));
            assert!(p.exists(), "missing {:?}", p);
        }
    }
}

#[test]
fn detect_expectations_oracle_well_formed() {
    // Sanity: the oracle file has the expected columns for every active
    // CI fixture, and `expected_class` is one of the documented values.
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("detect_expectations.tsv");
    let text = std::fs::read_to_string(&path).expect("read oracle");
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(cols, [
        "array_id", "expected_class", "expected_base_width_bp",
        "expected_hor_k", "base_width_tol_bp", "expected_reason_contains",
        "notes",
    ]);
    let mut n = 0;
    for line in lines {
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 7, "row `{}` has wrong column count", line);
        assert!(matches!(
            fields[1],
            "simple_TR" | "HOR" | "irregular_HOR" | "mixed" | "ambiguous"
        ), "bad class in row `{}`", line);
        n += 1;
    }
    // 22 active fixtures (T09 is .deferred.yaml).
    assert_eq!(n, 22, "expected 22 oracle rows, got {n}");
}

//! End-to-end smoke test for `kitehor rescore`.
//!
//! Pipeline: kite-periodicity on the smoke fixture → rescore the resulting
//! peaks.tsv. Asserts that:
//!
//! 1. The CLI runs to completion and emits `<prefix>.peaks.tsv`.
//! 2. The output schema is the input schema + the four identity columns,
//!    with row count preserved.
//! 3. The headline correctness claim: for HOR fixtures the HOR-unit
//!    period scores higher in `identity_med` than the monomer-sized
//!    period. This is the whole point of the stage.
//! 4. `--force` is required to overwrite an existing output file.

use std::collections::HashMap;
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

fn run(cmd: &mut Command) {
    let o = cmd.output().expect("failed to spawn kitehor");
    assert!(
        o.status.success(),
        "{:?} failed:\nstdout: {}\nstderr: {}",
        cmd,
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}

/// Build {(case_id, period) → identity_med} from a rescored peaks file.
fn parse_identity(path: &Path) -> HashMap<(String, usize), f64> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut lines = content.lines();
    let header = lines.next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let icase = cols.iter().position(|c| *c == "case_id").unwrap();
    let iper = cols.iter().position(|c| *c == "period").unwrap();
    let imed = cols.iter().position(|c| *c == "identity_med").unwrap();
    let mut out = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells[imed] == "NA" {
            continue;
        }
        let case = cells[icase].to_string();
        let per: usize = cells[iper].parse().unwrap();
        let med: f64 = cells[imed].parse().unwrap();
        out.insert((case, per), med);
    }
    out
}

#[test]
fn rescore_pipeline_appends_columns_and_separates_hor_from_monomer() {
    let dir = tempfile::tempdir().unwrap();
    let kite_prefix = dir.path().join("smoke.kite");
    let kite_peaks = {
        let mut p = kite_prefix.as_os_str().to_owned();
        p.push(".tsv.peaks.tsv");
        PathBuf::from(p)
    };
    let kite_out = {
        let mut p = kite_prefix.as_os_str().to_owned();
        p.push(".tsv");
        PathBuf::from(p)
    };
    let rescore_prefix = dir.path().join("smoke.rescore");
    let rescore_out = {
        let mut p = rescore_prefix.as_os_str().to_owned();
        p.push(".peaks.tsv");
        PathBuf::from(p)
    };

    // Step 1: kite-periodicity → peaks.tsv
    run(Command::new(kitehor_bin())
        .arg("kite-periodicity")
        .arg(smoke_fasta())
        .arg("-o")
        .arg(&kite_out));

    assert!(kite_peaks.exists(), "kite did not produce {:?}", kite_peaks);

    // Step 2: rescore. Use --top-n=5 + small --samples so the test
    // doesn't pay the O(P²) DP cost on kite's long-period tail; the
    // monomer (rank 2) and HOR-unit (rank 1) we assert on are well
    // within the top-5.
    run(Command::new(kitehor_bin())
        .arg("rescore")
        .arg(smoke_fasta())
        .arg("--peaks")
        .arg(&kite_peaks)
        .arg("-o")
        .arg(&rescore_prefix)
        .arg("--samples")
        .arg("50")
        .arg("--top-n")
        .arg("5"));

    assert!(rescore_out.exists());

    // Header check: must end with the thirteen appended columns.
    let content = std::fs::read_to_string(&rescore_out).unwrap();
    let header = content.lines().next().unwrap();
    assert!(
        header.ends_with(
            "\tidentity_med\tidentity_iqr\tidentity_p25\tidentity_n\tshift_med\tshift_consistency\tphantom\tsubrepeat\tcoverage_frac\tspatial_contrast\tfounder_period\tkmer_autocorr_founder\tkmer_phase_contrast"
        ),
        "unexpected header: {}",
        header
    );

    // Row count check: rescored row count == kite peaks row count.
    let kite_rows = std::fs::read_to_string(&kite_peaks)
        .unwrap()
        .lines()
        .count();
    let rescore_rows = content.lines().count();
    assert_eq!(rescore_rows, kite_rows, "row count must be preserved");

    // Headline correctness: HOR-unit identity > monomer identity for both
    // HOR fixtures.
    let id = parse_identity(&rescore_out);

    let h3_mono = id[&("hor_k3".to_string(), 100)];
    let h3_hor = id[&("hor_k3".to_string(), 300)];
    assert!(
        h3_hor > h3_mono + 0.05,
        "hor_k3: expected identity at HOR unit (300) > monomer (100); got hor={} mono={}",
        h3_hor,
        h3_mono
    );

    let h5_mono = id[&("hor_k5".to_string(), 150)];
    let h5_hor = id[&("hor_k5".to_string(), 750)];
    assert!(
        h5_hor > h5_mono + 0.05,
        "hor_k5: expected identity at HOR unit (750) > monomer (150); got hor={} mono={}",
        h5_hor,
        h5_mono
    );

    // Perfect tandem repeat: identity must be effectively 1.0.
    let pure = id[&("tandem_pure".to_string(), 60)];
    assert!(
        pure > 0.99,
        "tandem_pure period=60 expected identity ~1.0, got {}",
        pure
    );
}

#[test]
fn rescore_refuses_to_overwrite_without_force() {
    let dir = tempfile::tempdir().unwrap();
    let kite_prefix = dir.path().join("smoke.kite");
    let kite_peaks = {
        let mut p = kite_prefix.as_os_str().to_owned();
        p.push(".tsv.peaks.tsv");
        PathBuf::from(p)
    };
    let kite_out = {
        let mut p = kite_prefix.as_os_str().to_owned();
        p.push(".tsv");
        PathBuf::from(p)
    };

    run(Command::new(kitehor_bin())
        .arg("kite-periodicity")
        .arg(smoke_fasta())
        .arg("-o")
        .arg(&kite_out));

    // First rescore: succeeds, emits <prefix>.peaks.tsv next to existing.
    let prefix = dir.path().join("smoke.kite.tsv");
    run(Command::new(kitehor_bin())
        .arg("rescore")
        .arg(smoke_fasta())
        .arg("--peaks")
        .arg(&kite_peaks)
        .arg("-o")
        .arg(&prefix)
        .arg("--force")
        .arg("--samples")
        .arg("10"));

    // Second rescore with same prefix (same output path) — without --force
    // must fail because the output exists.
    let out = Command::new(kitehor_bin())
        .arg("rescore")
        .arg(smoke_fasta())
        .arg("--peaks")
        .arg(&kite_peaks)
        .arg("-o")
        .arg(&prefix)
        .arg("--samples")
        .arg("10")
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "expected failure without --force; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("already exists")
            || String::from_utf8_lossy(&out.stderr).contains("--force"),
        "expected mention of --force in error, got: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

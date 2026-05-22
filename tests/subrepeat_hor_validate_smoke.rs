//! Smoke test for `subrepeat-scan` and `hor-validate`. Builds a small
//! synthetic FASTA via `kitehor simulate`, runs kite + rule-classify
//! upstream, then exercises both stages and asserts their per-stage
//! TSV bundles are populated.

use std::path::PathBuf;
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

#[test]
fn subrepeat_and_hor_validate_run_on_synth() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    // Step 1: synth one HOR record (k=5, monomer=100 → tile=500).
    let fa = dir.join("hor.fa");
    let ok = Command::new(kitehor_bin())
        .args([
            "simulate",
            "--monomer-size",
            "100",
            "--multiplicity",
            "5",
            "--copies",
            "100",
            "--sub-rate-intra",
            "0.04",
            "--sub-rate-inter",
            "0.02",
            "--case-id",
            "test_hor",
            "--out",
        ])
        .arg(&fa)
        .status()
        .unwrap();
    assert!(ok.success(), "simulate failed");

    let kpre = dir.join("k.tsv");
    Command::new(kitehor_bin())
        .args(["kite-periodicity"])
        .arg(&fa)
        .arg("-o")
        .arg(&kpre)
        .status()
        .unwrap();
    let peaks = {
        let mut p = kpre.clone().into_os_string();
        p.push(".peaks.tsv");
        PathBuf::from(p)
    };

    let vpre = dir.join("v");
    Command::new(kitehor_bin())
        .args(["rule-classify"])
        .arg(&peaks)
        .arg("-o")
        .arg(&vpre)
        .status()
        .unwrap();
    let verdicts = {
        let mut p = vpre.into_os_string();
        p.push(".verdicts.tsv");
        PathBuf::from(p)
    };

    // subrepeat-scan
    let srpre = dir.join("sr");
    let ok = Command::new(kitehor_bin())
        .args(["subrepeat-scan"])
        .arg(&fa)
        .arg("-o")
        .arg(&srpre)
        .arg("--kite-peaks")
        .arg(&peaks)
        .status()
        .unwrap();
    assert!(ok.success(), "subrepeat-scan failed");
    let sr_summary = {
        let mut p = srpre.clone().into_os_string();
        p.push(".subrepeat.tsv");
        PathBuf::from(p)
    };
    assert!(sr_summary.exists());
    let body = std::fs::read_to_string(&sr_summary).unwrap();
    let n = body.lines().count() - 1;
    assert!(n >= 1, "subrepeat.tsv should have at least 1 row");

    // hor-validate
    let hvpre = dir.join("hv");
    let ok = Command::new(kitehor_bin())
        .args(["hor-validate"])
        .arg(&fa)
        .arg("--verdicts")
        .arg(&verdicts)
        .arg("--global-peaks")
        .arg(&peaks)
        .arg("-o")
        .arg(&hvpre)
        .status()
        .unwrap();
    assert!(ok.success(), "hor-validate failed");
    let hv_out = {
        let mut p = hvpre.into_os_string();
        p.push(".hor_within_tile.tsv");
        PathBuf::from(p)
    };
    assert!(hv_out.exists());
    let body = std::fs::read_to_string(&hv_out).unwrap();
    // Expect at least 1 row corresponding to the HOR call. The exact
    // density_hint depends on the seed, but it must be one of the
    // documented enum strings.
    assert!(
        body.lines().count() >= 2,
        "hor_within_tile.tsv had no data row"
    );
    let line = body.lines().nth(1).unwrap();
    let hint = line.split('\t').nth(14).unwrap();
    assert!(
        matches!(
            hint,
            "spatially_confirms_hor"
                | "localized_duplication"
                | "ambiguous"
                | "insufficient_phase_bins"
                | "NA"
        ) || hint.starts_with("k_too_low_for_test"),
        "unexpected density_hint = {hint}"
    );
}

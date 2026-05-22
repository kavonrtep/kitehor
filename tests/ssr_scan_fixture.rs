//! `ssr-scan` smoke test against the 10-record synthetic fixture.
//! Runs the Rust port and asserts the regions TSV contains hits whose
//! 1-based-inclusive start / 0-based-exclusive end convention is
//! preserved exactly, and that the consensus path fires for at least
//! one SSR-rich record.

use std::path::PathBuf;
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

#[test]
fn ssr_scan_runs_on_synthetic_fixture() {
    let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tools/rule_proto/subrepeat/synthetic.fasta");
    if !fasta.exists() {
        eprintln!("skipping: fixture not present {fasta:?}");
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    // Step 1: run kite-periodicity to get peaks.
    let kpre = tmp.path().join("k.tsv");
    let status = Command::new(kitehor_bin())
        .args(["kite-periodicity"])
        .arg(&fasta)
        .arg("-o")
        .arg(&kpre)
        .status()
        .unwrap();
    assert!(status.success(), "kite-periodicity failed");
    let peaks = {
        let mut p = kpre.into_os_string();
        p.push(".peaks.tsv");
        PathBuf::from(p)
    };

    // Step 2: run ssr-scan.
    let out_prefix = tmp.path().join("rs");
    let status = Command::new(kitehor_bin())
        .args(["ssr-scan"])
        .arg(&fasta)
        .arg("-o")
        .arg(&out_prefix)
        .arg("--kite-peaks")
        .arg(&peaks)
        .status()
        .unwrap();
    assert!(status.success(), "ssr-scan failed");
    let ssr_tsv = {
        let mut p = out_prefix.clone().into_os_string();
        p.push(".ssr.tsv");
        PathBuf::from(p)
    };
    let regions_tsv = {
        let mut p = out_prefix.into_os_string();
        p.push(".ssr.regions.tsv");
        PathBuf::from(p)
    };

    let ssr_body = std::fs::read_to_string(&ssr_tsv).unwrap();
    let n_records = ssr_body.lines().count() - 1; // minus header
    assert!(
        n_records >= 10,
        "expected ≥10 records in ssr.tsv, got {n_records}"
    );

    // At least one record should have hit the consensus_single or
    // consensus_multi path (the synthetic fixture includes 2 ssr_aware
    // variants per `make_fixtures.py`).
    let method_consensus = ssr_body
        .lines()
        .skip(1)
        .filter(|l| l.contains("consensus_single") || l.contains("consensus_multi"))
        .count();
    assert!(
        method_consensus >= 1,
        "expected ≥1 consensus_* method, got {method_consensus}"
    );

    // regions TSV: hits must obey 1-based start / 0-based-exclusive end.
    let regions_body = std::fs::read_to_string(&regions_tsv).unwrap();
    for line in regions_body.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        let start: usize = cells[5].parse().unwrap();
        let end: usize = cells[6].parse().unwrap();
        assert!(start >= 1, "start must be 1-based: {start}");
        assert!(end >= start, "end {end} < start {start}");
    }
}

//! End-to-end integration test for the FASTA-like periodogram bundle.
//!
//! Verifies both CLI surfaces:
//!   - `kitehor kite-periodicity --periodogram <PATH>`
//!   - `kitehor analyze --periodogram <PATH>`
//!
//! and asserts they produce byte-identical bundles on the smoke fixture.

use std::path::PathBuf;
use std::process::Command;

fn binary_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("target");
    p.push(if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    });
    p.push("kitehor");
    p
}

fn require_binary() -> PathBuf {
    let bin = binary_path();
    assert!(
        bin.exists(),
        "binary {:?} not found — run `cargo build --release` first",
        bin
    );
    bin
}

/// Parses a record's header line of the form
/// `>case_id|channel length=N kmer=K` and returns
/// `(case_id, channel, length, kmer)`.
fn parse_header(line: &str) -> (String, String, usize, usize) {
    assert!(
        line.starts_with('>'),
        "expected '>'-prefixed header: {line}"
    );
    let body = &line[1..];
    let mut toks = body.split_whitespace();
    let id_chan = toks.next().expect("missing id|channel");
    let (id, chan) = id_chan
        .split_once('|')
        .expect("expected `id|channel` in header");
    let mut length = None;
    let mut kmer = None;
    for t in toks {
        if let Some(v) = t.strip_prefix("length=") {
            length = Some(v.parse::<usize>().unwrap());
        } else if let Some(v) = t.strip_prefix("kmer=") {
            kmer = Some(v.parse::<usize>().unwrap());
        }
    }
    (
        id.to_string(),
        chan.to_string(),
        length.expect("missing length="),
        kmer.expect("missing kmer="),
    )
}

#[test]
#[ignore = "requires `cargo build --release` first; run with --ignored"]
fn periodogram_bundle_layout_and_cross_cli_byte_equality() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fasta = manifest.join("test_data/smoke/sequences.fasta");
    assert!(fasta.exists(), "fixture missing: {:?}", fasta);

    let tmp = std::env::temp_dir().join("kitehor_periodogram_bundle");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    let kite_out = tmp.join("kite.tsv");
    let kite_periodogram = tmp.join("kite.periodogram");
    let analyze_prefix = tmp.join("analyze");
    let analyze_periodogram = tmp.join("analyze.periodogram");

    let bin = require_binary();

    // 1. kite-periodicity surface.
    let status = Command::new(&bin)
        .args([
            "kite-periodicity",
            fasta.to_str().unwrap(),
            "-o",
            kite_out.to_str().unwrap(),
            "--periodogram",
            kite_periodogram.to_str().unwrap(),
        ])
        .status()
        .expect("spawn failed");
    assert!(
        status.success(),
        "kite-periodicity exited non-zero: {status}"
    );
    assert!(
        kite_periodogram.exists(),
        "kite-periodicity did not produce {:?}",
        kite_periodogram
    );

    // 2. analyze surface.
    let status = Command::new(&bin)
        .args([
            "analyze",
            fasta.to_str().unwrap(),
            "-o",
            analyze_prefix.to_str().unwrap(),
            "--periodogram",
            analyze_periodogram.to_str().unwrap(),
        ])
        .status()
        .expect("spawn failed");
    assert!(status.success(), "analyze exited non-zero: {status}");
    assert!(
        analyze_periodogram.exists(),
        "analyze did not produce {:?}",
        analyze_periodogram
    );

    // 3. Layout sanity on the kite-periodicity output. The smoke fixture
    // has 3 records, so we expect 3 × 2 = 6 record headers and 6 value lines.
    let body = std::fs::read_to_string(&kite_periodogram).unwrap();
    let lines: Vec<&str> = body.lines().collect();
    assert_eq!(
        lines.len(),
        6 * 2,
        "expected 12 lines (3 records × 2 channels × 2 lines/record), got {}",
        lines.len()
    );

    // Walk the bundle: each pair (header, values).
    let mut seen_ids: Vec<(String, String)> = Vec::new();
    for chunk in lines.chunks(2) {
        let (id, chan, length, kmer) = parse_header(chunk[0]);
        assert_eq!(kmer, 6, "default kite kmer is 6");
        assert!(
            matches!(chan.as_str(), "H" | "bg"),
            "unknown channel {chan:?}"
        );
        let vals: Vec<&str> = chunk[1].split_whitespace().collect();
        assert_eq!(
            vals.len(),
            length,
            "channel {} of {} has {} values, expected length={}",
            chan,
            id,
            vals.len(),
            length
        );
        if chan == "H" {
            // Integer-valued.
            for v in &vals {
                v.parse::<i64>()
                    .unwrap_or_else(|_| panic!("non-integer H value: {v:?}"));
            }
        } else {
            for v in &vals {
                v.parse::<f64>()
                    .unwrap_or_else(|_| panic!("non-float bg value: {v:?}"));
            }
        }
        seen_ids.push((id, chan));
    }
    // 3 records × (H, bg) — order is record-major, channel-minor.
    assert_eq!(seen_ids.len(), 6);
    for window in seen_ids.chunks(2) {
        assert_eq!(window[0].0, window[1].0, "channels not paired by record");
        assert_eq!(window[0].1, "H");
        assert_eq!(window[1].1, "bg");
    }

    // 4. Cross-CLI byte equality. The auto path through `analyze` calls
    // the same `kite::analyze` with the same defaults, so the bundles
    // must match exactly.
    let kite_bytes = std::fs::read(&kite_periodogram).unwrap();
    let analyze_bytes = std::fs::read(&analyze_periodogram).unwrap();
    assert_eq!(
        kite_bytes, analyze_bytes,
        "kite-periodicity and analyze produced different periodogram bundles"
    );
}

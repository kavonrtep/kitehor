//! Python prototype parity test for `kitehor tandem-validate`.
//!
//! Generates a small synthetic fixture in a temp dir, runs the upstream
//! kite + rule-classify stages to produce verdicts + peaks, then runs
//! both the Rust and Python implementations of `tandem_validate`
//! against the same inputs and asserts the `decision_hint` column
//! matches for every record.
//!
//! Marked `#[ignore]` because it needs `python3` + `pandas` + the
//! reference Python prototype at `tools/rule_proto/tandem_validate.py`
//! on the local filesystem. CI doesn't run it; engineers run it
//! manually before tagging a release that touches the detector:
//!
//! ```text
//! cargo test --release --test tandem_validate_python_parity -- --ignored
//! ```
//!
//! The fixture is generated deterministically (fixed seeds) so the
//! parity comparison is reproducible across runs.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn kitehor_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kitehor"))
}

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// One inline FASTA record: a `simple_tr` whose monomer carries an
/// embedded subrepeat block. The construction matches the canonical
/// "nested microsatellite inside a TR monomer" case the unified
/// detector is designed to flag (`localized_subrepeat`).
///
/// `monomer_len` = total monomer bp. `ms_motif` (≥ 20 bp to clear the
/// detector's default `cand_min_period = 20`) is repeated `ms_copies`
/// times to form a contiguous internal block at `ms_offset`. The
/// monomer is then repeated `n_copies` times.
///
/// The "random" backbone is a deterministic LCG-style PRNG over ACGT
/// so the test is reproducible without pulling in `rand`.
fn nested_subrepeat_record(
    rec_id: &str,
    monomer_len: usize,
    ms_motif: &[u8],
    ms_copies: usize,
    ms_offset: usize,
    n_copies: usize,
    seed: u64,
) -> (String, Vec<u8>) {
    let ms_block: Vec<u8> = ms_motif
        .iter()
        .cycle()
        .take(ms_motif.len() * ms_copies)
        .copied()
        .collect();
    assert!(
        ms_offset + ms_block.len() <= monomer_len,
        "MS block overflows monomer"
    );
    // Build one monomer: deterministic random backbone + embedded MS block.
    let bases = b"ACGT";
    let mut state = seed.wrapping_mul(0x9E3779B97F4A7C15);
    let mut monomer: Vec<u8> = (0..monomer_len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) as usize) & 3]
        })
        .collect();
    for (i, &c) in ms_block.iter().enumerate() {
        monomer[ms_offset + i] = c;
    }
    // Repeat the monomer with no inter-copy mutation — for parity-test
    // purposes we just need the kite + detector to find a structure.
    let mut seq: Vec<u8> = Vec::with_capacity(monomer_len * n_copies);
    for _ in 0..n_copies {
        seq.extend_from_slice(&monomer);
    }
    (rec_id.to_string(), seq)
}

/// Write a multi-record FASTA to `path`.
fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    use std::io::Write;
    let mut f = std::fs::File::create(path).expect("create fixture FASTA");
    for (id, seq) in records {
        writeln!(f, ">{id}").unwrap();
        // 80-bp lines (anyone who reads the file by eye will thank us).
        for chunk in seq.chunks(80) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
    }
}

/// Run a single `kitehor simulate` call, returning the simulator's FASTA
/// as `(case_id, seq)`. Truth metadata is discarded.
fn simulate(
    case_id: &str,
    monomer_size: usize,
    multiplicity: usize,
    copies: usize,
    seed: u64,
    tmp: &Path,
) -> (String, Vec<u8>) {
    let out = tmp.join(format!("{case_id}.fa"));
    let status = Command::new(kitehor_bin())
        .args(["simulate"])
        .args(["--monomer-size", &monomer_size.to_string()])
        .args(["--multiplicity", &multiplicity.to_string()])
        .args(["--copies", &copies.to_string()])
        .args(["--seed", &seed.to_string()])
        .args(["--case-id", case_id])
        .arg("--out")
        .arg(&out)
        .status()
        .expect("spawn simulate");
    assert!(status.success(), "simulate failed for {case_id}");
    let body = std::fs::read_to_string(&out).expect("read simulate FASTA");
    let mut id = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in body.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            id = rest.split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.extend_from_slice(line.as_bytes());
        }
    }
    assert!(!seq.is_empty(), "empty simulate output for {case_id}");
    (id, seq)
}

/// Parse a TSV into `record_id → decision_hint` map. Assumes the input
/// has a `record_id` column and a `decision_hint` column.
fn parse_decisions(path: &Path) -> HashMap<String, String> {
    let body = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {path:?}: {e}"));
    let mut lines = body.lines();
    let header: Vec<&str> = lines.next().expect("empty TSV").split('\t').collect();
    let id_idx = header
        .iter()
        .position(|c| *c == "record_id")
        .expect("missing record_id column");
    let dec_idx = header
        .iter()
        .position(|c| *c == "decision_hint")
        .expect("missing decision_hint column");
    let mut out = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        out.insert(cells[id_idx].to_string(), cells[dec_idx].to_string());
    }
    out
}

#[test]
#[ignore = "needs python3 + pandas + tools/rule_proto/tandem_validate.py"]
fn decision_hint_matches_python_prototype() {
    let prototype = manifest_dir().join("tools/rule_proto/tandem_validate.py");
    if !prototype.exists() {
        panic!(
            "Python prototype missing at {prototype:?}; this test needs the \
             reference implementation in tools/rule_proto/"
        );
    }
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // ----- 1. Build the fixture FASTA -----
    // Each entry exercises one of the detector's decision paths; per-entry
    // commentary at the call site explains intent.
    let records: Vec<(String, Vec<u8>)> = vec![
        // HOR k=2: tests the skip_k2 path.
        simulate("hor_k2", 200, 2, 60, 101, dir),
        // HOR k=3: founder candidate at the host/3 borderline.
        simulate("hor_k3", 150, 3, 70, 102, dir),
        // HOR k=5: founder well inside host/3, exercises the
        // founder-vs-tile presence test on a uniform tile.
        simulate("hor_k5", 130, 5, 60, 103, dir),
        // HOR k=7: large k where 3*founder may push the window toward host.
        simulate("hor_k7", 100, 7, 50, 104, dir),
        // Pure tandem: simple_tr verdict, host is the monomer itself, no
        // peaks below host/3 should qualify after the relative-floor filter.
        simulate("pure_tr", 200, 1, 100, 105, dir),
        // Nested subrepeat: simple_tr with an embedded ~29-bp periodicity
        // inside a 1500-bp monomer. Exercises the localized_subrepeat path
        // for kind=Other candidates.
        nested_subrepeat_record(
            "nested_sub",
            1500,
            b"AAGCTGACGTAGGCTACAAGCTAAGCTGA", // 29-bp motif, > cand_min_period=20
            12,                               // 12 copies → ~348 bp block
            500,                              // offset within monomer
            40,                               // 40 monomer copies → ~60 kb
            201,
        ),
    ];

    let fixture = dir.join("fixture.fa");
    write_fasta(&fixture, &records);

    // ----- 2. kite + rule-classify (shared inputs for both runs) -----
    let kite_tsv = dir.join("kite.tsv");
    let peaks_tsv = dir.join("kite.peaks.tsv");
    let status = Command::new(kitehor_bin())
        .args(["kite-periodicity"])
        .arg(&fixture)
        .arg("-o")
        .arg(&kite_tsv)
        .arg("--out-peaks")
        .arg(&peaks_tsv)
        .arg("--classify")
        .status()
        .expect("spawn kite-periodicity");
    assert!(status.success(), "kite-periodicity --classify failed");

    // ----- 3. Rust tandem-validate -----
    let rust_prefix = dir.join("rust");
    let status = Command::new(kitehor_bin())
        .args(["tandem-validate"])
        .arg(&fixture)
        .arg("--verdicts")
        .arg(&kite_tsv) // kite-periodicity --classify writes verdicts into the primary TSV
        .arg("--peaks")
        .arg(&peaks_tsv)
        .arg("-o")
        .arg(&rust_prefix)
        .status()
        .expect("spawn tandem-validate");
    assert!(status.success(), "Rust tandem-validate failed");
    let rust_out = {
        let mut p = rust_prefix.clone().into_os_string();
        p.push(".tandem_validate.tsv");
        PathBuf::from(p)
    };
    let rust_decisions = parse_decisions(&rust_out);

    // ----- 4. Python tandem_validate.py -----
    let py_out = dir.join("py.tandem_validate.tsv");
    let status = Command::new("python3")
        .arg(&prototype)
        .arg("--fasta")
        .arg(&fixture)
        .arg("--verdicts")
        .arg(&kite_tsv)
        .arg("--peaks")
        .arg(&peaks_tsv)
        .arg("--out")
        .arg(&py_out)
        .arg("--kite-binary")
        .arg(kitehor_bin())
        .status()
        .expect("spawn python3 tandem_validate.py");
    assert!(status.success(), "Python tandem_validate.py failed");
    let py_decisions = parse_decisions(&py_out);

    // ----- 5. Compare -----
    let mut diffs: Vec<String> = Vec::new();
    let mut all_ids: Vec<&String> = rust_decisions.keys().chain(py_decisions.keys()).collect();
    all_ids.sort();
    all_ids.dedup();
    for id in all_ids {
        let r = rust_decisions.get(id);
        let p = py_decisions.get(id);
        match (r, p) {
            (Some(r), Some(p)) if r == p => {}
            (Some(r), Some(p)) => diffs.push(format!("  {id}: rust={r}  python={p}")),
            (Some(_), None) => diffs.push(format!("  {id}: missing from python output")),
            (None, Some(_)) => diffs.push(format!("  {id}: missing from rust output")),
            (None, None) => {}
        }
    }
    assert!(
        diffs.is_empty(),
        "Rust vs Python decision_hint divergence on {} record(s):\n{}",
        diffs.len(),
        diffs.join("\n")
    );
    assert_eq!(
        rust_decisions.len(),
        records.len(),
        "rust emitted {} rows; expected {}",
        rust_decisions.len(),
        records.len()
    );
}

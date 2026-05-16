# Kite Emit-Periods Integration Review

Date: 2026-05-16  
Reviewed commit: `f9793b6 feat(kite): --emit-periods writes v2 detector periods.tsv`  
Context: `HEAD` is a docs follow-up (`f12e7dd`); the integration code in
`src/emit_periods.rs`, `src/main.rs`, `src/cli.rs`, and `src/lib.rs` is
unchanged since `f9793b6`.

## Summary

The commit adds a useful and mostly clean bridge from `kitehor
kite-periodicity` to the v2 detector by emitting a detector-compatible
`periods.tsv`. The implementation is well isolated in `src/emit_periods.rs`,
has mapper-level unit tests, and documents the score mapping clearly.

The main risks are semantic rather than mechanical: emitted "secondary" periods
can exceed the documented top-3 scope, ML-classifier output is not reflected in
the emitted periods file, and there is no true end-to-end test proving
`kite-periodicity --emit-periods` output is accepted by `detect`.

## Findings

### High: Secondary periods are not limited to Kite's top-3 peaks

The plan and CLI text say "other top-3 peaks", but `append_secondaries()` walks
the full `kr.peaks` list and adds up to three periods after excluding
founder/tile (`src/emit_periods.rs:120`). This means a HOR call can emit
founder + tile + three additional peaks, and those additional peaks may come
from rank 4, 5, or deeper.

This matters because the detector still evaluates low-score rows. A score of
`0.60` prevents the strong-period shortcut, but it does not prevent a candidate
from becoming a final class through the canonical column-IC path. Extra harmonic
or noisy Kite peaks can therefore increase false `mixed` or wrong-width calls.

Recommendation: decide whether the contract is "top-3 total Kite periods" or
"up to three extra secondaries". If it is top-3, filter secondaries to
`kr.peaks.iter().take(3)`. If extra secondaries are intentional, update the
plan/help text and add a detector-level regression showing they do not harm
classification.

### Medium: `--use-ml-classifier` predictions are ignored by `--emit-periods`

In `run_kite_periodicity()`, `--emit-periods` uses rule verdicts only when
`--classify && !--use-ml-classifier`; otherwise it falls back to raw Kite peaks
at `0.60` (`src/main.rs:423`). The ML classifier can produce `founder`, `tile`,
and recovered HOR calls (`src/classify.rs:116`), but those high-confidence
periods are not emitted.

The result can be confusing: `predictions.tsv` may report an ML `hor` with
founder/tile, while the emitted periods file contains only low-confidence raw
hints. Feeding that file into `detect` loses the classifier's strongest signal.

Recommendation: either explicitly reject `--emit-periods` with
`--use-ml-classifier`, or map ML `Verdict` values into the same periods schema
when `founder`/`tile` are available. If the fallback is intentional, document it
as a limitation in the CLI help and README.

### Medium: Same-FASTA pipeline can fail when Kite QC skips records

`kite-periodicity` filters loaded FASTA records to `LoadStatus::Ok` before
running Kite (`src/main.rs:244`). `--emit-periods` writes rows only for those
`results`. The v2 detector, however, loads FASTA permissively and attempts to
join every record to `periods.tsv`. A record skipped by Kite QC but still
present in the same FASTA will have no period rows, causing `detect` to fail
unless `--allow-missing-periods` is supplied.

This is similar to the documented `NoSignal` case, but not the same. The README
mentions `NoSignal`, not QC-skipped records.

Recommendation: document that `--allow-missing-periods` is also required when
the original FASTA contains records rejected by Kite QC, or emit a small
sidecar/count of skipped IDs so the user can see why `detect` needs that flag.

### Medium: Missing end-to-end coverage for the new bridge

`src/emit_periods.rs` has useful unit tests for row mapping and TSV formatting,
but there is no integration test that runs:

```bash
kitehor kite-periodicity input.fa --classify --emit-periods kite.periods.tsv
kitehor detect input.fa --periods kite.periods.tsv -o det
```

The wiring in `src/main.rs` is therefore not protected against regressions in
CLI parsing, output path handling, multi-record joins, or detector strictness.

Recommendation: add one small integration test using an existing synthetic
fixture. Assert that `--emit-periods` writes the expected header and at least
one high-score row, then run `detect` on the same FASTA and verify a
schema-valid properties file.

### Low: Source labels are slightly misleading for tandem calls

For `RuleVerdict::Tandem`, the emitted high-score monomer uses source
`kite_founder` (`src/emit_periods.rs:87`). The detector currently ignores
`source`, so this is not a behavioral bug. It may still confuse downstream
users because a simple tandem monomer is not a HOR founder.

Recommendation: consider `kite_monomer` or `kite_primary` for tandem calls, or
document that `kite_founder` means "primary trusted base period" in this file.

## Verification Notes

I could not run Rust tests locally:

```text
cargo test --release emit_periods --lib
error: failed to parse lock file ... Cargo.lock
Caused by: lock file version 4 requires `-Znext-lockfile-bump`
```

The installed local toolchain is Cargo/Rust `1.75.0`, while this repository
declares `rust-version = "1.85"` and pins toolchain `1.95`.

I also attempted a small CLI smoke with `target/release/kitehor`, but that
binary predates the reviewed commit:

```text
target/release/kitehor timestamp: 2026-05-16 15:22
reviewed commit time:            2026-05-16 20:24
```

As expected, the old binary rejected `--emit-periods`.

## Suggested Priority

1. Fix or clarify the secondary-period cap semantics.
2. Add an end-to-end integration test for `--emit-periods` followed by
   `detect`.
3. Decide how `--emit-periods` should behave with `--use-ml-classifier`.
4. Document the QC-skipped-record case alongside `NoSignal`.
5. Optionally rename the tandem source label for clarity.

## Bottom Line

The integration is directionally sound and small enough to maintain. Before it
is treated as a stable Kite-to-detector pipeline, I would tighten the secondary
period contract and add an actual CLI-level regression test. Those two changes
would cover the highest-risk behavior without complicating the design.

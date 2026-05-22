# Rule-based pipeline (`kitehor analyze`)

End-to-end Rust port of the `tools/rule_proto/*.py` prototype. Six
subcommands cover the five stages plus the orchestrator:

```
FASTA
  │
  ▼
kite-periodicity       k-mer pair-distance periodogram (existing)
  │
  ▼
rule-classify          HOR / simple_tr / unresolved verdict
  │
  ├──────────────┬──────────────┐
  ▼              ▼              ▼
subrepeat-scan   ssr-scan       hor-validate
(nested-TR)      (short motifs) (within-tile + density)
  │              │              │
  └──────────────┴──────────────┘
                 │
                 ▼
            summary-merge       outer-join + 8-rule combined_class
                 │
                 ▼
          <prefix>.summary.tsv
```

The eight `combined_class` values:

| Class | Fires when |
|---|---|
| `pure_ssr` | `ssr_flag = yes` AND `ssr_dominant_motif_coverage_pct ≥ 80` |
| `tr_with_nested_tr` | `subrepeat_flag = yes` |
| `tr_with_subrepeat` | `hor_verdict = hor` AND `density_hint = localized_duplication` |
| `hor_with_ssr` | `hor_verdict = hor` AND `ssr_flag = yes` |
| `hor` | `hor_verdict = hor` |
| `tr_with_ssr` | `hor_verdict = simple_tr` AND `ssr_flag = yes` |
| `tr` | `hor_verdict = simple_tr` |
| `unresolved` | none of the above |

`tr_with_nested_tr` and `tr_with_subrepeat` are the same biological
phenomenon at different scales: a TR whose monomer contains internal
repetition. The first is detected by the subrepeat-scan (≥ few sub-TR
copies per monomer); the second by the hor-validate density check
(localized 2–3-copy duplication within a HOR-mis-called monomer).

## Quick start

```bash
# End-to-end
kitehor analyze <fasta> -o <prefix>

# Per stage (debugging / partial rerun)
kitehor kite-periodicity <fasta> -o <prefix>.kite.tsv
kitehor rule-classify   <prefix>.kite.tsv.peaks.tsv -o <prefix>
kitehor subrepeat-scan  <fasta> --kite-peaks <prefix>.kite.tsv.peaks.tsv -o <prefix>
kitehor ssr-scan        <fasta> --kite-peaks <prefix>.kite.tsv.peaks.tsv -o <prefix>
kitehor hor-validate    <fasta> --verdicts <prefix>.verdicts.tsv \
                                --global-peaks <prefix>.kite.tsv.peaks.tsv -o <prefix>
kitehor summary-merge   --verdicts <prefix>.verdicts.tsv \
                        --subrepeat <prefix>.subrepeat.tsv \
                        --ssr <prefix>.ssr.tsv \
                        --within-tile <prefix>.hor_within_tile.tsv \
                        -o <prefix>
```

## Outputs (TSV-per-stage contract)

`analyze` always emits all eight per-stage TSVs under `<prefix>.*`:

| Stage | File(s) | Column count |
|---|---|---:|
| kite-periodicity | `.kite.tsv`, `.kite.peaks.tsv` | 9 / 9 |
| rule-classify | `.verdicts.tsv` | 10 |
| subrepeat-scan | `.subrepeat.tsv`, `.windows.tsv` | 13 / 7 |
| ssr-scan | `.ssr.tsv`, `.ssr.regions.tsv` | 17 / 8 |
| hor-validate | `.hor_within_tile.tsv` | 16 |
| summary-merge | `.summary.tsv` | 32 |

## Algorithm details

### `rule-classify`

Clusters kite peaks by relative period gap (single-linkage on
`(p_cur − p_prev) / p_cur ≤ tol`), then walks a first-match-wins
decision tree:

1. **Case A** — top cluster is `k × shorter cluster`. Picks the
   smallest-period qualifying divisor. HOR(founder, k, tile=top).
2. **Case B walk** — top is the founder; k = 2..k_max looks for a
   cluster at `k × top.period` that produces a non-monotonic bump
   above the running max (k ≥ 3) or exceeds `non_mono_ratio × top`
   (k = 2). Requires the harmonic confirmation
   `score(2k·p) ≥ score((k+1)·p)`.
3. **monotonic_multiples** — has at least one larger-multiple cluster
   but no qualifying bump → `simple_tr`.
4. **lone_significant_cluster** — exactly one cluster passes
   `lone_significant_frac × top` → `simple_tr`.
5. **unresolved** — none of the above.

All hardcoded constants exposed as CLI flags. See `--help`.

### `subrepeat-scan`

For each record:

1. From kite peaks pick `(sub_candidate, host_candidate)` where the
   sub is the shortest qualifying period in the top-N by score, and
   the host is the strongest-scored period at least
   `host_sub_ratio_min ×` the sub.
2. Slide windows of size `max(window_mult_sub × sub_candidate,
   min_window_bp)` across the array with step `window // step_frac`.
3. Run kite **in-process** on each window using a window-id-derived
   FNV-1a seed (byte-equivalent with the prototype's subprocess
   invocations).
4. Classify each window as `sub` if its rank-1 peak is within `tol`
   of `sub_candidate` AND score ≥ `window_score_floor`.
5. Morphological smoothing: absorb runs shorter than `min_run` into
   neighbours (tie-break: previous neighbour wins when
   `prev_len ≥ next_len`).
6. Build contiguous `sub`-blocks. Flag `yes` iff at least one block
   AND at least one `non_sub` window exist.

### `ssr-scan`

1. Raw scan — for each motif length L = 1..14, find non-overlapping
   greedy runs of an L-base motif repeating ≥ `min_reps` times. The
   prototype uses regex `(([gatc]{L})\2{min-1,})`; the Rust port
   hand-rolls the scanner (Rust's `regex` crate doesn't support
   backreferences). Output is `(motif_lower, start_1based,
   end_0based_exclusive, repeats, normalized_motif)`.
2. `normalize_motif` = lex-min over all rotations of `motif.upper()`
   AND of its reverse complement.
3. `get_unique_motifs` drops any motif that is exactly `k × shorter`.
4. **Consensus path** (when `--kite-peaks` is supplied): from the kite
   top peak's period `P`, extract up to `consensus_max_monomers`
   canonical-distinct P-mers from the sequence; validate each by
   building a dimer (`P × consensus_dimer_copies` or
   `consensus_dimer_min_bp`) and running the raw scan on it. The
   number of unique validated canonicals selects:
   - 0 → `raw_fallback`
   - 1 → `consensus_single` (use the dimer's summary as authoritative)
   - ≥ 2 → `consensus_multi` (per-motif coverage from the raw scan)

### `hor-validate`

For each `hor_verdict = hor` row in the verdicts TSV:

1. **Skip** if `k < min_k_for_density` (default 4) — k = 2 or 3 are
   geometrically degenerate and emit `density_hint =
   k_too_low_for_test(k=N)`.
2. **Within-tile** — slice `seq[0..tile]`, run kite, compute
   `within_founder_top_ratio`. Drive `decision_hint` against
   thresholds `0.5 / 0.2 / 0.05`.
3. **Density windows** — slide windows of size
   `max(tile / density_window_tile_frac, min_founder_mult × founder,
   min_density_window_bp)` (step capped to keep
   ≤ `max_density_windows`). For each, run kite and check whether
   `score_near(founder) ≥ density_rel_floor × top_score`.
4. **Phase fold** into `phase_fold_bins` buckets by
   `(window_midpoint mod tile) // bin_width`; compute `contrast =
   max(frac_per_bin) − min(frac_per_bin)`.
5. **Combined decision**:
   - `localized_duplication` iff `density ≤ density_dup_max OR
     contrast ≥ phase_contrast_dup_min`
   - else `spatially_confirms_hor` iff `density ≥ density_hor_min AND
     contrast ≤ phase_contrast_hor_max`
   - else `ambiguous`

### `summary-merge`

Outer-join on `record_id` (sorted lexicographically — matches
pandas's outer-merge behavior in the user's env). Defaults fill
missing values; the 8-rule decision tree determines
`combined_class`. Float columns use `%.4g`.

## Float formatting policy

Three policies coexist, matching the prototype's per-script
defaults so the Rust output is byte-equivalent with the Python
prototype:

- `rule-classify` → `%.6g`
- `summary-merge` → `%.4g`
- `subrepeat-scan`, `ssr-scan` → pandas default (≈ shortest
  roundtrip + `.0` for integer-valued floats)
- `hor-validate` → `%.6g`

The shared helper `rule_classify::io::fmt_g(precision, x)`
implements Python's `%g` semantics (significand-precision, scientific
fallback for very small/large magnitudes, trailing-zero stripping).

## Byte-equivalence with the Python prototype

All five stages have been validated byte-identical with the
prototype on the smoke fixture (`test_data/smoke/sequences.fasta`)
and on the synthetic regression fixtures
(`tools/rule_proto/subrepeat/synthetic.fasta`,
`tools/rule_proto/fixtures/*.peaks.tsv`).

The `rule-classify` stage is additionally validated against the six
hand-curated fixtures in `tools/rule_proto/fixtures/` via
`tests/rule_classify_fixtures.rs`.

## Replaced modules

`P1` retired `src/rule.rs` (the older 4-condition rule). `P6` removed
the legacy ML pipeline (`classifier.rs`, `classify.rs`, `features.rs`,
`hor_call.rs`, `coverage.rs`, `models/`, `examples/validate_rf.rs`,
`config/classifier.toml`) along with their CLI flags
(`--use-ml-classifier`, `--no-hor-call`, `--hor-*`, `--rule-qmax`,
`--rule-top-n`, `--coverage`, `--classifier-config`, `--hor-model`,
`--k-model`, `--no-homology`).

## Performance

Targets on test_590 (2779 records, ~280 MB FASTA), single machine,
default rayon parallelism:

| Stage | Python wall | Rust target |
|---|---:|---:|
| kite-periodicity | 0:15 | 0:15 |
| rule-classify | 0:10 | < 0:01 |
| subrepeat-scan | 2:50 | < 0:30 |
| ssr-scan | 7:00 | < 0:20 |
| hor-validate | 0:25 | < 0:05 |
| summary-merge | 0:01 | < 0:01 |
| **analyze (end-to-end)** | **~10:40** | **< 1:30** |

# Rule-proto pipeline — Rust port implementation plan

> **Historical document — partially superseded by v0.10.**
>
> The pipeline-shape decisions in this plan still describe the kitehor
> rule-proto pipeline accurately, with one exception: in v0.10 the
> two `subrepeat-scan` + `hor-validate` stages were merged into a
> single unified `tandem-validate` stage, retiring the
> `tr_with_nested_tr` combined class. The current 5-stage / 7-class
> layout is documented in [`../rule_proto.md`](../rule_proto.md);
> the v0.10 retirement rationale + rollout plan live in
> [`tandem_validate_spec.md`](tandem_validate_spec.md) and
> [`tandem_validate_port_plan.md`](tandem_validate_port_plan.md). Read
> those for the current architecture; read this doc for the original
> port intent and the user-confirmed `[Q-N]` design decisions.

**Status**: v2 — 2026-05-22. User answers folded in (§11). Implementation in progress.
**Author**: drafted from `docs/rust_port_onboarding.md`, the prototype
sources at `tools/rule_proto/*.py`, and a fresh audit of the kitehor
Rust crate.
**Scope**: port the five-script Python prototype at
`tools/rule_proto/` into the kitehor Rust crate as native subcommands
plus a single end-to-end orchestrator. Replace the existing
`src/rule.rs` classifier. Preserve the prototype's TSV interfaces
for debugging and downstream analysis. Promote every hardcoded
constant to a CLI flag (defaults from the prototype).

This document is the *plan*. Implementation happens after review.

---

## §0 Decisions snapshot (user-confirmed before drafting)

| # | Decision | Confirmed |
|---|---|---|
| 1 | Replace the existing `src/rule.rs` (do not keep alongside) | yes |
| 2 | Five Python scripts → five kitehor subcommands | yes |
| 3 | Add an end-to-end orchestrator (`kitehor analyze`) | yes |
| 4 | TSV-per-stage exchange format (one file per script, debugging + downstream) | yes |
| 5 | Every hardcoded prototype constant becomes a CLI flag with the prototype default | yes |

Items left open (numbered in §11) are flagged inline as **`[Q‑N]`**
the first time they arise and resolved later.

---

## §1 What we are porting

```
        ┌────────────────────────────┐
        │   kitehor kite-periodicity │   (already Rust)
        │   --classify --emit-…      │
        └─────────────┬──────────────┘
              kite.peaks.tsv
                      │
        ┌─────────────┴──────────────┐
        ▼                            │
  rule-classify                      │
  (verdicts.tsv)                     │
        │                            │
        ├──→ hor-validate ←──────────┤   (peaks + FASTA + verdicts)
        │      (hor_within_tile.tsv) │
        │                            │
        ├──── subrepeat-scan ←───────┤   (FASTA + kite peaks)
        │       (subrepeat.tsv,      │
        │        windows.tsv)        │
        │                            │
        └──── ssr-scan ←─────────────┘   (FASTA + kite peaks)
                (ssr.tsv,
                 ssr.regions.tsv)
                      │
                      ▼
                summary-merge
                (summary.tsv)
                      │
                      ▼
               combined_class
```

The orchestrator `kitehor analyze <fasta>` runs all six stages in
sequence on a single FASTA, emitting all six TSVs under one
`<prefix>.*` namespace.

### 1.1 What we are not porting (out of scope)

- `tools/rule_proto/eval/*.py` and `make_fixtures.py` — eval harness
  and fixture generator stay in Python.
- `tools/rule_proto/hor_density_profile.py` — a diagnostic helper
  not part of the pipeline. Out of scope.
- `tools/training/`, `tools/features/` — legacy ML training pipeline.
  Untouched by this port (but see §6.3 on whether the legacy ML
  code paths in `src/` come out at the same time).

---

## §2 Crate-level architecture

### 2.1 New module layout

Add **five** new top-level modules + **one** shared k-mer utility
+ **one** orchestrator:

```
src/
├── kite.rs                       (existing — pub items expanded; see §2.2)
├── kmer_pairs.rs        NEW     shared pair-distance histogram + window slicer
├── rule_classify/       NEW     §5.1
│   ├── mod.rs                    public API + RuleClassifyConfig
│   ├── cluster.rs                period clustering
│   ├── decide.rs                 Case A / Case B / harmonic check / fallbacks
│   └── io.rs                     peaks.tsv reader + verdicts.tsv writer
├── subrepeat/           NEW     §5.2
│   ├── mod.rs
│   ├── window.rs                 per-window classification + smoothing + blocking
│   └── io.rs
├── ssr/                 NEW     §5.3
│   ├── mod.rs
│   ├── find_ssrs.rs              TideCluster port + normalize_motif + dedup
│   ├── consensus.rs              kite-driven dimer consensus
│   └── io.rs
├── hor_validate/        NEW     §5.4
│   ├── mod.rs
│   ├── within_tile.rs            first-tile kite check
│   ├── density.rs                spatial density + phase folding
│   └── io.rs
├── summary/             NEW     §5.5
│   ├── mod.rs                    outer join + combined_class decision tree
│   └── io.rs
├── analyze.rs           NEW     §5.6 — orchestrator
├── cli.rs                        +6 subcommand variants
└── main.rs                       +6 dispatch arms
```

Rationale for separate top-level dirs (one per stage) rather than
weaving them into `src/detect/`:

- `detect/` is the v2 line-width detector — orthogonal feature.
  Mixing in the rule-proto pipeline would entangle two unrelated
  milestone tracks and make `detect/mod.rs::auto_periods` (which
  consumes the rule classifier) a load-bearing entry point for the
  new pipeline by accident.
- Each new dir owns its TSV reader, writer, config, and tests.
  Internal modules (`cluster.rs`, `window.rs`, …) stay private.

### 2.2 Shared utility: `kmer_pairs.rs` (new)

The single biggest win in this port is fusing the kite-binary
subprocess loops (459 k + 59 k invocations on test_590) into
in-process pair-distance computations.

The shared module exposes:

```rust
/// Per-record k-mer pair-distance histogram, indexed by base
/// position so it can be sliced to a window.
pub struct PairIndex {
    k: usize,
    seq_len: usize,
    /// For each k-mer, its list of start positions in the sequence.
    positions: AHashMap<u64, Vec<u32>>,
}

impl PairIndex {
    pub fn build(seq: &[u8], k: usize) -> Self;
    /// Global H[d] over all kmer-pairs in [0, seq_len).
    pub fn profile_global(&self, max_distance: usize) -> Vec<f64>;
    /// H[d] restricted to k-mer pairs where the *left* member starts
    /// inside [start, end).  Window-scoped pair-distance histogram.
    pub fn profile_window(&self, start: usize, end: usize,
                          max_distance: usize) -> Vec<f64>;
}
```

`subrepeat::window` calls `profile_window(s, e, max_d)` per window
to get the per-window H[d] without re-indexing. `hor_validate::density`
does the same for its density-window sweep. **`[Q‑1]` — should
windowed peak finding reuse kite's background-replicate + smoothing
+ peak-scoring pipeline byte-for-byte, or only the score2-style
extraction? Default decision below: yes, reuse kite's
`find_peaks_with_score` on the windowed profile, with a per-window
re-estimated background derived from the same FNV-1a seed regime as
the global call (see §5.2 and §5.4).**

The existing `src/kite.rs::compute_neighbor_profile` (currently
private, lines 112–149 per the audit) is the reference implementation
of the global histogram. Two clean options:

- **Refactor**: extract `compute_neighbor_profile`,
  `compute_background`, `gaussian_smooth`, and `find_peaks_with_score`
  into `kmer_pairs.rs` + a `kite::peaks` sub-module so they can be
  called on arbitrary `Vec<f64>` profiles (preferred — see §10).
- **Promote-and-reuse**: make those four `pub`, keep them inside
  `kite.rs`. Simpler diff but doesn't isolate the new pipeline from
  the existing kite-periodicity output formatter.

This plan assumes the **refactor** option. It is the only structural
change required in existing kite code.

### 2.3 Crate dependencies

Already present per audit (`Cargo.toml` confirmed):
`clap (derive)`, `needletail 0.7`, `rayon 1`, `csv 1`,
`serde 1 (derive)`, `serde_json 1 (float_roundtrip — keep this)`,
`ahash 0.8`, `anyhow 1`, `thiserror 2`, `log 0.4`.

To add:

- `regex = "1"` — for `find_ssrs` in `src/ssr/find_ssrs.rs`. The
  prototype uses pattern `(([gatc]{L})\2{min-1,})` per motif
  length 1..14. A hand-rolled scanner is viable (the pattern is
  trivial) but `regex` keeps the port faithful and is fast enough
  for L ≤ 14. **`[Q‑2]` — accept this dep, or roll our own?**
  Default: accept `regex`. Cost: small (≈1.5 MB binary, 200 kloc
  build); benefit: byte-identical scanning logic and zero risk of
  regex-edge-case drift.

Nothing else.

### 2.4 Concurrency

Rayon is already used in `kite::analyze_records`. Each rule-proto
stage is **per-record-parallel**: every record's analysis is
independent. Per-stage parallelism is identical:

```rust
records.par_iter().map(|rec| run_stage(rec, &cfg)).collect()
```

No cross-record synchronisation; rayon thread pool inherits from
the process-wide setting. Performance §7 estimates assume default
rayon parallelism (#CPUs).

---

## §3 Data exchange format

### 3.1 File names

Each subcommand emits one or two TSVs at `<prefix>.<stage>.tsv`,
matching the prototype's conventions verbatim (so existing eval
scripts at `tools/rule_proto/eval/` continue to work against Rust
output without modification):

| Stage | Files emitted |
|---|---|
| `kite-periodicity --emit-periods` (existing) | `<prefix>.kite.tsv`, `<prefix>.kite.peaks.tsv` |
| `rule-classify` | `<prefix>.verdicts.tsv` |
| `subrepeat-scan` | `<prefix>.subrepeat.tsv`, `<prefix>.windows.tsv` |
| `ssr-scan` | `<prefix>.ssr.tsv`, `<prefix>.ssr.regions.tsv` |
| `hor-validate` | `<prefix>.hor_within_tile.tsv` |
| `summary-merge` | `<prefix>.summary.tsv` |

`analyze <fasta> -o <prefix>` emits the union (8 TSVs) into one
namespace.

### 3.2 Column order and naming

The Rust port emits **the same column order and names** as the
prototype. The Python audit identified the exact orderings; they
are reproduced verbatim in the per-stage sections §5.1 – §5.5 below
so the implementer has them on one page.

### 3.3 Numerical formatting

The prototype's three Python writers use **three different**
float-formatting policies (a pandas accident — but downstream
analysis tooling may have come to depend on them):

| Writer | Float format |
|---|---|
| `rule_proto.py` | `%.6g` (set explicitly) |
| `summary.py` | `%.4g` (set explicitly) |
| `subrepeat_scan.py`, `ssr_scan.py`, `hor_within_tile_check.py` | pandas default (≈ `%g`-style, locale-aware, up to 12 sig digits) |

Two viable approaches:

- **Match per-writer** (proposed default): replicate each writer's
  format string. Highest fidelity, supports byte-diff regression
  testing against current Python output.
- **Unify on one policy** (e.g. `%.6g` everywhere): cleaner, but
  loses byte-equivalence and forces eval scripts to tolerate
  formatting drift.

Plan: **match per-writer**. Use Rust's `format!("{:.6}", x)` family
with a small helper `fmt_g(precision: usize, x: f64) -> String` that
reproduces Python's `%.{p}g` semantics (Rust's default `{:e}/{}`
formatting does **not** match `%g` for trailing-zero suppression —
needs a hand-rolled `%g` shim). **`[Q‑3]` — confirm the
match-per-writer policy or relax it?**

### 3.4 String conventions to preserve

The audit surfaced several non-obvious conventions that must be
preserved verbatim:

- `None` (Python) → empty cell. `float("nan")` → literal `nan`.
  Literal `"NA"` strings in some columns (subrepeat `host_period_bp`,
  ssr `dominant_motif`, within-tile `within_top_period` /
  `density_hint` / `phase_contrast`, etc.) — the audit lists which.
- 1-based-inclusive start, 0-based-exclusive end in
  `ssr.regions.tsv` (TideCluster convention).
- `top_motifs` format: `motif:XX.X%` joined by `;`, single decimal,
  e.g. `"GT:88.4%;AT:5.0%"`.
- `consensus_monomer` format: `;`-joined `KMER(count)` strings,
  **uppercase**, e.g. `"GT(4321);AT(800)"`.
- `blocks` (subrepeat): `;`-joined `start-end` pairs, or `"NA"`.

These all materialise inside `<stage>/io.rs::write_tsv` — small,
testable formatting routines.

---

## §4 Subcommand CLI surface

All six subcommands follow the existing kitehor `clap` derive
pattern (see audit §4 for the three-step recipe: enum variant in
`cli.rs`, match arm in `main.rs`, handler fn).

### 4.1 `rule-classify`

```
kitehor rule-classify <peaks.tsv> -o <prefix>
    [--tol 0.015] [--min-period 20] [--min-cluster-frac 0.01]
    [--k-max 30] [--non-mono-ratio 0.5]
    [--founder-floor 0.1] [--high-k-tile-floor 0.05]
    [--lone-significant-frac 0.1]
    [--dump-clusters <dir>]
```

Reads `<peaks.tsv>` (kite peaks long-format); writes
`<prefix>.verdicts.tsv`. The three FLOOR constants are promoted to
CLI flags here (the prototype hardcodes them).

### 4.2 `subrepeat-scan`

```
kitehor subrepeat-scan <fasta> -o <prefix> --kite-peaks <peaks.tsv>
    [--tol 0.05] [--window-mult-sub 5] [--step-frac 4]
    [--top-n-sub 3] [--top-n-host 10]
    [--sub-floor 0.05] [--window-score-floor 0.3] [--min-run 3]
    [--host-sub-ratio-min 3]
    [--min-window-bp 1000]
```

`--kite-peaks` is required (we no longer shell out to kite — the
global peaks must be supplied; in `analyze` the orchestrator wires
them through). `HOST_SUB_RATIO_MIN` and `MIN_WINDOW_BP` are
promoted from prototype hardcoded constants to CLI flags.

Writes `<prefix>.subrepeat.tsv` + `<prefix>.windows.tsv`.

### 4.3 `ssr-scan`

```
kitehor ssr-scan <fasta> -o <prefix>
    [--kite-peaks <peaks.tsv>]
    [--ssr-flag-threshold-pct 30.0]
    [--consensus-dimer-copies 4]
    [--consensus-dimer-min-bp 30]
    [--consensus-max-monomers 3]
    [--consensus-freq-ratio-min 0.3]
    [--motif-min-reps "1:20,2:9,3:6,4:5,5:5,6:5,7:5,8:5,9:5,10:5,11:5,12:5,13:5,14:5"]
```

`--kite-peaks` optional (raw-fallback when omitted). All five
prototype constants get flags. `--motif-min-reps` is a single
string flag of `length:min_reps` pairs to keep the CLI tidy
(parsed via clap's `value_parser`).

Writes `<prefix>.ssr.tsv` + `<prefix>.ssr.regions.tsv`.

### 4.4 `hor-validate`

```
kitehor hor-validate <fasta> --verdicts <verdicts.tsv> \
    --global-peaks <peaks.tsv> -o <prefix>
    [--min-k-for-density 4]
    [--density-window-tile-frac 3]
    [--min-founder-mult 3] [--min-density-window-bp 200]
    [--max-density-windows 1000]
    [--density-rel-floor 0.2]
    [--phase-fold-bins 10]
    [--density-dup-max 0.35] [--density-hor-min 0.7]
    [--phase-contrast-dup-min 0.4] [--phase-contrast-hor-max 0.15]
    [--period-match-tol 0.02]
    [--max-tile-bp 200000]
```

Every threshold from the prototype is exposed.

Writes `<prefix>.hor_within_tile.tsv`.

### 4.5 `summary-merge`

```
kitehor summary-merge --verdicts <verdicts.tsv> \
    --subrepeat <subrepeat.tsv> --ssr <ssr.tsv> \
    [--within-tile <hor_within_tile.tsv>] \
    -o <prefix> [--pure-ssr-pct-threshold 80.0]
```

Writes `<prefix>.summary.tsv`. Outer-joins on `record_id`,
applies the 8-rule first-match-wins classifier.

### 4.6 `analyze`

```
kitehor analyze <fasta> -o <prefix> [<all sub-flags forwarded>]
```

One-shot orchestrator. Internally runs:

1. `kite-periodicity --classify --emit-periods` style call → in-memory
   `Vec<KiteResult>` + writes `<prefix>.kite.tsv` + `.kite.peaks.tsv`.
2. `rule-classify` → `<prefix>.verdicts.tsv`.
3. In parallel (rayon):
   - `subrepeat-scan` → `<prefix>.subrepeat.tsv` + `.windows.tsv`
   - `ssr-scan` → `<prefix>.ssr.tsv` + `.ssr.regions.tsv`
   - `hor-validate` → `<prefix>.hor_within_tile.tsv`
4. `summary-merge` → `<prefix>.summary.tsv`.

CLI exposure: each stage's flag is forwarded with a stage-prefix,
e.g. `--rule-tol 0.015 --subrepeat-tol 0.05 --ssr-flag-threshold-pct
30.0 --hor-density-rel-floor 0.2`. This keeps the global namespace
unambiguous. **`[Q‑4]` — keep this verbose prefix convention, or
let the user pass a config TOML?**

---

## §5 Per-stage detailed plan

### 5.1 `rule-classify`

**Input**: `<prefix>.kite.peaks.tsv` — required columns
`case_id, rank, period, score2_norm`.

**Output**: `<prefix>.verdicts.tsv`. Column order (verbatim from
prototype audit):

```
case_id  verdict  founder  multiplicity  tile  founder_score
tile_score  confidence  n_clusters  reason
```

**Float format**: `%.6g`. `founder`/`tile` rounded to 2 decimals
before formatting. Empty cells where prototype writes `None`.

**Algorithm** (per case_id, preserving `groupby(sort=False)` order):

1. Drop rows with `period < min_period`.
2. **Cluster** by relative gap on ascending `period`. Cut at
   `(p_cur - p_prev) / p_cur > tol`. (**Divisor is `p_cur`, not
   the cluster mean** — the audit flagged this as easy to misread.)
3. Per cluster: `rep_period = weighted_mean(period, w=score2_norm)`,
   `total_score = sum(score2_norm)`, `n_peaks`, `min_rank`, members.
   Fallback to plain mean when all weights are zero.
4. Drop clusters with `total_score < min_cluster_frac × max_cluster_total`.
5. Sort clusters by `total_score` desc; `top = clusters[0]`.
6. Decision tree (first-match-wins):
   - **Case A**: for each shorter cluster `c`, check
     `k = round(top.period / c.period) ∈ [2, k_max]` AND
     `|top.period − k·c.period| / top.period ≤ tol` AND
     `c.total_score ≥ founder_floor × top.total_score` AND
     `harmonic_check(c.period, k)`. Among qualifiers pick the one
     with **smallest `rep_period`**. Verdict: HOR(founder=c, k, tile=top).
   - **Case B walk** (k = 2, 3, …, k_max): for the cluster at
     `k·top.period` within tol, require non-monotonic bump
     (`score > prev_seen_score`), `≥ high_k_tile_floor × top` (k≥3),
     existence of a tile-doubling cluster at `2k·top` with score>0,
     and harmonic check. First k satisfying all fires. Verdict:
     HOR(founder=top, k, tile=match). For k=2 the
     `non_mono_ratio × top` threshold applies (prototype's specific
     branching at k=2).
   - **Any larger multiple at all** → `simple_tr(top)` reason
     `monotonic_multiples`.
   - **Lone significant cluster**: exactly one cluster has
     `total_score ≥ lone_significant_frac × top.total_score`
     → `simple_tr(top)` reason `lone_significant_cluster`.
   - Else `unresolved`.
7. **Harmonic check** (Cases A and B): require
   `cluster_score_near(2k·founder) ≥ cluster_score_near((k+1)·founder)`,
   where `cluster_score_near(target)` returns the **max
   `total_score` of any cluster within tol of `target`**, or 0 if
   none. Both missing ⇒ passes.
8. `multiplicity=1, tile=founder, tile_score=founder_score` for
   `simple_tr`.
9. `confidence = clip((founder.total_score + tile.total_score) /
   sum_all_total_score, 0, 1)`. For `simple_tr` only the founder
   term contributes; for the `no_clusters` unresolved row
   `confidence` is empty.

**`reason ∈ {top_is_multiple_of_founder, non_monotonic_bump, k2_ratio,
monotonic_multiples, lone_significant_cluster, no_multiples, no_clusters}`**
(values from prototype audit).

**Implementation**:

- `rule_classify::cluster::cluster_peaks(&[PeakRow], tol)
  -> Vec<Cluster>` — pure function, easy to unit-test.
- `rule_classify::decide::decide(&[Cluster], &Config) -> Verdict`
  — pure function, ditto.
- `rule_classify::io::read_peaks` / `::write_verdicts` — TSV I/O.

**Validation**:

- Six band-tolerant fixtures at `tools/rule_proto/fixtures/*.peaks.tsv`
  + `expected.tsv`. Reproduce `test_rule_proto.py`'s logic in Rust
  integration test at `tests/rule_classify_fixtures.rs`. **Required**
  to pass before this stage is considered done.
- Byte-diff against Python verdicts.tsv on a small sim corpus
  (e.g. 50 records from `test_data/sim/`) as an additional gate.
COMMENT - do also walidation again full set docs/new/rule_proto_impl_plan.md - compare to python implentation

**Performance**: trivial — sub-second on test_590.

---

### 5.2 `subrepeat-scan`

**Input**: FASTA + `--kite-peaks <peaks.tsv>`.

**Output A** — `<prefix>.subrepeat.tsv` — column order verbatim:

```
record_id  length_bp  host_period_bp  subrepeat_period_bp
subrepeat_flag  reason  n_windows_total  n_windows_sub
n_windows_non_sub  n_subrepeat_blocks  subrepeat_coverage_bp
subrepeat_coverage_pct  blocks
```

**Output B** — `<prefix>.windows.tsv` — column order verbatim:

```
record_id  window_start  window_end  top_period  top_score2_norm
class_raw  class
```

`subrepeat_flag ∈ {yes, no, none}`,
`reason ∈ {blocks+non_sub, no_blocks, no_non_sub_windows, no_candidate_pair}`,
`class ∈ {sub, non_sub}`.

**Pseudo-NA**: `host_period_bp = "NA"` and `subrepeat_period_bp =
"NA"` and `blocks = "NA"` when no candidate pair.

**Float format**: pandas-default (≈ 12 sig digits). The only
non-integer field is `subrepeat_coverage_pct`, rounded to 2 decimals
before formatting.

**Algorithm** (per record):

1. **Pick candidates** from `--kite-peaks` filtered to this record:
   - `sub_cand` = shortest period in `top_n_sub` by `score2_norm`
     with `score2_norm ≥ sub_floor`.
   - `host_cand` = period at **max `score2_norm`** in `top_n_host`
     with `period ≥ host_sub_ratio_min × sub_cand`. **(Audit note:
     prototype docstring says "longest" but code picks
     strongest-scored — port matches the code.)**
   - If either missing → emit "none" row, no windows file rows.
2. **Window the sequence**:
   - `w = max(round(window_mult_sub × sub_cand), min_window_bp)`.
   - If `w ≥ L`, one window `(0, L)`.
   - Else step `s = w // step_frac`; sliding windows
     `(0, w), (s, s+w), …` plus a final flush `(L-w, L)` if the last
     end < L.
3. **Per-window kite call** — **in-process** via the new
   `kmer_pairs::PairIndex::profile_window` + the refactored
   `kite::find_peaks_with_score`. The per-window top peak's period
   + `score2_norm` go into `windows.tsv`. If no peaks, emit
   `top_period=-1`, `top_score2_norm=0.0`, `class_raw=class="non_sub"`.
4. **Classify window**: `class_raw = "sub"` iff
   `|top_period − sub_cand| / sub_cand ≤ tol` AND
   `top_score2_norm ≥ window_score_floor`. Else `"non_sub"`.
5. **Smoothing** (`smooth_runs`, fixed-point): each iteration find
   the **shortest** run with `length < min_run`; absorb into the
   longer neighbour; ties at interior runs go to the **previous**
   neighbour (`prev_len >= next_len`). Loop until no run is below
   `min_run`. Audit confirms this exact tiebreak rule — match it.
6. **Blocks**: consecutive `class="sub"` smoothed windows → one
   block `[first_start, max(end)]` (handle overlapping windows by
   taking `max` of `end`).
7. **Flag**: `yes` iff `n_blocks ≥ 1 AND n_non_sub ≥ 1`; else `no`.

**`[Q‑1] resolution — per-window peak finding**:

The prototype shells out to the kite binary per window, so the
per-window peak detection uses kite's standard pipeline
(`compute_neighbor_profile` → `compute_background` → `gaussian_smooth`
→ `find_peaks_with_score`) on the **window subsequence**. The Rust
port should reproduce this contract:

- `PairIndex::profile_window(start, end, max_d)` computes the
  H[d] histogram restricted to k-mer pairs whose left member
  starts in `[start, end)`. This matches kite's behaviour on a
  sliced subsequence **only if** the window's k-mer composition is
  the dominant input to the background model.
- Per-window background re-estimation: re-run kite's
  composition-matched FNV-1a-seeded background on the window's
  k-mer counts. This preserves byte-equivalence with the
  subprocess version because the FNV-1a seed is derived from the
  k-mer multiset, which is identical between subprocess and
  in-process invocations.
- Periphery cost: the per-window background dominates the
  per-window total runtime today. If even after fusing the calls
  the per-window background re-estimation is too slow,
  fall back to the global array-wide background passed through.
  **Plan**: implement re-estimation first (correctness); benchmark;
  switch to passing the global background only if needed and only
  with a documented loss of byte-equivalence.

**Validation**:

- Synthetic fixtures at `tools/rule_proto/subrepeat/synthetic.fasta`
  (10 records). The six positive cases must flag `yes`; the four
  negatives must return `no`/`none`. Reproduce the prototype's
  assertions in `tests/subrepeat_fixtures.rs`.
- Byte-diff against Python `subrepeat.tsv` + `windows.tsv` on
  ≥50 records from `test_data/sim/`.
CO


**Performance**: prototype 2:50 wall on test_590, dominated by
459 k subprocess invocations. Rust port target: ≤ 20 s on the same
corpus. The `PairIndex` build is `O(L)` per record; each window
scan is `O(window_kmer_count)`.

---

### 5.3 `ssr-scan`

**Input**: FASTA + optional `--kite-peaks`.

**Output A** — `<prefix>.ssr.tsv` — column order verbatim:

```
record_id  length_bp  ssr_flag  dominant_motif  dominant_motif_length
dominant_motif_repeats  dominant_motif_coverage_pct
total_ssr_coverage_pct  top_motifs  ssr_method  consensus_period_bp
consensus_monomer  ssr_raw_dominant_motif
ssr_raw_dominant_motif_coverage_pct  ssr_raw_total_coverage_pct
ssr_raw_n_regions  ssr_raw_top_motifs
```

**Output B** — `<prefix>.ssr.regions.tsv` — column order verbatim:

```
record_id  ssr_number  motif_length  motif_sequence  repeats
start  end  normalized_motif
```

`ssr_flag ∈ {yes, no}`,
`ssr_method ∈ {raw_fallback, consensus_single, consensus_multi}`.

**Quirks** (audit-confirmed, port verbatim):

- `motif_sequence` lowercase, `normalized_motif` uppercase.
- `start` 1-based inclusive, `end` 0-based exclusive (TideCluster).
- `dominant_motif = "NA"`, `dominant_motif_length = "NA"`,
  `top_motifs = "NA"` when no motifs found.
- `consensus_period_bp = "NA"` when kite has no entry for that
  record; `consensus_monomer = "NA"` when no validated dimer.
- `consensus_monomer` format: `KMER(count);KMER(count);…` uppercase.
- `top_motifs` format: `motif:XX.X%;…` (one decimal).

**Algorithm** (per record):

1. **`find_ssrs(seq)`**:
   - Lowercase the sequence.
   - For each `(L, min_reps) ∈ motif_min_reps_spec` in ascending L:
     - Regex `(([gatc]{L})\2{min-1,})` finditer.
     - Skip if motif itself is homopolymer (handle separately
       as `homopolymers_buf`; emit only if `results` empty
       at end).
     - Skip if `start` already claimed by a hit at shorter L
       (`locations` set, dedupe by start position).
   - Emit `{ssr_number, motif_length, motif_sequence (lowercase),
     repeats, start (1-based), end (0-based exclusive),
     normalized_motif (uppercase canonical)}`.
2. **`normalize_motif`**:
   `min(rotations(motif.upper()) ∪ rotations(reverse_complement.upper()))`.
3. **Per-canonical aggregation**: group hits by `normalized_motif`,
   sum `total_repeats`, sum `total_coverage_bp =
   sum(end - start + 1)`, count `n_regions`.
   **Note**: the `end - start + 1` formula combines 1-based start
   with 0-based-exclusive end so the bp count equals `end - start + 1`
   for a contiguous run. Port matches.
4. **`get_unique_motifs`**: drop any motif that is a `k × shorter`
   repeat (e.g. AT vs ATAT).
5. **Pick dominant motif** = max by `total_coverage_bp`.
   `top_motifs` = top 3 with `motif:pct%`.
6. **Raw summary**: `dominant_motif_coverage_pct =
   100 × dom_cov_bp / length_bp` (rounded 2 decimals).
   `ssr_flag = "yes"` iff `dominant_motif_coverage_pct ≥ 30`.
7. **Consensus path** (when `--kite-peaks` supplied):
   a. From `--kite-peaks` get this record's top-peak period `P`.
   b. `extract_consensus_monomers(seq, P)`: Counter over all
      P-mer sliding windows; walk most-common in order, dedupe by
      `normalize_motif`, break when `count < freq_ratio × top_count`
      or `consensus_max_monomers` collected.
   c. **Validate**: for each kmer build dimer of length
      `max(kmer_len × consensus_dimer_copies, consensus_dimer_min_bp)`;
      run `find_ssrs` on it; keep if `ssr_flag="yes"` in dimer.
   d. **Branch**:
      - 0 validated → `raw_fallback` (use step-6 result).
      - 1 unique canonical → `consensus_single` (use dimer's
        summary as authoritative).
      - 2+ unique canonicals → `consensus_multi` (use
        `build_multimotif_summary` with **raw** per-motif coverages
        for the validated canonicals, sum clamped to ≤ 100).
8. Authoritative columns come from this branch; raw step-6 values
   go into `ssr_raw_*` diagnostic columns.

**Implementation**:

- `ssr::find_ssrs::scan(&[u8], &MotifSpec) -> Vec<Hit>` — direct
  port. Use `regex::Regex` per L (cached, compiled once per
  process).
- `ssr::find_ssrs::normalize_motif(&[u8]) -> String` — pure.
- `ssr::consensus::extract_consensus_monomers(&[u8], P, &Cfg)
  -> Vec<(String, u64)>` — rolling hash over P-mer windows;
  `ahash::AHashMap<u64, u64>` counter; sort by count desc; dedupe.
  **Important**: P-mers containing `N` are skipped (audit-confirmed).
- `ssr::io::write_*` — TSV writers with the quirks above.

**Validation**:

- 10-record synthetic set + 11 hand-confirmed records spread across
  test_590/test1/test2. Reproduce in `tests/ssr_fixtures.rs`.
- Byte-diff vs Python `ssr.tsv` + `ssr.regions.tsv`.

**Performance**: prototype ~7 min on test_590 (slowest stage).
Rust target: ≤ 20 s. Bottleneck is P-mer Counter — rolling hash
gives ~50× speedup; `regex` Rust port gives ~5–10× over Python.

---

### 5.4 `hor-validate`

**Input**: FASTA + `--verdicts <verdicts.tsv>` + `--global-peaks
<peaks.tsv>`.

**Output** — `<prefix>.hor_within_tile.tsv` — column order verbatim:

```
record_id  global_founder_bp  global_tile_bp  global_founder_score
global_tile_score  global_founder_tile_ratio  within_top_period
within_top_score  within_founder_score  within_founder_top_ratio
decision_hint  founder_density  phase_contrast  density_n_windows
density_hint  skip_reason
```

**Float format**: `%.6g`.

**`density_hint ∈ {localized_duplication, spatially_confirms_hor,
ambiguous, insufficient_phase_bins, NA, k_too_low_for_test(k=N)}`**.

**`decision_hint ∈ {strongly_confirms_hor, weakly_confirms_hor,
ambiguous, suggests_within_monomer_duplication, NA}`**, gated by
`within_founder_top_ratio` against `0.5 / 0.2 / 0.05`. Carry through
unchanged.

**`skip_reason ∈ {"", tile_out_of_range, tile_exceeds_array,
kite_no_peaks_in_window}`**.

**Algorithm** (per record where `verdict == "hor"`):

1. **Skip if `k ≤ 3`** → `density_hint = "k_too_low_for_test(k=N)"`,
   no kite calls. (Audit confirmed: `k < min_k_for_density` means
   k ∈ {2, 3} are skipped by default.)
2. **Within-tile check** (first-tile only):
   - Slice `seq[:tile]`. Run kite in-process.
   - `within_top_period, within_top_score` from top peak.
   - `within_founder_score` = sum of `score2_norm` of peaks within
     `period_match_tol` of founder.
   - `within_founder_top_ratio = within_founder_score / top_score`.
   - `decision_hint` per the threshold table.
3. **Density windows**:
   - `w = max(round(tile / density_window_tile_frac),
     min_founder_mult × founder, min_density_window_bp)`.
   - `step = max(round(founder), max(1, (L-w)//max_density_windows + 1))`.
   - Window list: `(0, w), (step, step+w), …` capped at `≤ ~1000`.
4. **Per-window kite call (in-process)**: top peak + `f_score` =
   sum of `score2_norm` of peaks within `period_match_tol` of founder.
   `founder_present = (f_score / top_score) ≥ density_rel_floor`.
5. **Metrics**:
   - `founder_density = #present / #total`.
   - `phase_contrast = max(frac_per_bin) − min(frac_per_bin)`.
     `bin = (window_mid % round(tile)) // max(round(tile/phase_fold_bins), 1)`,
     clamped to `phase_fold_bins − 1`. Bins with zero windows
     excluded; need ≥ 2 populated bins or
     `density_hint = "insufficient_phase_bins"`.
6. **Combined decision** (audit-confirmed order, duplication wins):
   - `localized_duplication` if `density ≤ density_dup_max OR
     contrast ≥ phase_contrast_dup_min`.
   - Else `spatially_confirms_hor` if `density ≥ density_hor_min
     AND contrast ≤ phase_contrast_hor_max`.
   - Else `ambiguous`.
7. **Global columns**: `global_founder_score` and `global_tile_score`
   are **sums** (not max) of `score2_norm` of all global peaks
   within `period_match_tol` of founder/tile.

**Implementation**:

- `hor_validate::within_tile::run` — runs kite on `seq[:tile]`
  using shared `PairIndex` + kite peak finder.
- `hor_validate::density::run` — windowed sweep using
  `PairIndex::profile_window`.
- Phase folding: small, pure function over a `Vec<(midpoint, present)>`.

**Validation**:

- 11 user-confirmed `tr_with_subrepeat` records must classify as
  `localized_duplication` (per onboarding §5.4 and README's
  per-corpus tables). Add as `tests/hor_validate_duplications.rs`.
- Byte-diff vs Python on test_590 (or sim_v2) for ≥ 100 records.

**Performance**: prototype 25 s on test_590 (59 k subprocess kite
calls). Rust target: ≤ 5 s. Single-record cost = `O(L)` for the
`PairIndex` build + `O(n_windows · window_len)` for scans.

---

### 5.5 `summary-merge`

**Input**:
- `--verdicts` (rule-classify output; `case_id` renamed to
  `record_id`)
- `--subrepeat` (subrepeat-scan output)
- `--ssr` (ssr-scan output)
- `--within-tile` (optional)

**Output** — `<prefix>.summary.tsv` — column order verbatim
(32 columns):

```
record_id  length_bp  hor_verdict  hor_founder  hor_multiplicity
hor_tile  hor_confidence  subrepeat_flag  subrepeat_host_period_bp
subrepeat_period_bp  n_subrepeat_blocks  subrepeat_coverage_pct
ssr_flag  ssr_dominant_motif  ssr_dominant_motif_length
ssr_dominant_motif_repeats  ssr_dominant_motif_coverage_pct
ssr_total_coverage_pct  ssr_top_motifs  ssr_method
consensus_period_bp  consensus_monomer  ssr_raw_dominant_motif
ssr_raw_dominant_motif_coverage_pct  ssr_raw_total_coverage_pct
ssr_raw_n_regions  ssr_raw_top_motifs  density_hint
founder_density  phase_contrast  density_n_windows  combined_class
```

**Float format**: `%.4g`.

**Algorithm**:

1. Outer-join verdicts + subrepeat + ssr on `record_id`. Left-join
   within-tile if present.
2. Defaults: `hor_verdict ← "unresolved"`, `subrepeat_flag ← "none"`,
   `ssr_flag ← "no"`, `ssr_dominant_motif_coverage_pct ← 0.0`,
   `density_hint ← ""`.
3. `combined_class` first-match-wins:
   1. `pure_ssr` if `ssr_flag=="yes" AND ssr_dom_pct ≥ pure_ssr_pct_threshold`
   2. `tr_with_nested_tr` if `subrepeat_flag=="yes"`
   3. `tr_with_subrepeat` if `hor_verdict=="hor" AND density_hint=="localized_duplication"`
   4. `hor_with_ssr` if `hor_verdict=="hor" AND ssr_flag=="yes"`
   5. `hor` if `hor_verdict=="hor"`
   6. `tr_with_ssr` if `hor_verdict=="simple_tr" AND ssr_flag=="yes"`
   7. `tr` if `hor_verdict=="simple_tr"`
   8. `unresolved` otherwise

**Implementation**:

- `summary::join::merge(&[VerdictRow], &[SubrepRow], &[SsrRow],
  Option<&[WithinTileRow]>) -> Vec<SummaryRow>` — pure, fast.
- Use `BTreeMap<String, SummaryRow>` for the outer-join (sorted by
  `record_id` matches pandas merge ordering when `sort=True`; if
  prototype is `sort=False`, use `IndexMap` instead — verify
  against prototype audit's note on join ordering).

**Validation**:

- Exact per-class counts on test_590 (within ±1 per the onboarding
  doc's tolerance) for the **whole** pipeline output. Best run
  alongside §5.6.

**Performance**: trivial, sub-second.

---

### 5.6 `analyze` (orchestrator)

```rust
fn run_analyze(args: AnalyzeArgs) -> Result<()> {
    let records = io::load_fasta(&args.fasta, args.load_qc)?;

    // Stage 1: kite peaks (writes <prefix>.kite.tsv, .kite.peaks.tsv)
    let kite_results = kite::analyze_records(&records, &args.kite_cfg);
    write_kite_outputs(&args.prefix, &kite_results)?;

    // Stage 2: rule-classify (writes <prefix>.verdicts.tsv)
    let verdicts = rule_classify::run(&kite_results, &args.rule_cfg)?;
    rule_classify::io::write(&args.prefix, &verdicts)?;

    // Stages 3a, 3b, 3c in parallel (independent inputs)
    let (subrep_out, ssr_out, hvt_out) = rayon::join3(
        || subrepeat::run(&records, &kite_results, &args.subrep_cfg),
        || ssr::run(&records, &kite_results, &args.ssr_cfg),
        || hor_validate::run(&records, &kite_results, &verdicts,
                             &args.hvt_cfg),
    );
    subrep_out?; ssr_out?; hvt_out?;  // and write their TSVs

    // Stage 4: summary-merge
    summary::run_and_write(&args.prefix, ...)?;
    Ok(())
}
```

In-process data flow eliminates round-trip TSV serialisation
between stages; each stage **still writes its TSV** (for
debugging), but the next stage receives the in-memory typed
struct rather than re-parsing the TSV. This is the second major
perf win after the kite-fusion.

**`[Q‑5] — Should `analyze` also support a "skip-stage" flag**
(e.g. `--no-ssr-scan`) for partial reruns? Default: no, since each
stage is already a standalone subcommand for that purpose.

---

## §6 Migration & cleanup of existing code

### 6.1 Delete or rewrite — `src/rule.rs`

User-confirmed: replace, do not retain. Plan:

1. Introduce `src/rule_classify/` (the new classifier).
2. Update `src/emit_periods.rs` to consume the **new** verdict type
   (rename `RuleVerdict` → `LegacyRuleVerdict` short-term, or
   refactor the bridge to a trait — see §6.2). The
   `EmitPeriodsRow::from_verdict` mapping (founder/tile/k +
   "no signal" + "no HOR" + "unresolved") is preserved.
3. Update `src/detect/mod.rs::auto_periods` to call the new
   classifier in place of the legacy one. The score-mapping table
   in `emit_periods.rs:8-19` stays.
4. Update `src/main.rs::run_kite_periodicity` to dispatch the new
   classifier under `--classify`.
5. Delete `src/rule.rs`. Update `src/lib.rs` to remove the
   `pub mod rule` line (and the doc comments that name it).

### 6.2 Verdict-type compatibility — `src/emit_periods.rs`

`emit_periods.rs::build_rows` takes `Option<&RuleVerdict>`. The new
classifier's natural verdict enum is richer (carries `reason`,
`confidence`, etc.). Two options:

- **(A)** Add a small `From<&NewVerdict> for LegacyShape` adapter
  inside `emit_periods` so the bridge module's API doesn't change.
  Minimal churn.
- **(B)** Refactor `emit_periods::build_rows` to accept the new
  verdict directly and update the v2 detector's auto-periods path
  to feed it. Cleaner long-term.

Plan: **(A)** during port. Schedule **(B)** as a follow-up
once the new pipeline is stable.

### 6.3 Legacy ML pipeline — `[Q‑6]`

The audit identifies multiple legacy-ML modules:

- `src/classifier.rs` (366 LOC) — RF loader.
- `src/classify.rs` (314 LOC) — ML verdict orchestrator.
- `src/features.rs` (441 LOC) — feature builder.
- `src/hor_call.rs` (535 LOC) — independent rule layer used by
  `--no-hor-call` / `--hor-qmax` flags.
- `src/monomer_model.rs` (294 LOC) — `probe_period` cosine identity
  helper.

**`[Q‑6]` — should these be deleted as part of this port, or
left intact for now?** The onboarding doc says the ML path is
"over-sensitive on real centromeric arrays and under-sensitive on
real HORs … use only when the input is drawn from a similar
distribution to the synthetic training corpus." If the new
pipeline supersedes ML in practice, deleting all of `classifier.rs`,
`classify.rs`, `features.rs`, `hor_call.rs` removes ~1700 LOC of
dead code. `monomer_model.rs::probe_period` is independent — keep
unless audited.

**Default proposal** (subject to user input):
- **Keep**: `monomer_model.rs` (used by `probe_period`).
- **Delete in this port**: `rule.rs`.
- **Defer deletion** (until user confirms ML path is unused):
  `classifier.rs`, `classify.rs`, `features.rs`, `hor_call.rs`.
  These would be marked `#[deprecated]` and removed in a follow-up.

### 6.4 CLI cleanup — `src/cli.rs`

The legacy ML/rule flags on `KitePeriodicityArgs`
(`--use-ml-classifier`, `--no-hor-call`, `--hor-qmax`,
`--hor-min-family-share`, etc.) stay until §6.3 resolves. The new
`--classify` semantics: by default dispatch to the new classifier;
keep `--use-ml-classifier` as an opt-in to the legacy ML path
during the deprecation window.

### 6.5 Doc updates

- `CLAUDE.md` (kitehor) — update "What this is" §1 to describe the
  new pipeline. Add `analyze` to the workflow section.
- `docs/rule.md` — supersede with `docs/rule_proto.md` describing
  the new classifier; keep old `rule.md` as `docs/archive/`.
- `README.md` (if present) — point at `docs/rule_proto.md`.

---

## §7 Performance targets

End-to-end on test_590 (2779 records, ~280 MB FASTA), single
machine, default rayon parallelism. Targets are upper bounds —
the port should beat them.

| Stage | Python wall | Rust target | Rationale |
|---|---:|---:|---|
| kite-periodicity | 0:15 | 0:15 | already Rust |
| rule-classify | 0:10 | < 0:01 | trivial |
| subrepeat-scan | 2:50 | < 0:30 | kite-fusion ~50× |
| ssr-scan | 7:00 | < 0:20 | rolling hash + Rust regex ~20–50× |
| hor-validate | 0:25 | < 0:05 | kite-fusion ~100× |
| summary-merge | 0:01 | < 0:01 | trivial |
| **analyze (end-to-end)** | **~10:40** | **< 1:30** | |

Acceptance gate: full analyze on test_590 in under 90 s on the
project workstation. If we beat 60 s, great.

---

## §8 Validation strategy

### 8.1 Tier 1 — fixture-equivalent (mandatory)

| Stage | Fixture | Tolerance | Test file |
|---|---|---|---|
| rule-classify | `tools/rule_proto/fixtures/*.peaks.tsv` (6 records) | exact `verdict`, banded founder/tile from `expected.tsv` | `tests/rule_classify_fixtures.rs` |
| subrepeat-scan | `tools/rule_proto/subrepeat/synthetic.fasta` (10 records) | flag + block count per prototype | `tests/subrepeat_fixtures.rs` |
| ssr-scan | 10-record synthetic + 11 hand-confirmed | canonical motif + dominant pct exact | `tests/ssr_fixtures.rs` |
| hor-validate | 11 confirmed `tr_with_subrepeat` records | all classify as `localized_duplication` | `tests/hor_validate_fixtures.rs` |

These are **blocking** for each stage's PR.

### 8.2 Tier 2 — byte-diff vs Python (per-stage)

For each stage, dump Python output and Rust output side-by-side on
50–100 records sampled from `test_data/sim/` and `test_data/sim_v2/`.
Compare with `diff`. The diff should be empty modulo documented
exceptions (the floating-point formatting policy decisions in §3.3).

Helper script: `tools/rule_proto/diff_outputs.py` — wraps `diff`
with per-column tolerance. Reusable; lives alongside the prototype.

### 8.3 Tier 3 — per-class count match (whole pipeline)

The onboarding doc gives the reference per-class counts for
test1 (155), test2 (165), test_590 (2779). The Rust `analyze`
output must match these **within ±1 per class** as a release gate.

```
test_590 reference: hor=173 tr_with_nested_tr=703
                    tr_with_subrepeat=25 hor_with_ssr=0
                    pure_ssr=110 tr_with_ssr=22 tr=1529
                    unresolved=217
```

### 8.4 Tier 4 — wider corpus eval (informational)

Run the existing `tools/rule_proto/eval/*.py` scripts against the
Rust outputs. Treat the numbers as "should be very close" but
don't gate releases on micro-differences. **`[Q‑7]` — what's the
right size for a continuous-integration corpus?** test_590's
280 MB is too big for CI; the 6-fixture set is too small. Propose
a curated 200-record subset.

---

## §9 Phased rollout

Each phase is a separate PR; each phase ends with green CI + a
commit-message-style milestone tag.

| Phase | Deliverable | Gating tests |
|---|---|---|
| **P0** | Refactor: extract `kmer_pairs.rs`, promote kite internals (no behaviour change). | All existing tests still pass. |
| **P1** | `rule-classify` subcommand + replacement of `rule.rs`. Update `emit_periods.rs` + `detect/auto_periods`. | Tier 1 rule-classify + existing detect/auto-mode tests. |
| **P2** | `summary-merge` subcommand. | Smoke test on stubbed inputs. |
| **P3** | `ssr-scan` subcommand. | Tier 1 ssr fixtures + Tier 2 byte-diff. |
| **P4** | `subrepeat-scan` + `hor-validate` (paired — both share the kite-fusion infrastructure). | Tier 1 fixtures + Tier 2 byte-diff. |
| **P5** | `analyze` orchestrator. | Tier 3 per-class match on test_590. |
| **P6** | Deprecate / delete legacy ML modules per §6.3. | All prior tests still pass; no `--use-ml-classifier` callers in tree. |
| **P7** | Doc updates (`CLAUDE.md`, `docs/rule_proto.md`). | Manual review. |

P3 before P4 because: `summary-merge` is trivial and validates the
joint output shape early; `ssr-scan` is the longest single port
and benefits from early start; `subrepeat-scan` + `hor-validate`
share the `kmer_pairs::PairIndex::profile_window` API and are
better as a pair.

---

## §10 Open design notes (decided in this draft but worth flagging)

1. **`kmer_pairs.rs` is a new top-level module** (rather than a
   submodule of `kite`). Rationale: it'll be consumed by
   `subrepeat`, `hor_validate`, **and** the future detect-side
   block-native-width computation in M8 (per
   `docs/new/detect_m8_plan.md`). A neutral home prevents
   `kite::peaks` from becoming a kitchen sink.

2. **No serde-derive on the row structs initially.** The existing
   crate uses hand-formatted `writeln!` / `write!` for TSV output
   (per audit §6). The new pipeline matches that style for the
   three writers whose float-formatting needs are non-trivial
   (rule-classify `%.6g`, summary `%.4g`, custom `NA`/`nan`
   strings). serde-csv may be revisited later if it can model the
   per-writer quirks cleanly.

3. **Float-formatting helper** (`fmt_g`) is shared across writers.
   Place at `src/util/fmt.rs` or under `src/io_helpers.rs`. Unit
   test against Python output: emit a vector of edge-case floats
   (0, NaN, ±inf, very small, very large, integer-valued) in
   `%.4g` and `%.6g`, compare byte-for-byte with
   `python -c "for x in xs: print('%.6g' % x)"`.

4. **Verdict-type ABI for `emit_periods`** (§6.2) — pick option
   (A) for minimal churn during the port; option (B) is a
   follow-up.

5. **`hor_call.rs` is left intact** for now (it's tied to
   independent CLI flags). Marked for review at §6.3.

---

## §11 Answered open questions

All defaults stand unless overridden. Final decisions:

| Q | Decision |
|---|---|
| **Q‑1** | **Byte-equivalent**: per-window background re-estimation (matches prototype subprocess calls). |
| **Q‑2** | **Accept `regex`** crate dep for `find_ssrs`. |
| **Q‑3** | **Replicate per-writer** float formatting (rule-classify `%.6g`, summary `%.4g`, others pandas-default). |
| **Q‑4** | **Stage-prefixed flags** on `analyze`. TOML config deferred to v2. |
| **Q‑5** | **No** skip-stage flag on `analyze`. Use individual subcommands. |
| **Q‑6** | **DELETE** legacy ML modules in this port: `classifier.rs`, `classify.rs`, `features.rs`, `hor_call.rs` (+ `models/`). `monomer_model.rs::probe_period` survives (independent). Phase P6 of the rollout. |
| **Q‑7** | **Curated small diverse corpus** (replaces 200-record proposal). Selection method: run the prototype against test_590, pick a small set covering each `combined_class` value plus a few edge cases (HOR with marginal scores, simple_tr with lone-significant cluster, etc.). Committed under `test_data/ci_corpus/` with a `manifest.tsv` documenting provenance. Implemented in P7. |
| **Q‑8** | **Mirror** window-id-derived FNV-1a seed (e.g. `{rec_id}__SP_{s}_{e}`). Preserves byte-equivalence with prototype subprocess background. |
| **Q‑9** | **String flag** `--motif-min-reps "1:20,2:9,…"` is OK. |
| **Q‑10** | **Always emit all TSVs** from `analyze`. No `--emit-all-stages` flag. TSV-per-stage is the contract. |

Additional inline annotations from v1 review:

- **Rule-classify validation (§5.1)** — extend Tier 2 to a **full-set
  diff against Python**, not just a 50-record sample. The Python
  prototype is the reference oracle.

---

## §12 Suggested commit cadence

One PR per phase. Each PR description includes:
- The Tier-1 fixture diff (`diff -u prototype.tsv rust.tsv`,
  empty after the port).
- Tier-2 byte-diff summary line.
- A `cargo test --release` invocation log.

Final phase (P5: `analyze`) PR includes the Tier-3 per-class table
side-by-side with the onboarding §6 reference values.

---

## Appendix A — Prototype constants → CLI flags table

Quick reference. Default values from the prototype.

| Constant (Python) | CLI flag | Default | Stage |
|---|---|---|---|
| `FOUNDER_FLOOR` | `--founder-floor` | 0.1 | rule-classify |
| `HIGH_K_TILE_FLOOR` | `--high-k-tile-floor` | 0.05 | rule-classify |
| `LONE_SIGNIFICANT_FRAC` | `--lone-significant-frac` | 0.1 | rule-classify |
| `tol` | `--tol` | 0.015 | rule-classify |
| `min_period` | `--min-period` | 20 | rule-classify |
| `min_cluster_frac` | `--min-cluster-frac` | 0.01 | rule-classify |
| `k_max` | `--k-max` | 30 | rule-classify |
| `non_mono_ratio` | `--non-mono-ratio` | 0.5 | rule-classify |
| `HOST_SUB_RATIO_MIN` | `--host-sub-ratio-min` | 3 | subrepeat-scan |
| `MIN_WINDOW_BP` | `--min-window-bp` | 1000 | subrepeat-scan |
| `DEFAULT_MIN_RUN_WINDOWS` | `--min-run` | 3 | subrepeat-scan |
| `tol` | `--tol` | 0.05 | subrepeat-scan |
| `window_mult_sub` | `--window-mult-sub` | 5 | subrepeat-scan |
| `step_frac` | `--step-frac` | 4 | subrepeat-scan |
| `top_n_sub` | `--top-n-sub` | 3 | subrepeat-scan |
| `top_n_host` | `--top-n-host` | 10 | subrepeat-scan |
| `sub_floor` | `--sub-floor` | 0.05 | subrepeat-scan |
| `window_score_floor` | `--window-score-floor` | 0.3 | subrepeat-scan |
| `SSR_FLAG_THRESHOLD_PCT` | `--ssr-flag-threshold-pct` | 30.0 | ssr-scan |
| `DEFAULT_SPECS` | `--motif-min-reps` | "1:20,2:9,3:6,4:5,…,14:5" | ssr-scan |
| `CONSENSUS_DIMER_COPIES` | `--consensus-dimer-copies` | 4 | ssr-scan |
| `CONSENSUS_DIMER_MIN_BP` | `--consensus-dimer-min-bp` | 30 | ssr-scan |
| `CONSENSUS_MAX_MONOMERS` | `--consensus-max-monomers` | 3 | ssr-scan |
| `CONSENSUS_FREQ_RATIO_MIN` | `--consensus-freq-ratio-min` | 0.3 | ssr-scan |
| `MIN_K_FOR_DENSITY` | `--min-k-for-density` | 4 | hor-validate |
| `DENSITY_WINDOW_TILE_FRAC` | `--density-window-tile-frac` | 3 | hor-validate |
| `MIN_FOUNDER_MULT` | `--min-founder-mult` | 3 | hor-validate |
| `MIN_DENSITY_WINDOW_BP` | `--min-density-window-bp` | 200 | hor-validate |
| `MAX_DENSITY_WINDOWS` | `--max-density-windows` | 1000 | hor-validate |
| `DENSITY_REL_FLOOR` | `--density-rel-floor` | 0.2 | hor-validate |
| `PHASE_FOLD_BINS` | `--phase-fold-bins` | 10 | hor-validate |
| `DENSITY_DUP_MAX` | `--density-dup-max` | 0.35 | hor-validate |
| `DENSITY_HOR_MIN` | `--density-hor-min` | 0.7 | hor-validate |
| `PHASE_CONTRAST_DUP_MIN` | `--phase-contrast-dup-min` | 0.4 | hor-validate |
| `PHASE_CONTRAST_HOR_MAX` | `--phase-contrast-hor-max` | 0.15 | hor-validate |
| `PERIOD_MATCH_TOL` | `--period-match-tol` | 0.02 | hor-validate |
| `PURE_SSR_PCT_THRESHOLD` | `--pure-ssr-pct-threshold` | 80.0 | summary-merge |

---

## Appendix B — Subcommand-to-Python-script mapping (one-line cheat sheet)

| Rust subcommand | Python script | Source LoC |
|---|---|---:|
| `kite-periodicity` (existing) | n/a (already Rust) | — |
| `rule-classify` | `tools/rule_proto/rule_proto.py` | 372 |
| `subrepeat-scan` | `tools/rule_proto/subrepeat_scan.py` | 443 |
| `ssr-scan` | `tools/rule_proto/ssr_scan.py` | 527 |
| `hor-validate` | `tools/rule_proto/hor_within_tile_check.py` | 469 |
| `summary-merge` | `tools/rule_proto/summary.py` | 226 |
| `analyze` | (none — new orchestrator) | — |

Estimated Rust LoC, including tests: ~3500–4500 (Python prototype
total ≈ 2040 LoC).

---

*End of v1 draft. Please review and answer §11 — I'll fold the
answers into v2 before kicking off P0.*

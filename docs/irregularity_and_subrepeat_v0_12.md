# Irregularity detection + subrepeat gate (v0.12)

This document describes two related improvements to kitehor's structural
classification of tandem repeats, made during the v0.11 → v0.12 cycle:

1. **`tr_with_subrepeat` density gate** — already shipped in Rust
   (commit `bdaafe4`). Reduces false-positive `tr_with_subrepeat` calls
   in the summary cascade from 471 → 12 on the 3024-record IPIP corpus.
2. **Irregularity detection prototype** — Python-only at
   `tools/rule_proto/irregularity.py`. Quantifies indel-like
   phase-shift events in a tandem array using register-locked k-mers
   as phase markers. Standalone for now; planned for Rust port and
   later integration into `analyze`.

Both share a common motivation: distinguishing **structural** features
of a tandem array (real subrepeats, real indels) from **localized
events** (indel artifacts, natural per-monomer length drift) that the
existing pipeline was conflating.

Audience: someone picking up further development or optimization of
either stage. Assumes familiarity with `kitehor`'s pipeline
(`kite-periodicity → rule-classify → tandem-validate → summary-merge`)
and the existing documentation under `docs/`.

---

## Table of contents

1. [Background: pipeline recap](#background-pipeline-recap)
2. [Part 1 — `tr_with_subrepeat` density gate](#part-1--tr_with_subrepeat-density-gate)
   1. [Problem](#11-problem)
   2. [Fix](#12-fix)
   3. [Implementation](#13-implementation)
   4. [Validation & corpus impact](#14-validation--corpus-impact)
   5. [Tunables and how to override](#15-tunables-and-how-to-override)
   6. [Future work for this gate](#16-future-work-for-this-gate)
3. [Part 2 — Irregularity prototype (v1, phase-mod-P)](#part-2--irregularity-prototype)
   1. [Motivation: which signal we're after](#21-motivation-which-signal-were-after)
   2. [Algorithm overview](#22-algorithm-overview)
   3. [Algorithm step-by-step](#23-algorithm-step-by-step)
   4. [Output schema](#24-output-schema)
   5. [Tunables and their defaults](#25-tunables-and-their-defaults)
   6. [Calibration history](#26-calibration-history)
   7. [Known limitations](#27-known-limitations)
   8. [Performance notes](#28-performance-notes)
   9. [Files of interest](#29-files-of-interest)
   10. [Rust port plan](#210-rust-port-plan)
4. [Part 3 — Irregularity v2 (distance-residual + dropout/indel split)](#part-3--irregularity-v2-distance-residual)
   1. [Why a second prototype](#31-why-a-second-prototype)
   2. [Core algorithm](#32-core-algorithm)
   3. [Phase-bin clustering (the load-bearing fix)](#33-phase-bin-clustering)
   4. [Output schema](#34-v2-output-schema)
   5. [`ground_truth3` calibration corpus](#35-ground_truth3-calibration-corpus)
   6. [v1 vs v2 head-to-head](#36-v1-vs-v2-head-to-head)
   7. [Known limitations of v2](#37-known-limitations-of-v2)
   8. [Files of interest](#38-files-of-interest-v2)
5. [Glossary](#glossary)

---

## Background: pipeline recap

The current v0.12 pipeline on a FASTA input:

```
FASTA
  │
  ▼  kite-periodicity        (src/kite.rs)
  ├──→ kite.tsv              one row per record: rank-1/2/3 monomer_size
  └──→ kite.peaks.tsv        long format: every kept peak per record
       │
       ▼  rule-classify       (src/rule_classify/)
       └──→ verdicts.tsv     hor / simple_tr / unresolved, with founder,
                              multiplicity, tile, confidence, n_clusters
       │
       ├───→ tandem-validate  (src/tandem_validate/)
       │     emits tandem_validate.tsv → tv_density, tv_decision …
       │
       └───→ ssr-scan         (src/ssr/)
             emits ssr.tsv / ssr.regions.tsv → raw + consensus SSR coverage
       │
       ▼  summary-merge       (src/summary/)
       └──→ summary.tsv      joins everything; produces combined_class
                              from the 9-rule cascade
```

The **irregularity** scan is a *new, parallel* stage at the same level
as tandem-validate / ssr-scan: it consumes the FASTA + the kite top
periods and emits per-record irregularity metrics. Not yet wired into
`analyze` — runs standalone.

---

## Part 1 — `tr_with_subrepeat` density gate

### 1.1 Problem

Prior to v0.12, the summary cascade promoted **any** record with
`tv_decision == "localized_subrepeat"` to `tr_with_subrepeat` (or its
SSR partner). No check on **how much** of the array actually carried
the subrepeat signal — just "did some window hit a subrepeat".

This conflated two structurally distinct cases:

| pattern | `tv_density` typical | what it really is |
|---|---|---|
| Real `tr_with_subrepeat` (subrepeat present in every monomer at a phase-locked offset) | ≥0.7 | a structural property of the monomer — e.g. HSAT-1's 17-bp subrepeat in a 170-bp monomer |
| Indel / localized duplication artifact (a few windows of the array transiently look subrepeat-like due to a single event) | ≪0.5 | a rare local event, not array-wide structure |

On the 3024-record IPIP corpus, of the 473 records flagged
`localized_subrepeat`:

| `tv_density` bracket | n | % |
|---|--:|--:|
| < 0.1 | 297 | 63% |
| 0.1 – 0.5 | 149 | 31% |
| ≥ 0.5 | 27 | 6% |
| ≥ 0.7 (the new gate) | **12** | **2.5%** |

**~94 % of pre-v0.12 `tr_with_subrepeat` calls were misfires** on
records where only a handful of windows showed any subrepeat
signature.

### 1.2 Fix

Add a density gate to the cascade: a `localized_subrepeat` decision
hint must additionally have `tv_density ≥ subrepeat_density_min`
(default **0.4** as of v0.12.1; was 0.7 in v0.12) before it promotes
to `tr_with_subrepeat`. Otherwise the record falls through to its
natural class (`tr` / `unresolved` / `hor` / `*_with_ssr`).

Surface the underlying density as a new summary column,
`subrepeat_coverage_pct = tv_density × 100`, so the reader can see the
gate's input directly per record.

The default was **calibrated against the `ground_truth3` simulator's
`subcompound` stratum** (200 records with controlled compound-monomer
structure):

| `density_min` | true compound recall | FP on 1100 non-subrepeat controls |
|---|--:|--:|
| 0.7 (initial v0.12 default) | 44.5% | 0% |
| 0.6 | 60.0% | 0% |
| 0.5 | 75.5% | 0% |
| **0.4 (current)** | **80.5%** | **0%** |
| 0.3 | 81.0% | 0.18% (2/1100) |

0.4 nearly doubles recall vs 0.7 while still holding zero false
positives on a 1100-record control set. 0.3 produces one extra recall
point at the cost of 2 false positives, so 0.4 is the right balance.
On the IPIP corpus the v0.12.1 default produces 34 `tr_with_subrepeat`
calls (vs 12 under v0.12).

### 1.3 Implementation

Changes are entirely in `src/summary/`:

- **`src/summary/mod.rs`**:
  - `Config.subrepeat_density_min: f64` added with default 0.7.
  - `combined_class()` gains a `tv_density: Option<f64>` parameter.
    Treats `None` as failing the gate (safe default — without density
    data we can't certify the promotion).
  - Cascade order unchanged; the only logic change is the new conjunct
    `tv_density.is_some_and(|d| d >= cfg.subrepeat_density_min)` inside
    the `localized_subrepeat` branch.
- **`src/summary/io.rs`**:
  - Parses the `tv_density` column from `tandem_validate.tsv` and
    threads it into `combined_class()`.
  - Adds a `subrepeat_coverage_pct` column to `TARGET_COLUMNS` and
    `col_fmt()`. Always emitted (empty when `tv_density` is missing).
- **`src/cli.rs`**:
  - `SummaryMergeArgs` gains `--subrepeat-density-min` (default 0.7).
  - `AnalyzeArgs` gains the same as `--subrepeat-density-min`.
- **`src/main.rs`**: wires both flags into `summary::Config`.
- **Tests**: 5 new tests in `summary/mod.rs::tests` covering the gate's
  positive / negative / boundary / missing-density cases. The pre-existing
  smoke fixture (`tests/summary_merge_smoke.rs`) was bumped to
  `tv_density=0.8` so it still exercises the positive `tr_with_subrepeat`
  path.
- **CI manifest** (`test_data/ci_corpus/manifest.tsv`): 3 partial/short/
  small fixtures relabeled from `tr_with_subrepeat` to `unresolved`
  with notes documenting the gate-induced demotion.

### 1.4 Validation & corpus impact

Verified byte-equivalence with the Python prototype on both ground_truth
and IPIP corpora. Re-running `summary-merge` on the IPIP intermediates
(verdicts + tandem_validate + ssr) with the new code:

| class | baseline (pre v0.12) | new (gate=0.7) | Δ |
|---|--:|--:|--:|
| `pure_ssr` | 46 | 46 | 0 |
| `tr_with_subrepeat` | **471** | **12** | **−459** |
| `tr_with_subrepeat_with_ssr` | 2 | 0 | −2 |
| `hor` | 207 | **279** | +72 |
| `hor_with_ssr` | 3 | 5 | +2 |
| `tr` | 1950 | 2263 | +313 |
| `tr_with_ssr` | 35 | 35 | 0 |
| `unresolved` | 232 | 306 | +74 |
| `unresolved_with_ssr` | 78 | 78 | 0 |

Two effects, both in the right direction:

- **−459 demoted `tr_with_subrepeat`**: matches the p97.5 of the
  `tv_density` distribution for that class. By design.
- **+72 `hor`**: records where `hor_verdict == hor` but a *low-density*
  `localized_subrepeat` was cascade-shadowing them. The gate removes
  the shadowing.

### 1.5 Tunables and how to override

```bash
# rule_classify standalone
kitehor summary-merge --subrepeat-density-min 0.5 ...

# end-to-end analyze
kitehor analyze --subrepeat-density-min 0.5 ...
```

At `--subrepeat-density-min 0` the gate is effectively off (reproduces
pre-v0.12 behaviour modulo float-rounding noise).

The 0.4 default (v0.12.1) was chosen from the `ground_truth3`
calibration table in §1.2 — highest recall point where FP rate stays
at 0% on 1100 simulated controls. The original 0.7 default was
chosen from the IPIP density distribution alone (without simulated
controls) and turned out to be too strict on real compound subrepeats
where the internal tile fills only part of the monomer (e.g. HSAT-1's
17-bp repeat in 170-bp monomer where the simulator's region_frac
analogue is ~0.5-0.7, not 1.0).

### 1.6 Future work for this gate

- **A `_with_local_event` class** if you want to surface
  low-`tv_density` records as a distinct diagnostic rather than mixing
  them into plain `tr`. Currently the demoted records lose their
  subrepeat hint entirely from the combined_class — only the
  `tv_decision`, `tv_density`, and new `subrepeat_coverage_pct` columns
  carry it.
- **Phase coherence test inside `tandem-validate`**: currently
  `tv_decision == localized_subrepeat` doesn't check whether the
  subrepeat appears at the *same phase* across positive windows.
  Real subrepeats are phase-locked; indel artifacts are not.
  Adding a phase-coherence column there would let the cascade
  distinguish "real subrepeat at moderate coverage" from "patchy
  artifact at moderate coverage".

---

## Part 2 — Irregularity prototype

### 2.1 Motivation: which signal we're after

The kite + rule-classify + tandem-validate stages tell us what the
*structural period* of an array is and whether the array carries a
subrepeat or HOR ladder. They don't tell us **whether the period itself
is consistently maintained across the array**.

In real biology, indel events (insertions, deletions, small
duplications) at scales below a full monomer disrupt the tandem
register without destroying it. Two classes of disruption matter:

- **Phase-shifting indels** — an insertion of δ bp shifts the *register*
  of every monomer downstream by δ. Multiple such events accumulate
  along the array.
- **Whole-monomer copy gain/loss** — adding or removing exactly one full
  monomer (δ = P) leaves the downstream register unchanged modulo P
  and is **invisible** to this method. (This is a known limitation —
  a separate "expected vs observed copy count" metric would be needed.)

So the question we answer: "**how many phase-shifting indel-like events
are present in this array, of what sizes, and how unevenly are they
distributed?**".

The signal we use: **register-locked k-mers as phase markers**. A k-mer
that occurs at a stable phase within each monomer is a marker for that
position. If many such markers all shift by the same δ at a given point
in the array, that's an indel of size δ.

A substitution that knocks out one k-mer occurrence affects that single
marker only — the rest stay in register. So substitutions and indels
are naturally separable: substitutions raise the `dropout_rate`, indels
raise `step_count` / `phase_shift_burden_pct`.

### 2.2 Algorithm overview

```
For each FASTA record:
   1. Read kite-supplied period P.
   2. Extract top-100 6-mers + their position lists from the sequence.
   3. Per k-mer:
      a. Greedily pick up to 5 anchor phases (densest 5-bp windows over
         positions mod P).
      b. Filter: keep k-mers whose WINDOWED anchor coverage (per chunk
         of 50 consecutive occurrences) has median ≥ 0.4.
   4. Cluster surviving k-mers into independent groups via prefix/suffix
      overlap (union-find on shared (k−1)-mer).
   5. Refine period P̂ via joint least-squares slope fit through all
      occurrences in all surviving k-mers.
   6. For each group, for each occurrence x_i: compute deviation from
      the nearest anchor (circular), bin by monomer index j = ⌊x_i / P̂⌋.
   7. Aggregate: per monomer j, take median deviation across all groups
      that have signal there.  → d[j].
   8. Unwrap d[j] (wrap-around guard).
   9. Detrend: subtract any residual linear slope (a wrong P̂ produces a
      ramp, real same-sign indels also do — we report the slope as a
      diagnostic but operate on the residual).
  10. Compute adaptive noise floor σ̂ as the MAX of three sources:
        σ_iid    = MAD of (d − median_filter(d, 5)) × 1.4826 × 3
        σ_rw     = MAD of (Δd) × 1.4826 × 3 × √min_seg
        σ_P      = 0.05 × P̂
  11. Binary segmentation on detrended d[j] with reduction threshold
        σ̂² × min_seg, recurses top-down, min segment length = 10.
  12. Per-event filter: each candidate changepoint must clear
        |consensus_step|       ≥ σ̂
        n_supporting_groups    ≥ max(3, 0.5 × n_groups)
        same_sign_fraction     ≥ 0.7
  13. Aggregate metrics and write one row to the irregularity TSV.
```

### 2.3 Algorithm step-by-step

#### Step 1 — Period

The scan needs the structural period `P` of the array. We take it from
`kite.tsv::monomer_size` (the rank-1 period from kite). The prototype
loads this directly from disk; the Rust port can take it via
`kite::KiteResult::peaks[0].period`.

Records with no entry in `kite.tsv` → `flag = no_period`.

#### Step 2 — Frequent k-mers

`k = 6` (matches kite's k for register-marker consistency). Forward
strand only, ACGT only — any non-ACGT byte in a k-mer skips it. Top-100
by frequency per record. Below the top, k-mers are too rare to give
meaningful position counts.

The 4⁶ = 4096-element k-mer space is small relative to corpus
diversity, so 100 is generous (rarely is there genuine diversity beyond
that).

#### Step 3 — Register-lock filter (the key design decision)

For each k-mer with sorted positions `x[0..n−1]`:

**3a. Greedy multi-anchor identification.** Phase histogram with 1-bp
bins over `x[i] mod P`. Repeatedly pick the densest 5-bp circular
window, record its center as an anchor, zero it out from the histogram.
Stop after `MAX_ANCHORS = 5` anchors or once cumulative coverage of
removed mass reaches `MIN_ANCHOR_COVERAGE = 0.4`.

**Why multi-anchor?** A k-mer can legitimately occur at several phases
within a single monomer (e.g. a 6-mer that appears at positions 50 and
150 within a 200-bp monomer has two anchor phases). Single-anchor
filtering rejected all such k-mers in long-period arrays. Multi-anchor
keeps them.

**3b. Windowed coverage filter.** Compute anchor coverage in
non-overlapping chunks of 50 consecutive occurrences each. Keep the
k-mer if the **median** chunk-coverage is ≥ `MIN_ANCHOR_COVERAGE`.

**Why windowed?** A long array with natural ±3 bp per-monomer drift
will spread the **global** phase across all of [0, P) over thousands of
monomers (random-walk std `≈ 3 × √N`), even though each local chunk of
~50 consecutive monomers is still sharply anchored. Global-coverage
filtering wrongly rejected such arrays (e.g. TRC_10:chr4_…, a 1 Mb
real `tr_with_subrepeat`). Windowed-coverage filtering correctly keeps
them — the *local* register is intact even if the *global* phase has
walked. Real broken-register arrays (TRC_16-style spread indels) fail
locally too and are correctly rejected.

This is the most important design decision in the prototype.

#### Step 4 — Effective-support groups

K-mers offset by 1 bp in the genome share a (k−1)-bp prefix-suffix
overlap and are essentially the same phase marker viewed through two
adjacent windows — not independent observations. Union-find:

- Each k-mer contributes its `kmer[:k−1]` (prefix) and `kmer[1:]`
  (suffix) as keys.
- Any two k-mers sharing a key are unioned.
- Connected components = independent k-mer groups.

This produces the `n_kmer_groups` reported per record, which is what
the per-event support filter operates on (not raw k-mer count).

We need `n_kmer_groups ≥ 3` to produce a reliable signal; otherwise
`flag = no_register_lock`.

#### Step 5 — Refined period P̂

Joint robust slope fit. Assign each occurrence an initial integer copy
index `j = round((x − φ_initial) / P)` using kite's `P` and the
candidate k-mer's initial densest-window phase. Stack all
`(j, x − φ_k_initial)` pairs across all surviving k-mers and fit
`y = j · P̂` by least squares through the origin.

Result: refined estimate `P̂` that captures any per-monomer length bias
the k-mer ladder reveals. We use **`P̂` everywhere downstream**, not
the kite-supplied `P`.

If `|P̂ − P|/P > 0.05`, we still proceed but flag `period_unstable`.

#### Step 6 — Per-occurrence deviation

For each k-mer group, pool all occurrences from all member k-mers and
recompute anchors using `P̂`. For each occurrence `x`:

```
j      = floor(x / P̂)              # monomer index
d      = nearest-anchor distance(x mod P̂, anchors)  # in [−P̂/2, P̂/2]
```

`d` is positive when x is "later" in the monomer than the nearest
anchor, negative when earlier.

#### Step 7 — Per-monomer aggregate

For each monomer index `j` (over the range observed):

```
d[j]      = median of group medians at j
support[j] = number of groups with any occurrence at j
dropout[j] = n_kmer_groups − support[j]
```

#### Step 8 — Unwrap

Standard phase-unwrap: walk `j`, when `|d[j] − d[j−1]| > P̂/2`, add or
subtract P̂ to keep continuity. Required because cumulative drift can
wrap modulo `P̂` over long arrays.

#### Step 8b — Detrend

Linear regression `d[j] = a + b·j + ε`. Subtract the line. Report `b`
as `period_drift_bp_per_copy` — diagnostic for whether `P̂` is correct
(near 0) vs the array has same-direction cumulative drift (significantly
nonzero).

Operate on the residual `d_detrended[j]` for steps 9–12.

#### Step 9 — Adaptive noise floor σ̂

Three sources, take the max:

| component | what it captures | formula |
|---|---|---|
| σ_iid | per-point smoothness — substitution/jitter noise | `1.4826 · median(|d − median_filter(d, 5)|) · 3` |
| σ_rw  | random-walk variance over a min_seg-long segment — natural per-monomer drift accumulated | `1.4826 · median(|Δd|) · 3 · √min_seg` |
| σ_P   | fraction-of-period absolute floor — anything below 5 % of P is unlikely to be a meaningful structural event | `0.05 · P̂` |

The `σ_rw` term is the most important addition during calibration. Without
it, natural drift (a random walk in `d[j]`) produces extreme-value
"spike" pairs that binary segmentation interprets as up-then-down
indels. Including `σ_rw` makes the changepoint threshold scale with
the array's actual per-step drift variance.

#### Step 10 — Binary segmentation

Top-down recursive split on the detrended signal. For each segment, find
the split point that maximizes SSR reduction (sum of squared residuals
before vs after split). If the reduction exceeds `σ̂² · min_seg`, accept
the split and recurse on both halves. Minimum segment length 10 means
any event must be flanked by ≥ 10 copies of stable register on each
side.

Implementation uses prefix sums over `d` and `d²` so each `best_split`
call is O(N) — without that, binary segmentation runs O(N²) per
record and hangs on Mb-scale arrays.

A hard cap `MAX_CHANGEPOINTS = 500` is enforced via the recursion guard
as a safety against pathological signals.

#### Step 11 — Per-event same-sign support filter

Every candidate changepoint `j*` is independently evaluated:

```
For each group g with ≥ min_seg points on each side of j*:
  step_g = median(d_g in [j*, j*+min_seg]) - median(d_g in [j*-min_seg, j*])
consensus_step = median(step_g over groups)
n_supporting = #{g : sign(step_g) == sign(consensus_step) AND |step_g| ≥ σ̂}
same_sign_fraction = n_supporting / n_groups_present
```

Keep the event iff:

```
|consensus_step|        ≥ σ̂
n_supporting            ≥ max(3, 0.5 × n_kmer_groups)
same_sign_fraction      ≥ 0.7
```

This is the "independent observers agree" check. An accidental
changepoint from a few noisy groups gets rejected.

#### Step 12 — Aggregate metrics

```
step_count                 = count of kept events
phase_shift_burden_pct     = (Σ |step_size|) / array_length × 100
net_drift_bp_per_kb        = |Σ step_size| / array_length × 1000
jitter_bp                  = 1.4826 × median(|d_detrended − median_filter(d_detrended, 5)|)
dropout_rate               = mean over j of dropout[j] / n_kmer_groups
period_drift_bp_per_copy   = the linear slope b from step 8b
```

### 2.4 Output schema

`<prefix>.irregularity.tsv` — 13 columns, one row per record:

| # | column | type | description |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `period_P` | float | kite-supplied period (rank-1 monomer_size) |
| 4 | `refined_P` | float \| NA | step-5 refined period |
| 5 | `period_drift_bp_per_copy` | float \| NA | residual slope after refined-P fit |
| 6 | `n_kmer_groups` | int | independent k-mer groups after overlap clustering |
| 7 | `jitter_bp` | float \| NA | robust spread of detrended d[j] |
| 8 | `step_count` | int \| NA | number of detected indel-like events |
| 9 | `phase_shift_burden_pct` | float \| NA | accumulated structural shift per array-length, % |
| 10 | `net_drift_bp_per_kb` | float \| NA | net cumulative drift in bp/kb |
| 11 | `dropout_rate` | float \| NA | mean fraction of expected groups absent per monomer |
| 12 | `flag` | str | one of `ok`, `no_period`, `too_short`, `no_register_lock`, `period_unstable`, `highly_irregular`, `too_long` |
| 13 | `notes` | str | free-text diagnostic |

#### Flag semantics

| flag | metrics filled? | meaning |
|---|---|---|
| `ok` | all | scan completed, results are reliable |
| `no_period` | none | record has no kite period (missing from `kite.tsv`) |
| `too_short` | none | array_length < 10 × P (post-aggregation < 2 × min_seg copies) |
| `no_register_lock` | none | < 3 register-locked k-mer groups at kite P |
| `period_unstable` | filled, but trust diminished | refined P̂ disagrees with kite P by > 5% |
| `highly_irregular` | none | edge fallback: k-mers locked at intervals near P but anchor coverage too low even per-window. Rarely hit on the current IPIP corpus. |
| `too_long` | none | post-aggregation > 10,000 copy indices — safety cap |

### 2.5 Tunables and their defaults

Everything is a top-of-file constant in `irregularity.py`:

| constant | default | meaning |
|---|---|---|
| `K` | 6 | k-mer length (must match kite) |
| `TOP_KMERS` | 100 | top-N most frequent k-mers per record |
| `MAX_ANCHORS` | 5 | max anchor windows per k-mer |
| `MIN_ANCHOR_COVERAGE` | 0.4 | median windowed coverage required |
| `PHASE_WIN_BP` | 5 | width of each anchor window in bp |
| `MIN_KMER_GROUPS` | 3 | min independent groups to enable analysis |
| `MIN_COPIES_FOR_SCAN` | 10 | min copies for the early-gate `too_short` check |
| `MIN_SEG` | 10 | min copies on each side of a changepoint |
| `NOISE_FLOOR_MIN_BP` | 1.0 | minimum σ̂ in bp |
| `NOISE_FLOOR_K` | 3.0 | σ multiplier on MAD |
| `STEP_MIN_FRAC_OF_P` | 0.05 | period-relative absolute floor on step size |
| `SAME_SIGN_FRAC_MIN` | 0.7 | per-event same-sign requirement |
| `MIN_SUPPORT_FRAC` | 0.5 | per-event support fraction |
| `MAX_COPY_INDICES` | 10000 | hard cap to skip pathological arrays |
| `MAX_CHANGEPOINTS` | 500 | recursion guard in binary segmentation |
| `TOL_MODAL` | 0.05 | tolerance for "interval ≈ P" check in the highly_irregular fallback |

### 2.6 Calibration history

The algorithm went through several iterations during development.
Documented here because the design choices are non-obvious and the
"why not the simpler version" matters for further work.

**v0 — single-anchor phase**: each k-mer has one phase mod P, computed
as the densest 5-bp window's center. Failed on multi-occurrence-per-
monomer k-mers (common in long-period clean arrays).

**v1 — multi-anchor + global coverage gate at 0.7**:
greedy multi-anchor, gate on global coverage ≥ 0.7. Failed on:

- TRC_16-style arrays (real spread indels): k-mers are register-locked
  locally but the global phase has wrapped → coverage 0.28 → rejected.
- TRC_104-style arrays (small natural drift, short array): coverage
  ~0.56 just below the gate → rejected.

**v2 — gate at 0.4 + step floor `0.05 × P` + `min_seg = 5`**: better,
but TRC_104 still flagged with several spurious events from random-walk
extrema. TRC_100 (clean simple_tr) got 14 false-positive events at
auto-P=18 (a sub-monomer scale, see next).

**v3 — drop auto-P fallback + `min_seg = 10` + step floor remains**:
removed the auto-P "if kite P doesn't work, try a different period"
fallback because it dragged clean arrays into wrong-scale sub-period
analysis. Records that don't work at kite P now report
`no_register_lock` rather than producing meaningless metrics at a
sub-monomer scale.

**v4 — random-walk-aware σ̂**: the breakthrough for natural-drift cases.
Added `σ_rw = 1.4826 · median(|Δd|) · 3 · √min_seg` to the combined
noise floor. This explicitly models the fact that under natural
per-monomer length drift, `d[j]` is a random walk whose mean over a
min_seg-length segment has std `σ_step · √min_seg` — and binary
segmentation will find "step changes" of that scale even in pure
random-walk noise. With σ_rw included, TRC_104 → 0 events.

**v5 — windowed anchor coverage (current)**: changed the filter from
global to windowed (median over 50-occurrence chunks). Rescued long
arrays with natural drift like TRC_10 (1 Mb, ±3 bp/monomer drift)
from being categorically rejected as `highly_irregular`. Combined
with v4 σ̂, drift-only arrays now correctly report `flag=ok` with
`step_count=0`. The `highly_irregular` flag is now rarely (if ever)
triggered on real corpora.

The full corpus distribution evolved through these iterations:

| version | ok | no_register_lock | too_short | too_long | highly_irregular |
|---|--:|--:|--:|--:|--:|
| v3 | 1085 (35.9%) | 1194 (39.5%) | 479 (15.8%) | 157 (5.2%) | 109 (3.6%) |
| v5 | **1698 (56.2%)** | 586 (19.4%) | 582 (19.3%) | 158 (5.2%) | 0 |

The `too_short` rate increased from v3 → v5 because the v3 `min_seg = 5`
was bumped to 10. Records with 10–19 copies are now excluded from
analysis (need ≥ 2 × min_seg = 20). This is a conservatism trade-off
that helped on TRC_104.

### 2.7 Known limitations

1. **Whole-P copy gain/loss is invisible.** An insertion or deletion of
   exactly `P` bp leaves all downstream phases unchanged modulo `P`,
   so this method cannot detect it. A separate "expected vs observed
   tile copies" metric would be needed.
2. **Long-P clean arrays may fail `no_register_lock`.** If the
   structural period reported by kite is e.g. 1500 bp but the top-100
   k-mers actually operate at a finer scale (e.g. they occur every
   ~30 bp because of an internal subrepeat), no k-mer's modal interval
   is ≈ 1500 and the filter rejects all of them. v3 had an auto-P
   fallback for this, but it caused calibration problems and was
   dropped — see calibration history above.
3. **High-multiplicity HORs (e.g. 50 monomers per tile, P~8000 bp at
   monomer-scale of 170 bp)** confuse the period choice — kite may pick
   the monomer (170) or the tile (8000) as rank-1. We follow kite's
   choice. If kite picks the monomer, irregularity is measured at the
   monomer scale (which is fine). If kite picks the tile, we likely
   fail register-lock at that long scale.
4. **`too_short` is conservative.** Arrays with 10–19 copies are
   excluded because `min_seg = 10` needs ≥ 20 total. A separate
   "very-short-array" branch with `min_seg = 3` could rescue these,
   but those few-copy arrays generally don't have enough information
   for meaningful changepoint detection anyway.
5. **`step_count` doesn't tell you about event size distribution.**
   Two arrays with `step_count = 10` may have completely different
   biological meaning if one's events are all ~5 bp and the other's
   are all ~50 bp. Use `phase_shift_burden_pct` and `net_drift_bp_per_kb`
   alongside `step_count`.

### 2.8 Performance notes

Bottlenecks identified during development:

- **`best_split` in binary segmentation** was originally O(N²) per call
  because it used Python list slicing for the means and SSRs. Replaced
  with O(N) prefix-sum lookups in `cum[]` and `cum2[]` — multi-Mb
  records went from "hangs" to <1 s.
- **`MAX_COPY_INDICES = 10000` cap** skips records where the
  post-aggregation signal would be longer than 10000 monomers. These
  are typically multi-Mb arrays where the auto-P (now disabled)
  pulled the period down to a sub-monomer scale; with auto-P gone,
  hits to this cap are rare.
- **`MAX_CHANGEPOINTS = 500` recursion guard** ensures the binary
  segmentation can't run away on a pure-noise signal.

Full corpus (3024 records, total ~1.2 GB FASTA) currently runs in
~15 minutes on a single Python process. Acceptable for prototype; will
be 1–2 orders of magnitude faster in Rust.

### 2.9 Files of interest

```
tools/rule_proto/irregularity.py                            # the prototype
tmp/IPIP200579_2026-04-14/result_v011.irregularity.tsv      # latest corpus output
tmp/IPIP200579_2026-04-14/highly_irregular.fasta            # records flagged highly_irregular (currently empty for v5)
```

To run:

```bash
python3 tools/rule_proto/irregularity.py \
    --fasta <input.fasta> \
    --kite <prefix>.kite.tsv \
    -o <out_prefix>
# Optional: also writes <out_prefix>.irregularity.copies.tsv
python3 tools/rule_proto/irregularity.py ... --dump-copies
```

The TSV is written incrementally (every 25 records), so `wc -l` /
`head` / `grep` on the partial file during the run is fine.

### 2.10 Rust port plan

When ready to port:

1. **Module structure** mirroring other stages:
   ```
   src/irregularity/
       mod.rs      # public entry point + Config + run_subcommand
       scan.rs     # main algorithm: anchors, deviation, changepoint
       io.rs       # TSV writer
   ```
2. **CLI**: add `RuleIrregularityArgs` in `src/cli.rs`, `Command::Irregularity`
   variant, dispatch in `src/main.rs::run_irregularity`. Inputs: `--fasta`,
   `--kite` (kite TSV with period column), `-o`. Optional `--dump-copies`.
3. **Reuse**:
   - `crate::sequence::ArrayRecord` for FASTA loading
   - K-mer counting can borrow conceptually from `src/kite.rs`'s
     `compute_neighbor_profile` (but we need position lists, not just
     the histogram — likely a parallel function)
4. **Tests**:
   - Unit tests in `scan.rs::tests`: anchor finding on synthetic
     positions; circular unwrap correctness; binary segmentation min-
     step-size threshold; per-event same-sign filter
   - Integration test against ≥3 known cases (TRC_10 ok, TRC_104 ok,
     TRC_18:chr1_464573806_… ok with non-zero burden) like the existing
     `tests/analyze_ci_corpus.rs`
5. **Performance**: use prefix sums for binary segmentation from day
   one; precompute k-mer position lists once per record;
   parallelize over records via `rayon` like the rest of the pipeline.
6. **Integration into `analyze`**: once stable, add as the 6th per-stage
   TSV. Optionally feed `step_count` / `phase_shift_burden_pct` into
   the summary cascade as a derived column or new `combined_class`
   category (e.g. `tr_with_indels`) if the user wants it surfaced
   there.

---

## Part 3 — Irregularity v2 (distance-residual)

A second prototype at `tools/rule_proto/irregularity_v2.py` implements
the distance-residual approach (Petr's Approach 1) with explicit
dropout/indel separation (Approach 6). On the `ground_truth3`
simulated benchmark it cleanly outperforms v1.

### 3.1 Why a second prototype

v1 (phase-mod-P + multi-anchor + step-detection on monomer-indexed
signal) has two characteristics that turned out to limit it on the
simulated indel-event benchmark:

1. **Phase-mod-P signal accumulates** along the array. An indel of +δ
   shifts every downstream copy's deviation by +δ, so multiple discrete
   indel events appear as one big drift segment with merged step
   changes. **v1 saturates at ~1.5 detected events even when the
   simulator injected 20.**
2. **Sequence-overlap k-mer clustering** collapses all top-K k-mers
   from a short-period array (≤100 bp) into one effective group, since
   they share heavy 5-bp overlaps with each other. The
   `n_groups ≥ 3` gate then spuriously fails on subrepeat-style arrays.

v2's residual formulation (one local spike per indel) avoids the
saturation problem, and its **phase-bin clustering** avoids the
short-period collapse.

### 3.2 Core algorithm

```
For each FASTA record:
  1. Read kite-supplied period P.
  2. Extract top-100 6-mers + their positions (same as v1).
  3. Per-k-mer modal-distance d_k = mode of consecutive interval.
     Keep k-mer iff d_k is compatible with the structural period —
     i.e. d_k ≈ m·P for m ∈ {1,2,3} or d_k ≈ P/m for m ∈ {1..5}.
     Also require ≥60% of intervals near n·d_k for small integer n
     (handles missing-occurrence "2P/3P" gaps cleanly).
  4. Cluster k-mers by phase mod P into k-wide non-overlapping bins.
     K-mers in the same bin share the same monomer-internal position
     → not independent observations. K-mers in different bins ARE
     independent. (See §3.3 for why this beats sequence-overlap
     clustering.)
  5. Per k-mer, compute consecutive-pair residuals:
        n_skipped = round(d / d_k) - 1
        residual  = d - round(d / d_k) · d_k
     n_skipped > 0 signals a missing-occurrence event (drop-out).
     residual ≠ 0 signals coordinate displacement (indel).
  6. Aggregate by genomic bin (one bin = one P-wide monomer slot).
     For each bin, per group take MEDIAN residual; across groups
     compute consensus residual, support count, same-sign count.
  7. Per-array baseline jitter:
        baseline = max(1 bp, 1.4826 · MAD(all residuals) · 3)
  8. Per-bin event call: |consensus| ≥ max(baseline, 0.05·P)
                          AND support ≥ max(3, 0.5·n_groups)
                          AND same_sign_fraction ≥ 0.7
     Merge adjacent same-sign passing bins into one event.
  9. Aggregate metrics: n_events, burden_pct, max_shift_bp,
     drift_bp_per_kb, dropout_rate_per_pair.
```

The structural difference from v1: **each pair's residual is local**.
An indel of +δ between two specific consecutive occurrences shows up
as a single spike in one bin. The next pair (over the indel) has
residual ≈ 0. v1's phase signal would be shifted by +δ persistently
until the next indel; v2's signal stays at zero except at the spike.

### 3.3 Phase-bin clustering — the load-bearing fix

v1's k-mer-clustering criterion was "share ≥5-bp prefix or suffix"
(union-find). For long periods this gives ~30-50 independent groups
from top-100 k-mers, which is fine. For short periods (≤100 bp,
common in subrepeat-style arrays), all top-100 6-mers occupy a single
short sub-motif and share heavy overlap with each other → **single
union-find group**, fails the `n_groups ≥ 3` gate.

v2 instead groups k-mers by their **modal phase mod P** into
non-overlapping bins of width `k` bp (= one k-mer width). This gives
⌈P / k⌉ candidate bins; only the occupied ones are returned. K-mers at
the same phase position share the SAME monomer-internal location and
are correctly counted as a single observation; k-mers at different
phase positions are independent observations of the same indel.

Crucial implementation detail: **do not pool positions across k-mers
within a bin** before computing residuals — pooling creates a fake
modal distance (the within-bin spacing). Instead, compute residuals
per-k-mer with that k-mer's own `d_k`, then concatenate within the
group.

On the `ground_truth3` corpus this fix alone takes the
short-period stratum's analyse rate from 0% (v2 pre-fix) to
**100%** (v2 with phase-bin clustering) — every `sub_` record now
gets analyzed.

### 3.4 v2 output schema

`<prefix>.irregularity_v2.tsv` — 14 columns:

| # | column | meaning |
|---|---|---|
| 1 | `record_id` | FASTA record identifier |
| 2 | `array_length` | sequence length in bp |
| 3 | `period_P` | kite-supplied period |
| 4 | `n_kmer_groups` | independent phase-bin groups |
| 5 | `n_pairs_total` | sum of consecutive-pair observations |
| 6 | `baseline_jitter_bp` | per-array baseline residual scale |
| 7 | `indel_event_count` | discrete indel events detected |
| 8 | `indel_burden_pct` | `(Σ \|step\|) / array_length × 100` |
| 9 | `indel_max_shift_bp` | size of the largest detected event |
| 10 | `indel_drift_bp_per_kb` | net cumulative drift in bp/kb |
| 11 | `dropout_event_count` | sum of n_skipped across all pairs |
| 12 | `dropout_rate_per_pair` | dropouts / n_pairs_total |
| 13 | `flag` | `ok` / `no_period` / `too_short` / `no_register_lock` / `too_long` |
| 14 | `notes` | free-text diagnostic |

The split between `indel_event_count` (discrete localized shifts) and
`dropout_rate_per_pair` (missing-occurrence events that don't disturb
register) is **the explicit Approach-6 dropout/indel separation**, and
the corpus calibration below confirms the two metrics are orthogonal.

### 3.5 `ground_truth3` calibration corpus

Located at `ground_truth3/` — a focused 1400-case test set for
subrepeat + indel-irregularity calibration, generated by:

```bash
python3 ground_truth3/build_grid.py     # writes params.tsv
python3 ground_truth3/run_parallel.py \  # 6-worker driver (~10 min)
    -p ground_truth3/params.tsv \
    -o ground_truth3 \
    -s 20260601 -j 6
```

**12 strata × 100 cases** (300 for `sub`):

| family | strata | what it tests |
|---|---|---|
| Subrepeat (sub) | `sub` × 300 | Monomer 300-3000 bp carries 21-100 bp internal sub-motif. `hor_order=1`, zero indels. Should NEVER trigger irregularity events. Tests `tr_with_subrepeat` detection of the summary cascade. |
| Per-base indel rates (irreg) | `irreg00`, `irreg0001`, `irreg0005`, `irreg002`, `irreg005`, `irreg01` × 100 | Plain tandem repeat, monomer 100-1000 bp, with `indel_rate_inter` ∈ {0, 0.001, 0.005, 0.02, 0.05, 0.10}. Produces *random-walk-drift* noise (many tiny indels). Tests false-positive rate as indel noise increases. |
| Discrete indel events (event) | `event01`, `event03`, `event05`, `event10`, `event20` × 100 | Plain tandem repeat with N discrete localized indel events, each 5-30 bp insertion or deletion at random positions. Recorded in `events.tsv` (scope=`indel`). Tests recall of indel-event detection vs known true count. |

The simulator (`ground_truth/simulate_hor.py`) was extended in this
cycle with three new params: `indel_events`, `indel_event_size_min`,
`indel_event_size_max`. These trigger the new `apply_indel_events()`
function which acts on the assembled sequence AFTER monomer
construction, picking random positions, choosing insert/delete with
50/50 probability, and logging each event to `events.tsv` with
`source_idx=position`, `target_idx=signed_size`. Default `indel_events=0`
preserves backward compatibility.

### 3.6 v1 vs v2 head-to-head on `ground_truth3`

| stratum | n | true events | v1 mean | v1 ok% | v2 mean | v2 ok% | v2 recall |
|---|--:|--:|--:|--:|--:|--:|--:|
| `sub`         | 300 |  0.0 | 0.00 | 82.3% | 0.00 | **100%** | — |
| `irreg00`     | 100 |  0.0 | 0.00 | 93.0% | 0.00 | **100%** | — |
| `irreg0001`   | 100 |  0.0 | 0.00 | 90.0% | 0.00 | 100% | — |
| `irreg0005`   | 100 |  0.0 | 0.00 | 85.0% | 0.00 | 100% | — |
| `irreg002`    | 100 |  0.0 | 0.07 | 90.0% | 0.00 | 100% | — |
| `irreg005`    | 100 |  0.0 | 0.12 | 86.0% | 0.03 | 100% | — |
| `irreg01`     | 100 |  0.0 | 0.03 | 88.0% | 0.02 | 100% | — |
| **`event01`** | 100 |  1.0 | 0.23 | 91.0% | **0.43** | 100% | **43%** |
| **`event03`** | 100 |  3.0 | 0.61 | 89.0% | **1.22** | 100% | **41%** |
| **`event05`** | 100 |  5.0 | 1.04 | 89.0% | **1.88** | 100% | **38%** |
| **`event10`** | 100 | 10.0 | 1.45 | 84.0% | **3.86** | 100% | **39%** |
| **`event20`** | 100 | 20.0 | 1.47 | 90.0% | **7.23** | 100% | **36%** |

Three clean wins for v2:

1. **No false rejections** — every record gets analysable (`ok` for all
   strata). v1 rejects 7-18% of records as `no_register_lock` /
   `too_short` / `too_long`.
2. **No false positives on per-base indel noise** — v2 reports ~0
   events even at the highest per-base rates. v1 starts producing
   spurious 0.07-0.12 mean events as the rate climbs.
3. **Linear scaling on discrete events** — v2 recovers a stable
   **~36-43% of simulated events** with a near-constant slope. v1
   saturates at ~1.5 events regardless of true count (the saturation
   problem of the phase-mod-P signal).

The ~36-43% under-detection is the cost of the conservative
`0.05 × P` step floor: events of 5-15 bp are below threshold for P ≥
~200 bp. Per Petr's guidance, this trade-off (under-detection over
false-positives) is preferred.

The corpus also separately confirms the **Approach 6 dropout/indel
split** is working orthogonally: `dropout_rate_per_pair` cleanly
tracks the per-base indel rate (0.22 baseline → 0.69 at rate=0.10)
while staying at baseline for the discrete-event strata.

### 3.7 Known limitations of v2

1. **Under-detection of small events** at long periods. A 5-bp indel
   in a P=500 monomer is below the 0.05·P=25 bp step floor and gets
   filtered. Stable across event counts so the ranking is preserved.
2. **Discrete-event count vs simulated count is a slope, not a
   one-to-one match.** The ~0.36-0.43× factor is approximately linear;
   actual indel counts in real data should be inferred from
   `indel_event_count / 0.4` as a rough conversion.
3. **No event-position precision report.** `indel_max_shift_bp` is
   reported but the per-event genomic positions are not (would need
   the per-bin trace TSV).
4. **`MIN_NEAR_DK_FRAC = 0.6` is permissive.** A k-mer whose 40% of
   consecutive intervals are off-pattern still passes the filter.
   Tightening to 0.8 would reduce noise k-mers at the cost of fewer
   surviving k-mers per record.

### 3.8 Files of interest (v2)

```
tools/rule_proto/irregularity_v2.py             # the prototype
tools/rule_proto/irregularity.py                # v1, kept for comparison
ground_truth/simulate_hor.py                    # extended with apply_indel_events()
ground_truth3/build_grid.py                     # 1400-case grid generator
ground_truth3/run_parallel.py                   # multiprocessing driver
ground_truth3/sequences.fasta                   # 35 MB, regenerable
ground_truth3/truth.tsv                         # per-record params + identities
ground_truth3/events.tsv                        # per-event log (scope=indel)
ground_truth3/result.irregularity.tsv           # v1 output on the corpus
ground_truth3/result.irregularity_v2.tsv        # v2 output on the corpus
```

To run end-to-end on the calibration corpus:

```bash
# 1. regenerate inputs (if not present)
python3 ground_truth3/build_grid.py
python3 ground_truth3/run_parallel.py \
    -p ground_truth3/params.tsv -o ground_truth3 -s 20260601 -j 6

# 2. run kite
./target/release/kitehor kite-periodicity ground_truth3/sequences.fasta \
    -o ground_truth3/kite.tsv

# 3. run both irregularity prototypes (parallelisable)
python3 tools/rule_proto/irregularity.py \
    --fasta ground_truth3/sequences.fasta --kite ground_truth3/kite.tsv \
    -o ground_truth3/result
python3 tools/rule_proto/irregularity_v2.py \
    --fasta ground_truth3/sequences.fasta --kite ground_truth3/kite.tsv \
    -o ground_truth3/result
```

---

### Glossary

- **register-locked k-mer**: a k-mer whose positions within an array
  consistently land at the same phase (or one of a few phases) mod the
  structural period P. The marker we use to track phase shifts.
- **anchor (phase)**: a peak in the histogram of `positions mod P` for
  one k-mer. A k-mer can have multiple anchors when it occurs at
  multiple positions within each monomer.
- **windowed anchor coverage**: anchor coverage computed in non-
  overlapping chunks of 50 consecutive occurrences, then medianed
  across chunks. Robust to global phase wrap from natural drift.
- **copy index `j`**: integer monomer index along the array, `j = ⌊x / P̂⌋`.
- **deviation `d[j]`**: aggregate signal — median over k-mer groups of
  each group's nearest-anchor distance for occurrences at monomer `j`.
- **detrended `d[j]`**: `d[j]` with the residual linear slope subtracted.
  Real indels appear as step changes; the slope itself is reported
  separately as `period_drift_bp_per_copy`.
- **changepoint**: an index `j*` at which the mean of `d[j]` differs
  significantly between `[j*-min_seg, j*]` and `[j*, j*+min_seg]`.
- **σ_iid**: noise floor component from per-point smoothness
  (substitutions, k-mer drop-outs).
- **σ_rw**: noise floor component from random-walk variance (natural
  per-monomer drift accumulated over min_seg copies).
- **σ_P**: noise floor component from `0.05 × P̂` (period-relative
  absolute floor — anything below 5% of the period is not a meaningful
  structural event).
- **phase_shift_burden_pct**: `(Σ |step_size|) / array_length × 100`.
  Total *absolute* structural disruption — large even if individual
  events cancel out.
- **net_drift_bp_per_kb**: `|Σ step_size| / array_length × 1000`. Net
  cumulative drift — small when up- and down-shifts cancel.
- **dropout_rate**: fraction of expected k-mer groups absent on
  average per monomer. Primarily a substitution / divergence
  diagnostic — orthogonal to indel signal.
- **d_k (per-k-mer modal distance)**: in v2, the modal consecutive
  interval for one k-mer — typically `P` for register-locked k-mers,
  but may be `2P` (missed-occurrence dominant), `3P`, or `P/m` for
  multi-occurrence-per-monomer k-mers. Used as the reference for
  residual computation.
- **residual (v2)**: `d − round(d / d_k) · d_k` for a consecutive pair
  of occurrences. Zero for perfect-register or missing-occurrence
  gaps (2P → 0); nonzero for indel-shifted spacings.
- **phase-bin clustering (v2)**: grouping k-mers by `(modal phase mod P) // k`
  into non-overlapping k-bp wide bins. Each occupied bin = one
  independent observation of the array structure at that phase.
  Replaces the older sequence-overlap union-find which collapsed
  short-period arrays into a single group.

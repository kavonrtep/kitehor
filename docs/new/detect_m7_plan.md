# Detect M7 — per-segment recompute + same-width mixed detection

Status: **REVIEWED 2026-05-16**. All Q1–Q11 decisions resolved per
`docs/reviews/detect_m7_plan_review_2026-05-16.md`; scope tightened
to HOR-focused same-width mixed detection. Ready for M7.1 code.

## TL;DR

M7 closes the two architectural deferrals from
`docs/reviews/detect_implementation_review_completed_2026-05-16.md`
(findings #2 and #3), tracked as A16 in
`docs/new/detect_impl_plan.md` §0. After the 2026-05-16 review, M7
is scoped narrower than the original draft — it is **detector-side,
HOR-focused** same-width mixed detection. `simple_TR` stratification
is held out as a separate taxonomy decision (see §1.3 + Q5 below).

- **Internal per-block analysis** (NOT reported segments) —
  whole-array sequence is split into fixed analysis blocks for the
  consensus-identity computation only. `segments.tsv` rows continue
  to be emitted only for biologically meaningful events
  (phase-shift boundaries today; mixed sub-blocks added by M7).
  Clean single-family arrays keep `Properties.n_segments = 1`.
- **Same-width mixed detection** — adjacent satellite families that
  share base width and HOR multiplicity (e.g.,
  `mx_a200-08_b200-08_n050-050`) collapse to a single `HOR` call
  today. The kite → detect run on `ground_truth_v2/` showed
  `mixed` accuracy crash from 96 % to 18 % under realistic
  (non-oracle) period inputs. M7 adds segment-level **consensus
  identity** comparison gated by `DetectorConfig.stratification_*`
  thresholds (currently inert).

The two are strictly coupled — same-width mixed detection is a
function of per-segment consensuses — so they ship together.

## 1. Context and motivation

### 1.1 Why this isn't fixed yet

`src/detect/segment.rs::split()` (current ~50 LOC) splits the array
at `Properties.phase_shift_positions` and emits one `Segment` row
per inter-shift region. Every emitted segment inherits whole-array
fields. The MVP comment at the top of the file calls out the
deferral explicitly:

> For MVP M4, segments inherit the whole-array class because all
> CI fixtures with phase shifts share the same underlying repeat
> architecture on both sides of every shift. A future iteration
> can recompute per-segment widths once HMM-based segmentation
> lands (OQ6).

HMM-based segmentation is itself deferred to v2 (§0 Q5). M7 uses
a **threshold-based** segmentation that doesn't depend on HMM
and is enough to surface the consensus-identity signal we need
for same-width mixed.

### 1.2 What it unblocks

| Pain point                                | Fixed by                          |
|-------------------------------------------|-----------------------------------|
| `mx_*_a200-08_b200-08` calls (mixed→HOR)  | consensus-identity gate           |
| `stratification_same/diff_threshold` inert| consumed in `segment::split()`    |
| `segments.tsv` "all rows identical" lie   | per-segment wrap / R(k) / shift   |
| Inability to localise within-array failures| per-segment classification info  |
| Confidence on mixed/ambiguous still poor  | segment-level evidence in scoring |

### 1.3 What it does NOT fix

Listed explicitly so we don't sprawl:

- **Strand-aware inversion** (OQ3). Inversions still call `HOR`.
- **Nested HOR detection** (T09 fixture is `.deferred.yaml`).
- **HMM-based segmentation** (§0 Q5). M7 uses fixed-stride
  blocks + phase-shift boundaries.
- **Real `inter_monomer_identity`** (review #5). The relabel-only
  comment from A15 stays; the value remains R(1).
- **Per-segment classification.** Each reported segment still
  inherits the whole-array class for M7; only mixed-detection
  features (block consensus + identity) are computed per block.
- **`simple_TR` stratification → `mixed`** (review-blocking #2).
  `T15_stratification.yaml` is two same-period `SIMPLE_TR` blocks
  with different monomer templates and current expectations
  (`tests/detect_expectations.tsv`) keep it `simple_TR`. M7 does
  NOT change that taxonomy. The mixed override applies to
  `HOR` / `irregular_HOR` only. Reopening the simple_TR-as-mixed
  question requires updating `docs/new/taxonomy.md`, `T15`, and
  the expectations file together — a separate decision.
- **Reporting every fixed block as a `segments.tsv` row.** Review
  blocker #1: analysis blocks are internal; only meaningful
  reported segments (phase shifts, mixed sub-blocks) are written.

## 2. Open questions — resolved 2026-05-16

The original draft posed eleven open questions. Review
(`docs/reviews/detect_m7_plan_review_2026-05-16.md`) resolved all
of them. Each subsection now records **(decided)** answer + the
short rationale. Where the decided answer differs from the
original proposal, the difference is called out explicitly.

### Q1. Segmentation criterion — how do we split when there are no phase shifts?

**(decided)** Fixed analysis blocks with adaptive size:
```
block_rows = max(block_size_rows_min, ceil(n_rows / max_segments_per_array))
```
Phase-shift positions are added as extra split points. Blocks
with fewer than `min_segment_rows` informative rows are merged
into a neighbour or skipped from the identity computation. **The
blocks are an internal `AnalysisBlock` only — they do NOT
correspond 1:1 to `segments.tsv` rows.** Reported segments stay
defined as today (phase-shift boundaries), with mixed-call
sub-blocks added by M7 only when the override fires.

For HOR / irregular_HOR arrays, block boundaries are **aligned to
full HOR units** (`hor_length_bp` boundaries) so that partial
units at block edges don't drive a `mixed` call by themselves.

Difference from original draft: the original proposed reporting
every fixed block in `segments.tsv` and using a hardcoded
`block_size_rows_min`. Adaptive sizing + analysis-vs-reported
separation come from the review.

### Q2. Consensus computation per segment

**(decided)** Majority-vote consensus, with the comparison width
chosen by the array's class:

- **`simple_TR`**: `base_width_bp` (only width available).
- **`HOR` / `irregular_HOR`**: prefer `hor_length_bp` consensus
  when enough complete HOR units fit in each block (default ≥ 3
  complete units per block); fall back to `base_width_bp` when
  unit-level consensus is unreliable. **Reason:** at base_width
  alone, two HOR families with the same `(base_width, k)` may
  share many monomer-slot consensuses by chance — the unit-level
  consensus is where the structural difference shows up.

Difference from original draft: the original used `base_width_bp`
unconditionally. Comparing at `hor_length_bp` for HORs is the
review's recommendation and is necessary for the
`mx_a200-08_b200-08` style cases we're trying to fix.

### Q3. Identity metric for comparing segment consensuses

**(decided)** Hamming identity, **ignoring positions where
either consensus has `N`**. Two values per pair:

- `identity`: matches / (matches + mismatches) over non-N positions
- `coverage`: (matches + mismatches) / length

If `coverage < 0.70` (configurable as `min_identity_coverage`),
return `None` for that pair — the comparison is uninformative.

Difference from original draft: the original used raw Hamming.
Without an N guard, a block with many indels (collapsed to N
columns by the wrap) would look perfectly identical to anything,
producing false negatives on mixed; without the coverage floor,
a block with one informative column would dominate the decision.

### Q4. Two thresholds, two decisions

**(decided)** Defaults unchanged:
- `stratification_same_threshold` = 0.90
- `stratification_diff_threshold` = 0.80

Three-way per pair:
- identity ≥ 0.90 → same family
- identity ≤ 0.80 → different family
- 0.80 < identity < 0.90 → borderline (diagnostic only,
  not a mixed trigger)

**Pair aggregation:** evaluate **all valid pairwise comparisons**
across analysis blocks (not only "to seg 1"). If *any* valid pair
has identity ≤ `stratification_diff_threshold` → mixed override
fires. Pairs with `coverage < min_identity_coverage` are excluded
from the test.

Difference from original draft: the original mixed pseudocode
compared only to block 1; review caught the mismatch with the
prose ("any pair"). All-pairs is the right rule for arrays with
≥ 3 blocks where the first block isn't the divergent one.

### Q5. Class-level interaction order — when does mixed fire vs irregular_HOR?

**(decided)** M7 order:
1. classify (HOR / simple_TR / mixed / ambiguous)
2. **mixed override via consensus identity** — **HOR /
   irregular_HOR only** — *(NEW)*
3. irregular_HOR demotion
4. wobble-dominance guard

**`simple_TR` is NOT eligible for the mixed override in M7.**
Reason: `T15_stratification` is the explicit fixture for
two-block same-period simple_TR with different monomers, and
the current expectation in `tests/detect_expectations.tsv` is
`simple_TR`. Enabling the override on simple_TR would silently
change the taxonomy. If we later decide simple_TR
stratification should become `mixed`, that's a separate joint
change to `taxonomy.md` + `T15` + the expectations file.

Difference from original draft: original proposed simple_TR
eligibility as "default yes"; review correctly rejected that as
a silent taxonomy change.

### Q6. Per-segment class assignment

**(decided)** Light option. Every reported segment inherits the
final whole-array class. Heavy per-segment `classify::decide_array`
deferred to M8+. Note in the Segment row docstring that
`Segment.class` is **not** an independent per-segment
classification result; it's the array's class applied to this
sub-range.

Difference from original draft: alignment with review — same
decision, just made explicit in docstring/header.

### Q7. Output schema impact

**(decided)** Three schema changes:

1. **New `segments.tsv` columns:**
   - `consensus_identity_to_reference` — Hamming identity to the
     reference block (medoid block, or block 0 if all-pairs
     median is ambiguous). Renamed from the original
     `consensus_identity_to_seg1` per review.
   - `consensus_identity_coverage` — fraction of non-N positions
     used (review #79: low-coverage identity must be
     distinguishable from a high-confidence comparison).
2. **`consensus.fa`:** emit per-segment consensus records **only
   when `Properties.class == Mixed`**. Naming:
   `<array_id>_seg<N>_monomer` (and `_hor_unit` if a unit-level
   consensus was the basis for the decision). Single-family
   arrays continue to emit just the whole-array monomer + optional
   hor_unit; no segment proliferation.
3. **`diagnostics.json` `schema_version`** bumps 1 → 2. The new
   per-array block exposes `analysis_blocks` (internal) + the new
   identity columns. Schema-drift tests bumped.

Difference from original draft: column name (review preference),
explicit coverage column (review must-have), explicit
"reported segments are not analysis blocks" boundary.

### Q8. Performance budget

**(decided)** New config knobs (added to `DetectorConfig`):

```rust
pub max_segments_per_array: usize = 32,
pub min_segment_rows: usize        = 20,
pub min_identity_coverage: f64     = 0.70,
pub min_complete_units_per_block: usize = 3,  // HOR-unit consensus floor
```

(Currently inert: `stratification_same_threshold`,
`stratification_diff_threshold`.)

Implementation rules:
- **Skip per-segment R(k)** — segments inherit class (Q6 light),
  so only wrap + column IC + consensus are needed per block →
  O(N) total per array.
- **No nested rayon** — array-level batch is already parallel via
  rayon at `detect-batch`. Don't add segment-level parallelism
  until profiling demands it.
- **Adaptive block size** — see Q1.

With these, the per-array overhead should be < 100 ms even on
the largest v2 arrays.

Difference from original draft: explicit config fields (review
flagged the implicit constants); no nested rayon (review).

### Q9. Backward compatibility with the M6 baseline (94.4 %)

**(decided)** Acceptance gate tightened per review:

- `mixed` on **kite-derived periods** ≥ 70 % (current 18 %; the
  78 pp regression is the load-bearing motivation).
- `mixed` on **oracle periods** ≥ 94 % (i.e., no more than a 2 pp
  drop from the current 96 %). This is the new constraint the
  review added — we must not lose existing wins.
- No other category drops by > 2 pp from the M6 baseline.
- Core CI fixtures (T01, T05, T06, T07, T10, T13, T17, T18) still
  pass exactly.
- Explicit false-mixed guards (Tests §):
  - clean HOR with high within-family divergence (`hor_clean d40`)
    must not become mixed
  - `hor_wobble` and `hor_shift` HORs must not become mixed
  - inversion (`hor_event_inversion`) — see Risks §
  - `T15_stratification` locked to `simple_TR`

If the thresholds can't meet these gates, we re-tune
`stratification_*` before merging.

Difference from original draft: oracle-side mixed floor and
fixture-locked false-mixed list both come from the review.

### Q10. Same-width mixed via kite-derived periods

**(decided)** Ship M7 with **no kite-side changes**. The new
mixed test runs on whole-array sequence at the chosen comparison
width — independent of the period candidate list — so the
regression should fix from the detector side alone. Re-run
kite → detect after M7 lands, write an updated report, and only
then decide whether `src/emit_periods.rs` needs further tuning.

### Q11. Diagnostics JSON schema version bump

**(decided)** Bump `diagnostics.json::schema_version` 1 → 2.
Schema-drift tests bumped accordingly. Downstream report
readers (`tools/detect_eval/eval.py`, `tools/detect_eval/report.py`)
must keep working — they read array_id + class + a small set of
properties columns, all stable.

## 3. Implementation plan

### 3.1 Modules touched

| File | Change |
|---|---|
| `src/detect/analysis_blocks.rs` *(new)* | Internal `AnalysisBlock { start_row, end_row, consensus, n_complete_units }`; builder respects `max_segments_per_array`, `min_segment_rows`, HOR-unit alignment. |
| `src/detect/segment.rs` | Keep `split(&Properties)` for phase-shift segments (today's behaviour). M7 adds a separate path that emits segment rows only when `Properties.class == Mixed`, one row per divergent analysis block. |
| `src/detect/types.rs` | Add two `Segment` fields: `consensus_identity_to_reference: Option<f64>` and `consensus_identity_coverage: Option<f64>`. Bump `SEGMENTS_HEADER`. |
| `src/detect/mod.rs::run_array_m4` | Insert step "mixed-override via analysis blocks" between classify and irregular demotion. Apply only to `Class::HOR` / `Class::IrregularHOR`. |
| `src/detect/classify.rs` | Expose helper `mixed_override_via_blocks()` so the override is callable from `mod.rs`. |
| `src/detect/consensus.rs` | Add `consensus_on_slice(seq, width, start_row, end_row)`. Reuse the existing majority-vote logic. |
| `src/detect/config.rs` | Add `max_segments_per_array=32`, `min_segment_rows=20`, `min_identity_coverage=0.70`, `min_complete_units_per_block=3`. Document their use. |
| `src/detect/io.rs` | Emit per-segment consensus to `consensus.fa` ONLY when class=mixed. Bump diagnostics `schema_version` to 2. |
| `tests/detect_*` | New fixtures + assertions (see §6). |
| `docs/new/detect_impl_plan.md` | New A19 amendment; mark M7 in §10. |
| `docs/reports/` | Post-M7 rerun report. |

### 3.2 Sub-milestones (PR-sized, after review)

**M7.1 — Analysis blocks + per-block consensus + identity stats (no class change).**
- New module `src/detect/analysis_blocks.rs` with the builder
  (adaptive sizing, HOR-unit alignment, `min_segment_rows` skip).
- New `consensus::consensus_on_slice` helper.
- Per-block identity computation with N-skip + coverage tracking.
- Wire into `run_array_m4` to compute the stats but **not** act
  on them.
- Tests:
  - block builder respects `max_segments_per_array` on huge arrays
  - HOR-unit alignment drops partial units at edges
  - Hamming-with-N-skip behaviour
  - low-coverage pair returns `None`
- **Acceptance**: existing 316 tests still pass; new tests pass;
  no benchmark regression (class behaviour unchanged).

**M7.2 — Mixed override for HOR / irregular_HOR.**
- Add `consensus_identity_to_reference` and
  `consensus_identity_coverage` to `Segment`; bump
  `SEGMENTS_HEADER`.
- In `run_array_m4`, after classify and before irregular demotion:
  if `Class::HOR | Class::IrregularHOR` AND any valid pairwise
  identity ≤ `stratification_diff_threshold` → rewrite class to
  `Mixed` with reason citing the divergent pair.
- Emit one `Segment` row per analysis block ONLY when the
  override fires (mixed path); keep clean arrays at
  `n_segments = 1`.
- Tests:
  - positive: 2-block same-`(base_width, k)` mixed → `mixed`
  - negative-clean-hor-d40: high within-family divergence stays
    `HOR`
  - negative-wobble: stays `HOR`
  - negative-phase-shift: stays `HOR`
  - negative-T15: stays `simple_TR` (override doesn't fire on
    simple_TR)
  - low-coverage pair: stays HOR
- **Acceptance**: `mixed` on v2 corpus ≥ 70 % (kite-derived) AND
  ≥ 94 % (oracle); no other category drops > 2 pp; M4 fixtures
  pass exactly.

**M7.3 — Schema + diagnostics + segment-consensus FASTA.**
- `consensus.fa`: emit per-segment consensus records ONLY when
  class=mixed, naming `<array_id>_seg<N>_monomer`.
- Bump diagnostics `schema_version` 1 → 2; update schema-drift
  tests.
- Update CLAUDE.md + README to mention the new `segments.tsv`
  columns + the mixed-only consensus rows.
- **Acceptance**: schema-drift tests bumped, `tools/detect_eval/`
  still parses every output (verified end-to-end).

**M7.4 — Docs + rerun reports + impl-plan amendment.**
- Add A19 to `detect_impl_plan.md` §0; mark M7 done in §10.
- Re-run kite → detect on v2; write
  `docs/reports/kite_to_detect_v2_<post-m7>.md` (or amend the
  existing report).
- Re-render the dashboard (`tools/detect_eval/report.py`).

### 3.3 New `Segment` row

```rust
pub struct Segment {
    pub array_id: String,
    pub segment_id: usize,
    pub start_bp: usize,
    pub end_bp: usize,
    pub class: Class,            // inherited from final array class
    pub base_width_bp: Option<usize>,
    pub hor_k: Option<usize>,
    pub column_conservation: Option<f64>,                  // per-segment IC
    pub phase_separation: Option<f64>,                     // per-segment phase_sep
    pub wobble_amplitude_bp: Option<f64>,                  // per-segment wobble
    pub irregularity_score: Option<f64>,                   // existing
    pub consensus_identity_to_reference: Option<f64>,      // NEW (M7.2)
    pub consensus_identity_coverage: Option<f64>,          // NEW (M7.2)
}
```

Only the last two fields are genuinely new. The reference block
is the medoid of the analysis-block set (block with highest sum
of pairwise identities to all others); ties broken by smallest
segment_id.

### 3.4 Pseudocode for the mixed override

```rust
fn mixed_override_via_blocks(
    seq: &[u8],
    blocks: &[AnalysisBlock],
    comparison_width: usize,   // hor_length_bp if available, else base_width_bp
    cfg: &DetectorConfig,
) -> Option<MixedDecision> {
    if blocks.len() < 2 {
        return None;
    }
    // Per-block consensus at the comparison width.
    let consensuses: Vec<Vec<u8>> = blocks
        .iter()
        .filter_map(|b| consensus::consensus_on_slice(
            seq, comparison_width, b.start_row, b.end_row,
        ))
        .collect();
    if consensuses.len() < 2 {
        return None;
    }
    // All-pairs identity with N-skip + coverage gate.
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..consensuses.len() {
        for j in (i + 1)..consensuses.len() {
            if let Some((ident, cov)) = hamming_identity_n_skip(
                &consensuses[i], &consensuses[j],
            ) {
                if cov >= cfg.min_identity_coverage {
                    pairs.push((i, j, ident));
                }
            }
        }
    }
    if pairs.is_empty() {
        return None;  // insufficient coverage anywhere
    }
    let min_pair = pairs
        .iter()
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        .copied()
        .unwrap();
    if min_pair.2 <= cfg.stratification_diff_threshold {
        return Some(MixedDecision {
            reason: format!(
                "mixed — block-consensus identity {:.3} ≤ diff_threshold {:.3} \
                 (blocks {} vs {})",
                min_pair.2, cfg.stratification_diff_threshold,
                min_pair.0, min_pair.1,
            ),
            // Reference block = medoid (highest sum of pairwise identities
            // to others); ties broken by smallest index.
            reference_block: pick_medoid(&pairs, consensuses.len()),
            divergent_blocks: pairs.iter()
                .filter(|p| p.2 <= cfg.stratification_diff_threshold)
                .flat_map(|p| [p.0, p.1])
                .collect::<HashSet<usize>>(),
        });
    }
    None
}
```

## 4. Risks

| Risk | Mitigation |
|---|---|
| False mixed on `hor_clean d40+` / `hor_wobble` | `diff_threshold` defaulted to 0.80; HOR-unit comparison + N-skip + coverage gate (Q3) keep noise off the test. Locked-in negative tests (§6) catch regressions. |
| Per-segment recompute slow on large arrays | Q8 mitigations: cap `max_segments_per_array=32`, skip per-block R(k), no nested rayon. |
| Inversion segments look mixed (RC consensus ≠ forward consensus) | Inversions defer strand-aware to v2 (OQ3). Two options: (a) explicit "skip mixed override when array is hor_event_inversion" path (needs a category-aware test — not available at detect time); (b) accept the false-mixed count as the cost of deferring strand awareness, document in `reason` field. **Default: (b)** — let benchmark numbers tell us if we need (a). |
| `simple_TR` stratification taxonomy change creeps in | Hard-coded class gate in `mod.rs`: mixed override fires only for `Class::HOR` and `Class::IrregularHOR`. `T15_stratification` test asserts `simple_TR`. |
| Borderline pairs (0.80 < id < 0.90) that should be mixed | Q4: borderline is diagnostic-only, not a trigger. If benchmark shows this is the common confusion, M8 candidate. |
| Schema bump in diagnostics.json breaks downstream eval | `tools/detect_eval/eval.py` and `report.py` read `array_id` + `class` + a fixed set of `Properties` columns; schema_version bump should be transparent. Verify via M7.3 acceptance. |
| Reference block instability (medoid changes when one block is borderline) | Tie-break by smallest segment_id; document in Segment struct comment. |

## 5. Acceptance criteria (M7 done when)

1. **`mixed` category on `ground_truth_v2/`:**
   - oracle periods: **≥ 94 %** (no more than 2 pp drop from
     today's 96 %)
   - kite-emitted periods: **≥ 70 %** (current 18 %)
2. No other category drops by > 2 pp from the M6 baseline (across
   both oracle and kite-derived runs).
3. Locked-in false-mixed negatives (each its own test):
   - clean HOR with `d40` divergence stays `HOR`
   - `hor_wobble` cases stay `HOR`
   - `hor_shift` cases stay `HOR`
   - `T15_stratification` stays `simple_TR`
   - `hor_event_inversion` — documented expected behaviour (see
     Risks); test asserts either current behaviour or the
     fixed-up behaviour explicitly
4. Positive case: new fixture for 2-block same-`(base_width, k)`
   HOR with distinct monomer templates → `mixed`.
5. Low-coverage pair test: blocks dominated by `N` columns must
   NOT trigger `mixed`.
6. Core CI fixtures (T01, T05, T06, T07, T10, T13, T17, T18) still
   pass exactly via `cargo test --release --test detect_m4`.
7. `cargo test --release` green; new unit tests for the analysis-
   block builder, Hamming-with-N-skip, and the mixed override
   pseudocode.
8. `tests/detect_kite_emit.rs` still passes; no schema-related
   regression.
9. `docs/new/detect_impl_plan.md` updated with A19 amendment and
   M7 marked done in §10.
10. `docs/reports/kite_to_detect_v2_<post-m7-date>.md` written
    with the new numbers and an updated dashboard.

## 6. Tests to add (locked, per review)

Each lives under `tests/detect_*.rs`:

| Test | Fixture | Expected |
|---|---|---|
| Positive: same-`(base_width, k)` 2-block HOR | new YAML (`tests/synth_configs/T20_same_width_mixed.yaml`) — two `HOR` blocks, identical base_width + k, different monomer templates | `Class::Mixed`, identity ≤ 0.80 in reason |
| Negative: clean HOR `d40` | exists (`hor_clean d40` cases) | stays `Class::HOR` |
| Negative: wobble HOR | exists (T03 / `hor_wobble`) | stays `Class::HOR` |
| Negative: phase-shift HOR | exists (T10 / `hor_shift`) | stays `Class::HOR` |
| Negative: `T15_stratification` | exists | stays `Class::SimpleTR` — override doesn't apply to simple_TR |
| Inversion (documented expectation) | exists (T12 / `hor_event_inversion`) | either stays `HOR` (false-mixed cost documented) or document the new expectation explicitly |
| Low-coverage: N-dominated blocks | synthetic (large N regions) | stays `Class::HOR` — coverage gate kicks in |
| Schema test: `segments.tsv` column count | post-M7 | 13 columns (was 11) |
| Schema test: `diagnostics.json::schema_version` | post-M7 | 2 |

## 7. Sub-milestones, sized for review

| | Description | Code LOC est | Wall est | Risk |
|---|---|---:|---:|---|
| M7.1 | Analysis blocks + per-block consensus + identity stats (no class change) | ~200 | 0.5 d | Low |
| M7.2 | Mixed override + new `Segment` columns + locked tests | ~150 | 0.5 d | Med — tuning |
| M7.3 | Schema bump + segment-consensus FASTA for mixed | ~80 | 0.5 d | Low |
| M7.4 | A19 amendment + rerun reports + dashboard | ~0 (md only) | 0.5 d | Low |

Total ≈ 2 days of focused work, four PRs. Each PR runs the full
test suite + benchmark eval; M7.2 has a calibration loop on
`stratification_*` thresholds.

## 8. Out of scope for M7 (deferred to M8+)

- Per-segment class assignment (Q6 option b — heavy per-block
  `classify::decide_array`)
- Adaptive segment boundaries via changepoint detection (Q1 option b)
- HMM-based segmentation (§0 Q5, deferred to v2)
- Strand-aware inversion (OQ3)
- Real `inter_monomer_identity` computation (review #5)
- Nested-HOR (T09) support
- `simple_TR`-as-mixed taxonomy change (review-blocking #2;
  needs joint update to `taxonomy.md` + `T15` + expectations)

## 9. Review history

| Date | Document | Outcome |
|---|---|---|
| 2026-05-16 | `docs/new/detect_m7_plan.md` (DRAFT) | Original 11-question plan posted |
| 2026-05-16 | `docs/reviews/detect_m7_plan_review_2026-05-16.md` | All Q1–Q11 resolved; two scope changes: (a) analysis blocks vs reported segments split; (b) `simple_TR` mixed override held out |
| 2026-05-16 | This document | DRAFT → REVIEWED, decisions folded in |

## 10. References

- `docs/reviews/detect_implementation_review_completed_2026-05-16.md`
  findings #2 (segment recompute) and #3 (same-width mixed).
- `docs/reviews/detect_m7_plan_review_2026-05-16.md` — review of this
  plan's first draft; resolutions inlined above.
- `docs/new/detect_impl_plan.md` §0 A16 (the explicit deferral).
- `docs/new/taxonomy.md` — pin for the `simple_TR` stratification
  policy referenced in Q5.
- `docs/reports/kite_to_detect_v2_2026-05-16.md` (the 78 pp mixed
  regression under kite-derived periods).
- `tests/detect_expectations.tsv` — pins the `T15_stratification`
  expectation referenced in Q5.
- `src/detect/segment.rs` (current MVP — extended, not replaced).
- `src/detect/config.rs` (where `stratification_*_threshold` live;
  new fields land here in M7.1).

# Detect M7 Plan Review

Date: 2026-05-16  
Plan reviewed: `docs/new/detect_m7_plan.md`

## Executive Take

M7 is the right next detector milestone. The plan targets a real
failure mode: same-width mixed arrays collapse to one `HOR` or
`simple_TR` call when upstream period evidence does not expose the
block structure. I would proceed with M7, but I would tighten the
design before coding. The main changes I recommend are:

1. Keep fixed blocks as internal analysis blocks, not automatically
   reported biological segments.
2. Resolve the `simple_TR` stratification policy before enabling the
   mixed override for `simple_TR`.
3. Compare the right consensus width for HORs, ideally the HOR-unit
   consensus when enough complete units exist.
4. Make pair aggregation, low-coverage identity handling, and output
   schema explicit.
5. Add the missing segment-cap configuration instead of treating it as
   an implicit implementation detail.

## Blocking Design Decisions

### 1. Do Not Conflate Analysis Blocks With Reported Segments

The plan proposes fixed 100-row blocks when no phase shifts exist. That
is fine for internal detection, but reporting every fixed block in
`segments.tsv` would change the meaning of `n_segments` and create many
segments for clean single-family arrays.

Recommendation:

- Introduce an internal `SegmentBlock` or `AnalysisBlock` for M7 mixed
  detection.
- Emit rows in `segments.tsv` only for biologically meaningful reported
  segments: phase-shift segments, mixed-call blocks, or other localized
  events.
- Keep `Properties.n_segments = 1` for clean single-family arrays with
  no reported segmentation.
- If diagnostics need the internal block count, add it to
  `diagnostics.json` rather than overloading `n_segments`.

This preserves current output semantics while still enabling the
same-width mixed gate.

### 2. `simple_TR` Eligibility Conflicts With Current Taxonomy

The plan proposes that `simple_TR` should be eligible for the same
consensus-identity mixed override. That conflicts with the current
fixture and taxonomy:

- `tests/detect_expectations.tsv` expects `T15_stratification` to remain
  `simple_TR`.
- `tests/synth_configs/T15_stratification.yaml` is two same-period
  `SIMPLE_TR` blocks with different monomer templates.
- `docs/new/taxonomy.md` treats same-period stratification as a property,
  not a class change.

Recommendation: for M7, do not enable the mixed override for
`simple_TR` by default. Apply it to `HOR` and `irregular_HOR` only,
where same `base_width` plus same `k` can still hide different HOR
families. If the desired taxonomy is changing so that same-period
divergent `simple_TR` blocks become `mixed`, update `taxonomy.md`,
`T15_stratification`, and `detect_expectations.tsv` before M7 code lands.

## Recommended Answers To Q1-Q11

| Question | Recommendation |
|---|---|
| Q1 segmentation | Use deterministic fixed analysis blocks plus phase-shift boundaries, but make block size adaptive: `block_rows = max(block_size_rows_min, ceil(n_rows / max_segments_per_array))`. Merge or skip blocks below the minimum informative row count. Do not emit every block as a reported segment. |
| Q2 consensus | Majority-vote consensus is appropriate for M7. For `simple_TR`, use `base_width_bp`. For `HOR` / `irregular_HOR`, prefer `hor_length_bp` consensus when enough complete HOR units exist; fall back to `base_width_bp` only when unit-level consensus is not reliable. |
| Q3 identity metric | Hamming identity is the right MVP metric, but it should ignore positions where either consensus has `N`. Return `None` when informative coverage is too low, for example `< 70%` of positions. |
| Q4 thresholds | Keep `same >= 0.90` and `different <= 0.80` as starting defaults. Use all valid pairwise comparisons, not only segment 1. The plan text says "any pair", but the pseudocode compares only to segment 1; align those. Borderline pairs should be diagnostic only, not an automatic `mixed`. |
| Q5 ordering | Run the mixed override after initial classification and before irregular-HOR demotion. Apply it to `HOR` / `irregular_HOR` in M7. Keep `simple_TR` out unless the taxonomy decision above is changed. |
| Q6 per-segment class | Use the light option for M7. Do not run full per-segment classification. If rows are emitted, set `Segment.class` from the final array class and document that it is not an independent per-segment classifier result. |
| Q7 schema | Add an identity column, but prefer `consensus_identity_to_reference` over `consensus_identity_to_seg1` if the implementation uses a medoid/reference block. Also consider `consensus_identity_coverage` because low-coverage identity must be distinguishable from high-confidence identity. Bump diagnostics schema to 2. Emit segment consensus FASTA only for mixed calls. |
| Q8 performance | Add explicit config fields for the cap, for example `max_segments_per_array = 32` and possibly `min_segment_rows`. Skip per-segment R(k) in M7. Since array-level detect already uses Rayon, do not add nested segment-level parallelism until profiling shows it helps. |
| Q9 acceptance | The `mixed >= 70%` kite-derived floor is reasonable for shipping, but also require oracle-period mixed performance to stay near the current 96%, for example `>= 94%` or no more than a 2 pp drop. Add false-mixed guards for clean HOR, wobble HOR, inversion, and the chosen T15 stratification policy. |
| Q10 kite changes | Agree: ship M7 without kite-side changes. Re-run kite -> detect after M7 and only then decide whether kite period emission still needs adjustment. |
| Q11 diagnostics | Agree: bump `diagnostics.json` `schema_version` to 2 and update schema-drift tests plus downstream report readers. |

## Implementation Notes

The current `segment::split(&Properties)` API cannot implement M7 because
it has no sequence, config, selected width, or candidate context. The
mixed override also needs to run before the final class, confidence, and
reason are finalized. I would refactor this as a small analysis step
inside `run_array_m4` or just after the first classification pass:

1. Build internal analysis blocks from the chosen comparison width.
2. Compute consensus and identity stats for those blocks.
3. Decide whether to override to `mixed`.
4. Then run irregular demotion, wobble guard, confidence, and output
   segment-row construction from the final `Properties`.

Avoid recomputing the same consensus work again in `run_one` when writing
outputs. Either return the segment analysis alongside `Properties` or
store enough final segment records to write once.

For HORs, block boundaries need extra care. If comparing HOR-unit
consensuses, blocks should be aligned to full HOR units or at least drop
partial units at block edges. A boundary-straddling fixed block should
not be allowed to drive a `mixed` call by itself.

## Test Coverage To Add

- Same-width, same-`k` two-block HOR that currently collapses to `HOR`
  and should become `mixed`.
- A clean same-width HOR with high within-family divergence that must not
  become `mixed`.
- Wobble and phase-shift HOR negative controls.
- Inversion negative control, or an explicit documented expected failure
  if strand-awareness remains deferred.
- `T15_stratification` locked to the chosen policy.
- Low-coverage consensus identity where many positions are `N`; this
  should not trigger `mixed`.
- Schema tests for `segments.tsv`, `consensus.fa`, and diagnostics schema
  version 2.

## Suggested Revised M7 Scope

I would define M7 as detector-side, HOR-focused same-width mixed
detection:

- M7.1: add internal analysis blocks, per-block consensus, and identity
  stats without changing class calls.
- M7.2: add the `HOR` / `irregular_HOR` mixed override using all valid
  pairwise identities.
- M7.3: add reported segment rows and segment consensus FASTA only for
  final `mixed` calls; bump diagnostics schema to 2.
- M7.4: rerun oracle and kite-derived v2 evaluations and update reports.

This keeps the milestone focused on the measured 78 pp mixed regression
while avoiding a silent taxonomy change for same-period `simple_TR`
stratification.

# Detect M7 — per-segment recompute + same-width mixed detection

Status: **DRAFT** (2026-05-16). Sign-off needed before any code lands.

## TL;DR

M7 closes the two architectural deferrals from
`docs/reviews/detect_implementation_review_completed_2026-05-16.md`
(findings #2 and #3), tracked as A16 in
`docs/new/detect_impl_plan.md` §0:

- **Per-segment recompute** — `segments.tsv` rows currently inherit
  the whole-array class / base_width / k / IC / phase_sep / wobble.
  M7 replaces that with genuine per-segment feature extraction so
  segment fields mean what they say.
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
- **Non-MVP segment classes** (per-segment regime A/B/C). Each
  segment still inherits the whole-array class for M7; only the
  *features* are per-segment.

## 2. Open questions

Each question has a **proposed default** I'd ship with. Items in
**bold** want explicit sign-off before code lands.

### Q1. Segmentation criterion — how do we split when there are no phase shifts?

| Option | Description | Cost | When right |
|---|---|---|---|
| **(a)** | Fixed-stride blocks at `block_size_rows_min` (100 rows). Append phase-shift boundaries as extra splits. | Cheap | MVP — deterministic, no thresholds to tune |
| (b) | Detect changepoints in row-to-row R(1). Adaptive segment boundaries. | Medium | If (a) misses architectural transitions inside a block |
| (c) | First-half vs second-half only. | Cheapest | Too coarse — misses 3-block mixed cases |
| (d) | HMM segmentation. | High | Already deferred to v2 (§0 Q5) |

**Proposed: (a).** Fixed 100-row blocks + phase-shift boundaries.
Simple, deterministic, exposes the stratification thresholds.
Downside: blocks straddling a transition will dilute both sides.
We can add (b) later as a refinement without changing the
interface.

### Q2. Consensus computation per segment

Two interpretations:

| Option | Description |
|---|---|
| **(a)** | Majority-vote consensus at `base_width`, ignoring N (current `consensus::consensus()` behaviour, just called per-segment slice). |
| (b) | Per-column profile (4-vector probabilities), enabling profile-based similarity. |

**Proposed: (a).** Hamming identity on majority-vote consensus is
both interpretable and matches the spec's "consensus identity"
language. (b) is a refinement we can take later if (a) is too
coarse on noisy segments.

### Q3. Identity metric for comparing segment consensuses

| Option | Description |
|---|---|
| **(a)** | Hamming identity over the consensus bytes (1 − mismatch/length). |
| (b) | Edit distance (handles indel-drifted segments). |
| (c) | k-mer composition (R(1)) at consensus level. |

**Proposed: (a).** All segments share the same `base_width`, so
their consensuses are the same length. Hamming is O(L), trivial
to test, and what the user reads when they look at the values.

### Q4. Two thresholds, two decisions

`DetectorConfig` already exposes:
- `stratification_same_threshold` (0.90) — "considered same family"
- `stratification_diff_threshold` (0.80) — "considered different family"

Three-way decision per pair:
- identity ≥ 0.90 → same family
- identity ≤ 0.80 → different family
- 0.80 < identity < 0.90 → undecided

**Proposed pair-aggregation rule for mixed call:** if *any* pair of
segments has identity ≤ `stratification_diff_threshold` → mixed.
Borderline (between same/diff thresholds) doesn't itself trigger
mixed — too noisy. We could revisit if benchmark shows borderline
cases as the common confusion.

### Q5. **Class-level interaction order — when does mixed fire vs irregular_HOR?**

Today's order in `mod.rs::run_array_m4`:
1. classify (HOR / simple_TR / mixed / ambiguous)
2. irregular_HOR demotion if `irregularity_score ≥ 0.50`
3. wobble-dominance guard (kept since 2026-05-16)

Proposed M7 order:
1. classify (as today)
2. **segment recompute + consensus identity → mixed override** *(NEW)*
3. irregular_HOR demotion
4. wobble-dominance guard

i.e., a single-family HOR with high block-IC variance keeps the
old "irregular_HOR" demotion path; an HOR whose segments
genuinely disagree on consensus becomes `mixed` first and bypasses
irregularity entirely.

**Open:** should `simple_TR` also be eligible for the mixed
override? Two simple_TR families at different monomer sequences
but same period is the simple_TR analogue of the same-width HOR
mixed case. **Proposed: yes**, gated by the same identity test.

### Q6. Per-segment class assignment

Plan §6.9 nominally calls for per-segment regime A/B/C
classification. Two interpretations:

| Option | Description | M7 fit |
|---|---|---|
| **(a)** | Light: every segment inherits the whole-array class. Per-segment fields (IC, phase_sep, wobble, irregularity) are real, but `Segment.class` ≡ `Properties.class`. | Yes |
| (b) | Heavy: each segment runs its own `classify::decide_array` with its slice of periods. Segment classes can differ. | M8 |

**Proposed: (a) for M7.** Heavy per-segment classify multiplies
the per-array runtime cost and we don't have a clear use case yet
— the array-level class is usually right. The light option still
delivers the headline win (same-width mixed detection) and the
honest per-segment features.

### Q7. **Output schema impact**

The `segments.tsv` schema already has the per-segment fields
(`base_width_bp`, `hor_k`, `column_conservation`, `phase_separation`,
`wobble_amplitude_bp`, `irregularity_score`). Today they're filled
with whole-array values; M7 fills them with per-segment values.
That's a **semantic** change to existing columns, not a schema bump.

New questions:

1. Do we add `consensus_identity_to_seg1` to `segments.tsv`? It's
   the load-bearing M7 number, useful for diagnostics.
2. Do we write per-segment consensus FASTAs to `consensus.fa`?

**Proposed:**
- Add `consensus_identity_to_seg1` as a new column. Requires
  `SEGMENTS_HEADER` and the test that asserts column count to be
  bumped. Bump `diagnostics.json schema_version` from 1 → 2.
- Write per-segment consensus to `consensus.fa` *only when the
  array is called `mixed`*. Naming: `<array_id>_seg<N>_monomer`.
  Single-family arrays continue to emit only the whole-array
  monomer + hor_unit (no segment proliferation in the output).

### Q8. Performance budget

Per-segment recompute is roughly O(n_segments × L_segment × max_k).
A 50 000-row array at 100-row segments = 500 segments, each ~17 kb,
each R(k) ~O(L × 30) → ~25 ms × 500 = 12 s per array. Times 1600
arrays = 5 h. Unacceptable.

**Proposed mitigations** (need sign-off):
- Cap `n_segments_per_array` (default 32). Coarser blocks for
  large arrays.
- Skip per-segment R(k) when segments inherit the whole-array
  class anyway (Q6 option a). We only need wrap + column IC +
  consensus per segment for the identity test — that's O(L) per
  segment, total O(N) per array. Should run in tens of ms.
- Rayon-parallelise the segment loop within each array.

If we adopt (a) + (a) of Q6 + skip per-segment R(k), the
per-array overhead should be sub-100 ms.

### Q9. Backward compatibility with the M6 baseline (94.4 %)

Adding mixed detection may flip currently-correct calls:

- HOR cases with mild within-array divergence may newly cross
  the `diff_threshold` and become mixed (false mixed).
- mixed cases currently called HOR may correctly become mixed
  (true mixed — the intended win).

**Proposed acceptance gate**:
- `mixed` category ≥ 70 % on the v2 corpus (vs 96 % oracle / 18 %
  kite). Bumping to 90 % is the target but we accept 70 % as the
  "M7 ships" floor.
- No other category drops by > 2 pp from the M6 baseline.
- Core CI fixtures (T01, T05, T06, T07, T10, T13, T17, T18) still
  pass exactly.

If those bounds aren't reachable, we re-tune `stratification_*`
thresholds before merging.

### Q10. **Same-width mixed via kite-derived periods**

The kite → detect run report shows the mixed regression is
*structural* — kite emits ≤ 2 strong periods per array. M7 fixes
the mixed regression *without changing the kite emitter*, because
the new test runs on whole-array sequence, not on the period
candidate list.

**Proposed:** ship M7 with no kite-side changes; re-run kite →
detect after M7 lands and update `docs/reports/`.

### Q11. **Diagnostics JSON schema version bump**

`detect/io.rs::write_diagnostics()` currently writes
`"schema_version": 1`. M7 changes per-segment fields semantically
and adds at least one column. **Proposed: bump to 2** and document
the change in the impl plan.

## 3. Implementation plan (draft)

### 3.1 Modules touched

| File | Change |
|---|---|
| `src/detect/segment.rs` | Rewrite. Add per-segment wrap + column IC + consensus + identity to seg1. |
| `src/detect/types.rs` | Add `Segment.consensus_identity_to_seg1: Option<f64>`. Bump `SEGMENTS_HEADER`. |
| `src/detect/mod.rs::run_array_m4` | Insert step "segment recompute + mixed override" between classify and irregular demotion. |
| `src/detect/classify.rs` | Expose a small helper `consider_mixed_via_segments()` so the override is callable from `mod.rs`. |
| `src/detect/consensus.rs` | Add `consensus_on_slice()` for per-segment use (or reuse `consensus()` with byte slices — already supports it). |
| `src/detect/config.rs` | (No change — `stratification_*` already exist.) |
| `src/detect/io.rs` | Add per-segment consensus to `consensus.fa` when class=mixed. Bump diagnostics schema_version to 2. |
| `docs/new/detect_impl_plan.md` | New A19 amendment; mark M7 as in-progress in §10. |
| `tools/detect_eval/` | Re-run + update reports. |

### 3.2 Sub-milestones (PR-sized)

**M7.1 — Per-segment features, no class change.**
- Refactor `segment::split()` to compute per-segment wrap + IC + R(1) + consensus.
- Fill `Segment` fields with per-segment values (no schema change).
- Tests: per-segment IC matches whole-array IC when only one segment;
  per-segment IC differs from whole-array IC when two segments differ.
- **Acceptance**: existing 264 lib + integration tests still pass;
  no benchmark regression.

**M7.2 — Consensus identity column + mixed override.**
- Add `consensus_identity_to_seg1` to `Segment` and `SEGMENTS_HEADER`.
- In `run_array_m4`, after classify and before irregularity demotion,
  if min pairwise identity ≤ `stratification_diff_threshold` →
  rewrite class to `Mixed` (with reason explaining which segments
  diverged).
- Tests: synthesise a 2-block mixed config (same width, same k,
  different divergence seed); assert detector calls `mixed`.
- **Acceptance**: mixed category accuracy on `ground_truth_v2/`
  improves to ≥ 70 % (kite-driven baseline 18 %); no other
  category drops > 2 pp.

**M7.3 — Wire stratification thresholds + segment-consensus FASTA.**
- Use `stratification_same_threshold` and `_diff_threshold` from
  config (currently inert).
- Emit per-segment consensus records to `consensus.fa` only when
  class=mixed.
- Bump diagnostics `schema_version` to 2.
- **Acceptance**: schema-drift tests bumped, regen of CLAUDE.md +
  README + impl plan as a new A19 amendment.

**M7.4 — Documentation + rerun reports.**
- Update `docs/new/detect_impl_plan.md` with A19 and mark M7 done
  in §10.
- Rerun kite → detect on v2 corpus; update
  `docs/reports/kite_to_detect_v2_*.md` with the new numbers (or
  write a fresh report dated post-M7).
- Re-render the dashboard.

### 3.3 New `Segment` row

```rust
pub struct Segment {
    pub array_id: String,
    pub segment_id: usize,
    pub start_bp: usize,
    pub end_bp: usize,
    pub class: Class,           // M7.1: stays whole-array class
    pub base_width_bp: Option<usize>,
    pub hor_k: Option<usize>,
    pub column_conservation: Option<f64>,        // NEW: per-segment IC
    pub phase_separation: Option<f64>,           // NEW: per-segment phase_sep
    pub wobble_amplitude_bp: Option<f64>,        // NEW: per-segment wobble
    pub irregularity_score: Option<f64>,         // already had this
    pub consensus_identity_to_seg1: Option<f64>, // NEW: M7.2
}
```

(The first 8 fields already exist in the struct — only
`consensus_identity_to_seg1` is genuinely new.)

### 3.4 Pseudocode for the mixed override

```rust
fn mixed_override_via_segments(
    seq: &[u8],
    base_width: usize,
    boundaries: &[usize],          // segment boundaries in rows
    cfg: &DetectorConfig,
) -> Option<MixedDecision> {
    if boundaries.len() < 2 { return None; }
    let consensuses: Vec<Vec<u8>> = boundaries
        .windows(2)
        .filter_map(|w| consensus::consensus_on_slice(seq, base_width, w[0], w[1]))
        .collect();
    if consensuses.len() < 2 { return None; }
    let seg1 = &consensuses[0];
    let identities: Vec<f64> = consensuses[1..]
        .iter()
        .map(|c| hamming_identity(seg1, c))
        .collect();
    let min_id = identities.iter().cloned().fold(f64::INFINITY, f64::min);
    if min_id <= cfg.stratification_diff_threshold {
        return Some(MixedDecision {
            reason: format!(
                "mixed — segment-consensus identity {:.3} ≤ diff_threshold {:.3}",
                min_id, cfg.stratification_diff_threshold
            ),
            n_segments_diverging: identities.iter().filter(|&&i| i <= cfg.stratification_diff_threshold).count() + 1,
        });
    }
    None
}
```

## 4. Risks

| Risk | Mitigation |
|---|---|
| False mixed on noisy single-family arrays (hor_clean d40+, hor_wobble) | `diff_threshold` defaulted to 0.80 — well below `same_threshold` 0.90 — so only genuinely divergent consensuses trigger. Tune in M7.2 acceptance. |
| Per-segment recompute slow on large arrays | Q8 mitigations. Cap n_segments per array. |
| Inversion segments look mixed (RC consensus ≠ forward consensus) | Inversions already deferred to v2 strand-aware. Add an explicit "skip mixed override when array is hor_event_inversion" path? Probably not — let the false-mixed count be the cost of deferring strand awareness. Document in the reason field. |
| Borderline (0.80 < id < 0.90) cases that should be mixed | Q4. Tune thresholds in M7.2 if benchmark shows this is the common confusion. |
| Schema bump in diagnostics.json breaks downstream eval | Eval uses array_id + class only; schema bump should be transparent. Verify with `tools/detect_eval/eval.py` before/after. |

## 5. Acceptance criteria (M7 done when)

1. `mixed` category accuracy ≥ 70 % on `ground_truth_v2/` under both:
   - Oracle periods (target ≥ 90 %; current 96 %, expect ~unchanged
     or slightly higher with the new signal).
   - Kite-emitted periods (target ≥ 70 %; current 18 %, the kite →
     detect report regression is the load-bearing motivation).
2. No category drops by > 2 pp from the M6 baseline.
3. Core CI fixtures (T01, T05, T06, T07, T10, T13, T17, T18) still
   pass exactly via `cargo test --release --test detect_m4`.
4. `cargo test --release` green; new unit tests for
   `consensus_on_slice`, `hamming_identity`, segment recompute,
   and the synthesised 2-block mixed fixture.
5. `tests/detect_kite_emit.rs` still passes; no schema-related
   regression.
6. `docs/new/detect_impl_plan.md` updated with A19 amendment and
   M7 marked done in §10.
7. `docs/reports/kite_to_detect_v2_<post-m7-date>.md` written
   with the new numbers and an updated dashboard.

## 6. Sub-milestones, sized for review

| | Description | Code LOC est | Wall est | Risk |
|---|---|---:|---:|---|
| M7.1 | Per-segment feature recompute (no class change) | ~150 | 0.5 d | Low |
| M7.2 | Identity column + mixed override | ~120 | 0.5 d | Med — tuning |
| M7.3 | Stratification thresholds + segment-consensus FASTA + schema_version bump | ~80 | 0.5 d | Low |
| M7.4 | Docs + rerun reports | ~0 (md only) | 0.5 d | Low |

Total ≈ 2 days of focused work, four PRs. Each PR runs the full
test suite + benchmark eval; M7.2 has a calibration loop on
`stratification_*` thresholds.

## 7. Out of scope for M7 (deferred to M8+)

- Per-segment class assignment (Q6 option b)
- Adaptive segment boundaries via changepoint detection (Q1 option b)
- HMM-based segmentation (§0 Q5, deferred)
- Strand-aware inversion (OQ3)
- Real `inter_monomer_identity` computation (review #5)
- Nested-HOR (T09) support

## 8. Things I want sign-off on before coding

The bolded **(B)** open questions:

- **Q5**: irregularity ordering and whether simple_TR is mixed-eligible.
- **Q7**: schema bump scope — new column + segment consensuses in FASTA.
- **Q8**: performance mitigations (cap n_segments, skip per-segment R(k)).
- **Q10**: ship M7 with no kite-side changes; reassess emit_periods after.
- **Q11**: diagnostics.json `schema_version` 1 → 2.

For each I've proposed a default. If any of those defaults are
wrong, I'd rather find out before M7.1 lands.

## 9. References

- `docs/reviews/detect_implementation_review_completed_2026-05-16.md`
  findings #2 (segment recompute) and #3 (same-width mixed).
- `docs/new/detect_impl_plan.md` §0 A16 (the explicit deferral).
- `docs/reports/kite_to_detect_v2_2026-05-16.md` (the 78 pp mixed
  regression under kite-derived periods).
- `src/detect/segment.rs` (current MVP — replaced in M7.1).
- `src/detect/config.rs` (where `stratification_*_threshold` live).

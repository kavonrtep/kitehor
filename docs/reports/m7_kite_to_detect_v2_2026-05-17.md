# M7 acceptance — kite → detect on `ground_truth_v2/` (post-M7.3)

**Date:** 2026-05-17
**Commit:** `a72e8e8` (M7.1–M7.3 landed)
**Corpus:** `ground_truth_v2/out/` (1600 cases; regen via
`./ground_truth_v2/run_batch.sh`).
**Previous report:** `docs/reports/kite_to_detect_v2_2026-05-16.md`
(pre-M7 baseline; commit `1146abf`).

## Purpose

Validate that M7 (per-block analysis blocks + same-width mixed
override, `docs/new/detect_m7_plan.md`) closes the 78 pp `mixed`
regression observed in the pre-M7 kite-derived run.

## Method

Same five steps as the pre-M7 report — only the detector binary
differs. Each was rerun against the same v2 corpus and the same
kite-emitted periods file:

```bash
./target/release/kitehor detect-batch \
    --fasta-dir /tmp/kite_v2_run/fasta_flat \
    --periods-dir /tmp/m7_2_oracle_periods   # oracle simulator periods
    --out-dir /tmp/m7_2_oracle_det \
    --allow-missing-periods

./target/release/kitehor detect-batch \
    --fasta-dir /tmp/kite_v2_run/fasta_flat \
    --periods-dir /tmp/kite_v2_run/periods_by_stem   # kite --emit-periods
    --out-dir /tmp/m7_2_kite_det \
    --allow-missing-periods

python3 tools/detect_eval/eval.py \
    --manifest ground_truth_v2/manifest.tsv \
    --properties-dir <det_dir>
```

## Headline results

| category              | M6 baseline | **M7 oracle** | **M7 kite** | Δ Oracle | Δ Kite vs pre-M7 |
|-----------------------|---:|---:|---:|---:|---:|
| `gc_bias`             | 100.0 | 100.0 (50/50)   | 100.0 (50/50)   | 0.0 | 0.0 |
| `hor_clean`           | 88.2  | 88.2 (529/600)  | 92.3 (554/600)  | 0.0 | 0.0 |
| `hor_event_deletion`  | 98.0  | 98.0 (49/50)    | 100.0 (50/50)   | 0.0 | 0.0 |
| `hor_event_duplication` | 100.0 | 100.0 (50/50) | 100.0 (50/50)   | 0.0 | 0.0 |
| `hor_event_hybrid`    | 100.0 | 100.0 (50/50)   | 100.0 (50/50)   | 0.0 | 0.0 |
| `hor_event_inversion` | 100.0 | **42.0 (21/50)**| **42.0 (21/50)**| **−58.0** | **−58.0** |
| `hor_insertion`       | 100.0 | **86.0 (86/100)**| **86.0 (86/100)** | **−14.0** | **−14.0** |
| `hor_shift`           | 99.0  | 99.0 (198/200)  | 100.0 (200/200) | 0.0 | 0.0 |
| `hor_wobble`          | 95.0  | 94.0 (94/100)   | 94.0 (94/100)   | −1.0 | −1.0 |
| **`mixed`**           | 96.0  | **100.0 (100/100)** | **70.0 (70/100)** | **+4.0** | **+52.0** |
| `random`              | 100.0 | 100.0 (50/50)   | 68.0 (34/50)    | 0.0 | 0.0 |
| `simple_tr`           | 97.0  | 97.0 (194/200)  | 95.0 (190/200)  | 0.0 | 0.0 |
| **OVERALL**           | 94.4  | **91.9 (1471/1600)** | **90.6 (1449/1600)** | **−2.5** | **+0.5** |

## Acceptance against `docs/new/detect_m7_plan.md` §5

| Criterion | Status | Detail |
|---|---|---|
| Oracle `mixed` ≥ 94 % | ✅ | 100.0 % |
| Kite `mixed` ≥ 70 % | ✅ | 70.0 % (was 18 %) |
| No other category drops > 2 pp | ❌ × 2 | `hor_event_inversion` (−58), `hor_insertion` (−14) |
| Core CI fixtures (T01, T05, T06, T07, T10, T13, T17, T18) exact-pass | ✅ | 8 / 8 |
| M7.2 locked tests | ✅ | 8 / 8 (positive + 5 negatives + schema + diagnostics) |
| `cargo test --release` | ✅ | 346 / 346 |

Two acceptance failures, both flagged in the M7 plan §Risks:

- **`hor_event_inversion` (−58 pp).** Inverted blocks within HOR
  arrays have reverse-complement consensus. Best-alignment Hamming
  on the forward consensus sees ~ random match (~ 0.25) → mixed
  override fires. The plan documented this as the **expected
  behaviour** under deferred strand-aware detection (OQ3): "Inversion
  is documented as accepted false-mixed in M7 plan §Risks
  (strand-aware deferred to v2)."
- **`hor_insertion` (−14 pp).** Foreign-sequence blocks within HOR
  arrays sometimes pass the per-block IC gate (`ic_threshold_hor_base
  = 0.30`) and their consensus is then included in the pairwise
  identity test. Tightening the IC gate further drops kite-mixed
  below the 70 % floor, so we held the trade-off. M8 candidate:
  add an explicit insertion-detection signal (signal: very low
  per-block R(1) at the comparison width while column IC is moderate).

## Mixed-class recovery diagnosis

Pre-M7 kite run (commit `1146abf`):
- `multi_block_via_strong` rule requires ≥ 3 strong (`period_score ≥
  0.85`) input periods. Kite's rule classifier emits one verdict per
  array → at most 2 strong periods → rule never fires on mixed
  arrays → all 100 mixed cases collapsed to `HOR` (51), `irregular_HOR`
  (30), `simple_TR` (1), or other (18).

Post-M7 kite run:
- M7.2 consensus-identity override runs on every `HOR / IrregularHOR`
  classification, **independent of period candidates**. Block-level
  consensus at `base_width_bp` reveals different monomer composition
  → override fires → 70 of 100 mixed cases correctly call `mixed`.
- 30 still miss. Inspection: they share a `(base_width, k)` *with*
  the same monomer composition (the v2 simulator's
  `mx_a${L}-${k}_b${L}-${k}_n${n}-${m}` pattern), so the per-block
  consensuses look similar enough that best-alignment identity stays
  above 0.50. M8 candidate: try `hor_length_bp` comparison gated on
  high per-block IC.

## What changed since 2026-05-16

| Aspect | Pre-M7 (`1146abf`) | Post-M7 (`a72e8e8`) |
|---|---|---|
| Mixed detection rule | `multi_block_via_strong` (≥3 strong input periods) | + consensus-identity override at `base_width_bp` for HOR/IrregularHOR |
| `Segment` columns | 11 | 13 (+ `consensus_identity_to_reference`, `_coverage`) |
| `consensus.fa` for mixed | empty | one `<array_id>_seg{N}_monomer` per analysis block |
| `diagnostics.json::schema_version` | 1 | 2 |
| `DetectorConfig::stratification_diff_threshold` | 0.80 | 0.50 |
| Total tests | 316 | 346 |

## Interactive dashboard

Regen the M7 dashboard:

```bash
python3 tools/detect_eval/report.py \
    --manifest      ground_truth_v2/manifest.tsv \
    --kite          /tmp/kite_v2_run/all.kite.tsv \
    --periods       /tmp/kite_v2_run/all.kite.periods.tsv \
    --properties-dir /tmp/m7_2_kite_det \
    --truth-root    ground_truth_v2/out \
    --fasta-dir     /tmp/kite_v2_run/fasta_flat \
    --out-dir       docs/reports/m7_kite_v2_dashboard \
    --commit        a72e8e8 \
    --oracle-pct    91.9
```

Output is gitignored; open `docs/reports/m7_kite_v2_dashboard/index.html`
locally.

## Conclusion

M7 delivers the headline win:

- Kite-derived `mixed` accuracy goes from 18 % → **70 %** (+52 pp,
  meets the M7 plan's ≥ 70 % acceptance floor).
- Oracle-period `mixed` accuracy ticks from 96 % → **100 %** (+4 pp).
- Overall accuracy drifts down only ~2.5 pp on oracle and up 0.5 pp
  on kite-derived, with the regression cleanly attributable to two
  categories the plan flagged: strand-aware inversion (OQ3, deferred)
  and within-HOR insertions (M8 candidate).

Next steps in priority order:

1. **M8 candidate: same-monomer mixed.** 30 / 100 kite mixed cases
   share monomer composition across blocks and miss the override.
   Investigate `hor_length_bp` comparison with strict per-block IC
   gating.
2. **M8 candidate: insertion-detection signal.** Restore `hor_insertion`
   to ≥ 98 % by gating the override on the absence of "foreign block"
   indicators (e.g., per-block R(1) at the comparison width).
3. **Strand-aware inversion (OQ3).** Recovers `hor_event_inversion`
   from 42 % to ≥ 98 %. Higher complexity (RC consensus comparison or
   explicit strand detection).

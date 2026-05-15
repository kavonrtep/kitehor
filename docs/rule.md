# Rule-based HOR classifier

Last updated: 2026-05-15.

`kitehor`'s default HOR call (`--classify`) is a four-condition rule
applied directly to the kite peak output. Every peak in the kite output
has already passed kite's `peak > background` and `score2_norm > 0.001`
filters, so we trust each peak as a real periodicity and only ask
whether the structure of peaks forms a Higher-Order Repeat.

## The rule

```text
HOR ⟺ d1 (strongest kite peak)
       = k × p_n
       for some integer k ∈ [2, qmax]
       and p_n a kite peak in the top-N by score
       and |d1 − k·p_n| ≤ max(tol_bp, tol_rel × k·p_n)
       and d1   ≥ lo_period
       and p_n  ≥ lo_period
```

If the rule fires: `verdict = hor`, `founder = p_n`, `tile = d1`,
`multiplicity = k`. If no match exists: `verdict = tandem` with
`monomer_bp = d1`. If kite returned no peaks: `verdict = no_signal`.
If d1 is below `lo_period`: `verdict = unresolved`.

## Defaults

| Parameter | Default | CLI flag |
|---|---:|---|
| `top_n`     | 3   | `--rule-top-n` |
| `qmax`      | 30  | `--rule-qmax` |
| `tol_bp`    | 5   | — |
| `tol_rel`   | 0.02 | — |
| `lo_period` | 15  | — |

`top_n = 3` is the user-validated value: it filters deep sub-period
harmonics out of the founder-candidate pool (e.g. lag-87/91 inside
AT178), which were the dominant source of false positives. Higher
values re-admit those spurious matches; lower values lose real HORs.

## Why unidirectional (`d1` must be the tile)

The rule only fires when `d1` is the *tile* and the founder is a
*weaker* peak at `d1/k`. The converse arrangement (`d1` = founder, with
a weaker peak at `d1·k`) is **not** tested.

This is a deliberate choice driven by the empirical study (see
`/tmp/empirical/` outputs from the 2026-05-15 session): bidirectional
search re-introduces the TRC_1 over-sensitivity pattern (small dominant
monomer with weak higher-order harmonics), which the user explicitly
asked to avoid. The cost is that two real-data HOR shapes are missed:

1. Very high multiplicity HORs (e.g. AT TRC_2's k=29-30 cases) where
   the founder peak is structurally stronger than the tile because the
   tile sits at a very long period that kite's score normalisation
   deflates.
2. Single TRC_14_drapa record (`2282061_2308341`) where d1=375
   (founder) and d2=744 (tile) — same inverted shape as above.

User explicitly accepts these misses ("missing weak / unusual HOR is
OK; over-calling on weak signal is not").

## Why no `tile_share` (score-ratio) floor

An earlier version of the rule required
`min(s_founder, s_tile) / max(s_founder, s_tile) ≥ 0.20`. This was
dropped (2026-05-15) for two reasons:

1. **Real HORs with strong inter-position divergence naturally have low
   share.** When founders differ a lot across positions, the tile peak
   dominates and the founder peak is barely above background — yet
   biologically the HOR is *more*, not less, pronounced. Filtering on
   share penalises exactly the strongest evidence.
2. **It changed exactly one verdict** across the verification corpus.
   The `top_n` constraint, not the share floor, was actually doing the
   work of filtering spurious sub-period harmonics. Confirmed
   empirically on synth CV (n=1,204), AT full-length (n=155), and
   TRC_14_drapa (n=23).

## Comparison with the legacy ML classifier

The legacy random-forest + Platt + verdict-orchestrator pipeline is
still available via `--classify --use-ml-classifier`. It produces a
richer column set (`hor_score`, `hor_score_raw`, `k_pred`,
`recovered`, `h_d1`, `h_founder`) but was found to be **over-sensitive
on real centromeric arrays** (most notably AT TRC_1: model called HOR
on 9 records where the user judged ~7 should be `unresolved`) and
**under-sensitive on real HORs with strong inter-position divergence**
(TRC_14_drapa: 9 of 23 missed). The ML model was trained on synthetic
data with founder-dominant kite peaks, which is the *opposite* shape
of typical real centromeric HOR.

| | AT TRC_1 HOR | TRC_14_drapa correct | Notes |
|---|---:|---:|---|
| ML default | 9 (over) | 14/23 | Model trained on founder-dominant synth |
| Rule | **2** | **22/23** | Only the two user-confirmed legit cases on TRC_1 |

The ML path is retained for back-compat and is the right choice if you
are working with data drawn from the same distribution as the
synthetic training set (clean k=2–6 HOR, founder-dominant kite output,
moderate divergence).

## Output columns under the rule

When `--classify` runs with the default rule path, the output TSV adds
five columns:

| Column | Meaning |
|---|---|
| `verdict`       | `hor` / `tandem` / `unresolved` / `no_signal` |
| `founder`       | Inferred founder period (bp), only for `hor` |
| `multiplicity`  | k (1 for tandem) |
| `tile`          | HOR tile period (bp) |
| `share`         | `min(s_founder, s_tile) / max(...)`. Diagnostic only, NOT a filter. |

## Supplementary coverage QC (`--coverage`)

Pass `--coverage` together with `--classify` to add nine more columns
for records called `hor`. Records not called `hor` get `NA`.

For each HOR call the array is split into `floor(L / tile)` tile-length
windows; the first window is taken as the reference HOR unit and every
subsequent window is compared to it using **Levenshtein-based identity**
(`1 - edit_distance / max_len`). The aggregate columns are:

| Column            | Meaning |
|---|---|
| `cov_mean`        | mean identity to the first tile |
| `cov_pass_70/80/90` | fraction of windows with identity ≥ {0.70, 0.80, 0.90} |
| `cov_first_half`  | mean identity in the first half of windows |
| `cov_second_half` | mean identity in the second half — large gap flags mosaic arrays |
| `cov_min` / `cov_max` | worst / best window; `cov_min` flags single-tile dropouts |
| `cov_n_tiles`     | number of comparison windows |

What the score does and does not catch:

- ✓ **mosaic / partial-array** patterns (first_half vs second_half asymmetry)
- ✓ **tile dropouts** (`cov_min` outliers)
- ✓ **wrong-period calls** (`cov_mean` collapses toward 0.25 random)
- ✗ **founder = sub-period of a longer real monomer** — both produce
  high tile-aligned identity by construction. Distinguishing this case
  requires an external constraint such as a minimum-founder-length
  floor (not yet exposed as a flag).

This is purely supplementary — it does **not** enter the rule's
HOR/non-HOR decision. Cost: ~ms per record for typical centromeric
tiles; longer for kb-scale tiles. Empirically: synth seed 201 (1,204
records, 234 HOR calls) runs in **~8 s on 6 threads**.

## Provenance

The rule emerged from the 2026-05-15 architectural review after the ML
pipeline was found to mis-classify real centromeric arrays. The
preceding empirical study, post-processing experiments, and per-record
review are documented in this repo as a series of `/tmp/empirical/`
artefacts (not committed) and in the conversation history.

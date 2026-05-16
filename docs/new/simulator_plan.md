epeat Array Simulator — Implementation Plan

A developer-facing plan for building the simulator that produces labelled synthetic tandem repeat arrays for testing the line-width characterisation pipeline. This document covers what to build; the structural taxonomy and pipeline spec define what each construct means biologically.

## 1. Companion documents

These are the authoritative references; this plan supplements rather than replaces them.

- `tr_structural_taxonomy_v1.md` — structural categories, grammar (§2), simulator YAML schema (§5.1), truth file schema (§5.2), wobble realisation (§5.3), test matrix (§5.4).
- `hordetect_image_spec_v2.md` — the detector that consumes simulator output. The simulator must produce ground truth compatible with the detector's property output (§4) and class set (§3).
- `simulator_schema.json` — formal JSON Schema for config validation.

## 2. Purpose and scope

**Purpose.** Generate single-array FASTA test cases with known structural parameters, for evaluating a tandem-repeat line-width classifier at both class level (does the classifier assign `HOR` to a HOR?) and property level (does it recover `k=12` from a `k=12` ground truth?).

**In scope.**
- Generation of arrays from one of the 15+ structural variants in the taxonomy.
- Ambient point mutations and indels.
- Continuous wobble (periodic or aperiodic).
- Discrete phase shifts, non-tandem insertions, hybrid monomers, local inversions, duplications, deletions.
- Truth file matching the detector's property schema.
- Period-candidate file simulating upstream period-generator output.

**Out of scope.**
- Sequence evolution simulation (no phylogeny, no realistic substitution models — uniform random mutations are sufficient).
- Population genetics.
- Discovery of tandem repeats in genomes — input is a configuration, not a genome.
- Variable-length arrays per config (one array per config in MVP; batch mode handles many configs).

## 3. Inputs and outputs

### 3.1 Input

One YAML config file per array, validated against `simulator_schema.json`. Structure described in taxonomy §5.1; formal schema in companion JSON Schema file.

### 3.2 Outputs

Per config (with `--out PREFIX`):

| File | Purpose |
|---|---|
| `PREFIX.fa` | Generated array as a single FASTA record. |
| `PREFIX.truth.tsv` | One-row truth file matching the detector property schema, plus `structural_expression` and `events_json`. See taxonomy §5.2. |
| `PREFIX.periods.tsv` | Period candidates as the detector would receive from upstream. |
| `PREFIX.diagnostics.json` | Optional. Full provenance: coordinate mappings, RNG sub-stream seeds, per-step state. Emitted when `--diagnostics` is set. |

### 3.3 Period candidate file

Simulates what an upstream period generator (TideHunter, TideCluster, NTRprism, etc.) would produce. Schema:

```
array_id    period_bp    period_score    source
```

For each generated array, emit:

- True `base_width` with high score (e.g. 0.9).
- True HOR length if `k > 1`, with high score.
- 2–3 distractor periods (a near-miss within ±5 bp of the true period, a harmonic at 2× or 3× the base width, and optionally one random false positive at low score).

The exact score values don't matter to the detector — it should be robust to candidate ranking — but distinct values prevent the test set from being unrealistically clean.

## 4. Algorithm

### 4.1 Pipeline (ordered)

```
1.  Load and validate config against JSON Schema
2.  Seed RNG; create named sub-streams
3.  Instantiate templates (slot consensuses, monomer templates)
4.  Expand `structure` into a base sequence, tracking coordinates
5.  Apply `modifiers` (wobble → integer edits via residual accumulator)
6.  Apply `post_generation` events (hybrid, inversion, duplication, deletion)
7.  Apply final mutation + indel pass
8.  Emit FASTA, truth.tsv, periods.tsv, optional diagnostics.json
```

Order matters. Mutations and indels go last so that structural truth (which template a base came from, which HOR copy, which slot) is established on a clean sequence. Wobble goes before post-generation events so that hybrid positions stated in logical coordinates still resolve cleanly. Post-generation goes before the noise pass so that inverted segments inherit the same mutation noise as the rest.

### 4.2 RNG and sub-streams

Use a deterministic seeded RNG (the language's default suffices; Python's `random.Random` or numpy `default_rng` are both fine). Create named sub-streams so each generation stage is independently reproducible:

```
top_seed = config.seed
streams = {
    "templates":   derive_seed(top_seed, "templates"),
    "structure":   derive_seed(top_seed, "structure"),
    "wobble":      derive_seed(top_seed, "wobble"),
    "events":      derive_seed(top_seed, "events"),
    "noise":       derive_seed(top_seed, "noise"),
}
```

`derive_seed(s, name) = hash((s, name)) mod 2^32` or equivalent. The benefit: changing the noise level doesn't change which templates were generated.

### 4.3 Template instantiation

Two template types:

**`HOR_slots`.** Generate `k` slot consensus sequences:
1. Draw slot 1 as a random sequence of length `monomer_length_bp` with the requested `gc_content`.
2. Generate slots 2..k by mutating slot 1 at per-base rate `inter_slot_divergence / 2`. This produces between-slot pairwise divergence approximately equal to `inter_slot_divergence` (small under-shoot acceptable; see §11.1 for the calibration note).
3. Cache the resulting list of k sequences under the template ID.

**`monomer`.** As above but `k = 1`.

When the `source` field is `sequence` or `file`, take the provided sequence verbatim and treat it as slot 1; if `k > 1`, derive other slots from it.

Templates are instantiated once. Two structure blocks referencing the same template name receive the *exact same* slot sequences. This is the structural distinction between "one phase-shifted HOR" and "two unrelated HORs with the same k".

### 4.4 Block expansion

Walk through `structure` list. Maintain:
- `output_sequence` (growing string or byte-array; bytes preferred for speed)
- `coordinate_map`: list of (bp_position, block_index, copy_index, slot_index) tuples, one per emitted base

For each block:

**`HOR`.** Append `n_copies` copies of `slot_1 + slot_2 + ... + slot_k` from the named template.

**`SIMPLE_TR`.** Append `n_copies` copies of the monomer sequence. (A `SIMPLE_TR` block referencing an `HOR_slots` template uses slot 1 only.)

**`SHIFT`.** Append `offset_bp` random bases (positive offset) or delete `|offset_bp|` bases from the end of `output_sequence` (negative offset). For positive offset, draw bases from the local composition (last 1 kb of `output_sequence`); if `output_sequence` is empty, fall back to uniform random. Negative offsets must not exhaust the previous block.

**`INSERTION`.** Append `length_bp` bases according to `kind`:
- `random`: Bernoulli(0.5) for each base.
- `AT_rich`: Bernoulli(0.2) G+C.
- `GC_rich`: Bernoulli(0.8) G+C.
- `retro_like`: synthetic LTR–internal–LTR structure. Suggested: two ~200 bp random sequences at the ends (the LTRs), identical to each other; a middle region of length `length_bp - 400` with random composition. Mark coordinates of the LTRs in diagnostics.
- `segdup_like`: copy a random subsequence of the requested length from `output_sequence` (if it exists and is long enough) or fall back to `random`.

Record the insertion's start coordinate and length for the truth file.

### 4.5 Modifier application (wobble)

Apply each modifier in order. For wobble:

1. Determine the affected base range (whole sequence if `target: all`, else from the named block range).
2. Compute the desired smooth shift curve at row resolution:
   - `sinusoidal`: `δ(r) = amplitude_bp · sin(2π · r / period_rows)`
   - `random_walk`: cumulative sum of N(0, σ²) noise with σ chosen to give std ≈ `amplitude_bp` after smoothing with a Gaussian of width `period_rows / 4` (default `period_rows=200` if absent).
3. Walk the affected range row-by-row using the dominant period of the underlying structure (HOR unit length or monomer length — take from the surrounding block). At each row boundary:

   ```
   residual += δ(r) - δ(r-1)
   while residual >= 1:
       insert 1 base near this row boundary
       residual -= 1
   while residual <= -1:
       delete 1 base near this row boundary
       residual += 1
   ```

   Inserted bases drawn from local composition (the surrounding 20–50 bp). Deletions remove the base immediately at the row boundary.

4. Update `coordinate_map` to reflect inserted/deleted positions.

Record the realised wobble amplitude (std of the realised `δ` series) and periodicity in the truth file.

### 4.6 Post-generation events

Apply each event in order. Coordinates in events refer to **logical positions** (e.g. `at_copy: 27` means HOR copy 27 of the surrounding block). Translate to realised positions via `coordinate_map`.

**`HYBRID`.** Locate the bp range corresponding to slot `slot` in copy `at_copy`. Replace this range with a chimera: first half from `source_slots[0]`, second half from `source_slots[1]`. The two halves should have lengths that sum to the original slot length (default: 50/50 split, rounded).

**`INVERSION`.** Locate the bp range corresponding to `length_copies` consecutive HOR units starting at `start_copy`. Reverse-complement that range in place.

**`DUPLICATION`.** Locate the bp range; insert a copy immediately after. Update downstream coordinates.

**`DELETION`.** Locate the bp range; remove it. Update downstream coordinates.

After each event, the `coordinate_map` must be kept consistent so subsequent events resolve correctly and the truth file reports correct realised positions.

### 4.7 Final mutation + indel pass

Single linear walk through `output_sequence`. For each base:
- With probability `mutation_rate`, substitute with a different base (uniform over the other three).
- With probability `indel_rate / 2`, insert a random base before this position.
- With probability `indel_rate / 2`, delete this position.

This step does not update `coordinate_map` for the truth file (mutation positions are not reported); but indel positions do shift downstream realised coordinates, so the truth-file generator must walk through the original positions and apply the same offset adjustments before writing.

A cleaner approach: do this pass *while writing the FASTA*, and write the truth file from the pre-noise sequence with coordinate offsets computed by the indel positions. This avoids needing a full coordinate map post-noise.

### 4.8 Coordinate bookkeeping

The single most error-prone area. Two coordinate systems:

| System | Used for | When valid |
|---|---|---|
| **Logical** | YAML config inputs (`at_copy`, `start_copy`, etc.) | Before any indels/wobble/noise |
| **Realised** | Final FASTA positions, `truth.tsv` coordinate fields | After all steps |

Translation happens at truth-file generation time. The simulator maintains the logical-to-realised mapping incrementally:
- Block expansion: trivial, no mapping needed.
- Wobble: maintain a per-row offset (sum of insertions − deletions up to this row).
- Post-generation events: each event updates a running offset for positions after the event.
- Final noise: positions shift by the cumulative indel net at each downstream position.

For the truth file, every reported position is realised. Diagnostics may include both.

## 5. Output formats

### 5.1 FASTA

Standard FASTA. One record per file. Header: `>{array_id}` where `array_id` is derived from the output prefix or specified explicitly in `global.array_id`. Line width: 80 characters.

### 5.2 truth.tsv

One row. Columns and schema in taxonomy §5.2. Critical fields:

- `structural_expression`: a string in the §2 grammar, reconstructed by walking the original `structure` list (e.g. `H([M_1..M_12],100,div=0.15)|shift(85)|H([M_1..M_12],100,div=0.15)`). Whitespace and newlines collapsed to single spaces.
- `events_json`: JSON array of dicts. One entry per insertion, hybrid, inversion, duplication, deletion. Schema:
  ```json
  {"type":"INSERTION","start_bp":INT,"length_bp":INT,"kind":STR}
  {"type":"HYBRID","at_bp":INT,"copy":INT,"slot":INT,"source_slots":[INT,INT]}
  {"type":"INVERSION","start_bp":INT,"length_bp":INT}
  {"type":"DUPLICATION","start_bp":INT,"length_bp":INT}
  {"type":"DELETION","start_bp":INT,"length_bp":INT}
  ```
- For arrays with no HOR (`k=1`), `hor_k` and `hor_length_bp` are `NA`.
- `n_segments`: number of segments expected after detector-side segmentation. Equals `1 + n_phase_shifts` when shifts produce real segment boundaries.

### 5.3 periods.tsv

Tab-separated, headered:

```
array_id    period_bp    period_score    source
array1      171          0.94            true_base
array1      2052         0.88            true_hor_unit
array1      173          0.71            near_miss
array1      342          0.65            harmonic
array1      895          0.42            false_positive
```

The `source` field is for documentation; the detector should not depend on it.

### 5.4 diagnostics.json

Optional, gated by `--diagnostics`. Structure:

```json
{
  "config": {...},                        // echoed config
  "rng_seeds": {...},                     // per-stream seeds
  "templates": {                          // realised template sequences
    "alpha_A": {"slots": ["ACGT...", ...], "realised_inter_slot_divergence": 0.147}
  },
  "blocks": [                             // per-block provenance
    {"index": 0, "type": "HOR", "start_bp_logical": 0, "start_bp_realised": 0, "end_bp_logical": 205200, "end_bp_realised": 205211}
  ],
  "wobble_realised": {                    // post-realisation wobble stats
    "amplitude_bp_std": 1.43,
    "periodicity_bp": 85510,
    "n_insertions": 47,
    "n_deletions": 39
  },
  "events": [...],                        // applied events
  "noise": {"n_substitutions": 8203, "n_insertions": 412, "n_deletions": 401}
}
```

## 6. CLI

```
simulator simulate --config CONFIG.yaml [--out PREFIX] [--seed SEED]
                   [--diagnostics] [--validate-only]

simulator validate --config CONFIG.yaml
    # Schema-validate without generating.

simulator batch --config-dir DIR --out-dir DIR [--seed-offset INT]
    # Run every *.yaml in DIR. Each gets prefix from filename. Seeds offset
    # per config by (filename hash + seed-offset) to ensure reproducible
    # but distinct runs across the matrix.

simulator schema --print
    # Dump the JSON Schema to stdout for tooling integration.
```

Conventions: CLI flags override config fields where they overlap. `--seed` overrides `seed` in YAML. `--out PREFIX` overrides `global.output`.

## 7. Suggested module structure (Python reference)

Python recommended for the simulator: faster iteration, ample bioinformatics libraries, workload is light (single megabase-class arrays per run). The detector remains Rust (per spec §14).

```
simulator/
  __init__.py
  cli.py              # argparse / click entry points
  config.py           # YAML loading + jsonschema validation
  rng.py              # seeded RNG + sub-streams
  templates.py        # HOR_slots and monomer template instantiation
  blocks.py           # structure-list expansion
  wobble.py           # wobble curves + residual-accumulator integer edits
  events.py           # hybrid, inversion, duplication, deletion
  coords.py           # logical ↔ realised coordinate maps
  noise.py            # final mutation + indel pass
  truth.py            # truth.tsv + structural_expression + events_json
  periods.py          # period candidate generation
  io_fasta.py         # FASTA writing
schema/
  simulator.schema.json
tests/
  test_validate.py
  test_templates.py
  test_simple_tr.py
  test_hor_clean.py
  test_hor_divergence_sweep.py
  test_wobble_sinusoidal.py
  test_wobble_random_walk.py
  test_phase_shift.py
  test_insertion.py
  test_inversion.py
  test_hybrid.py
  test_truth_grammar.py
  test_round_trip.py
configs/
  T01_simple_tr.yaml
  ...
  T18_at_rich.yaml
```

## 8. Testing strategy

Test pyramid:

**Unit tests.** Per module. Most important:
- `templates`: verify `inter_slot_divergence` calibration on multiple `k` values.
- `wobble`: verify std of realised `δ` series approximates requested amplitude (within 10%); verify periodicity recovery via FFT.
- `coords`: verify `logical→realised→logical` round-trip is identity on representative configs.

**Integration tests.** Run each config in `configs/T*.yaml`, verify truth file populates correctly. Spot-check FASTA length matches expectation (block sums + insertions − deletions).

**Property-level tests.** For each test config, recompute properties from the FASTA + truth (independent of the simulator's internal accounting) and verify they match the truth file. This catches coordinate-bookkeeping bugs.

**Determinism tests.** Same seed + same config → byte-identical FASTA. Different seeds → different FASTA but same property values within tolerance.

**Visual sanity check** (manual). Run a clean HOR config, plot the array as a line-width raster at the true `base_width`. Vertical stripes should be visible. This is the integration check that the simulator is producing structurally valid output.

## 9. Implementation milestones

Staged delivery; each milestone produces a useful working tool.

| # | Milestone | Deliverable |
|---|---|---|
| M1 | Schema + validator | Reject invalid configs cleanly; `validate` subcommand works |
| M2 | Templates + HOR + SIMPLE_TR blocks | Generate clean FASTA without noise/wobble; truth files for these configs |
| M3 | Noise pass (mutations + indels) | Truth files reflect noise parameters; FASTA shows expected divergence |
| M4 | SHIFT + INSERTION blocks; periods.tsv | All compositional cases except wobble/inversion working |
| M5 | Wobble (residual accumulator, both models) | T03, T04 configs produce correct realised amplitude/periodicity |
| M6 | Post-generation events (HYBRID, INVERSION, DUP, DEL) | All taxonomy variants covered; `events_json` populated |
| M7 | Diagnostics output + batch mode + test matrix | Full T01–T18 generated reproducibly; ready for detector evaluation |

Each milestone closes a set of test configs from §5.4 of the taxonomy.

## 10. Recommended dependencies (Python)

```
pyyaml >= 6.0          # YAML parsing
jsonschema >= 4.0      # config validation
numpy >= 1.24          # wobble curves, FFT
biopython >= 1.80      # FASTA writing (or roll your own — it's simple)
click >= 8.0           # CLI (argparse fine too)
pytest >= 7.0          # tests
```

If the developer prefers Rust to keep one language across the project: `serde_yaml`, `serde_json`, `jsonschema` (Rust crate), `rand`, `rand_distr`, `clap`, plus the FASTA crate of choice. Either is fine; the workload doesn't benefit from Rust's speed.

## 11. Open implementation choices

Items the developer should decide pragmatically. None are blockers.

### 11.1 Inter-slot divergence calibration

Mutating slot 1 at rate `d/2` to derive slot 2 produces between-slot identity ≈ `1 − d` only asymptotically. At small monomer lengths the realised divergence will deviate. Two options:

(a) Accept the small bias; report realised divergence in diagnostics.
(b) Mutate at slightly higher rate, calibrated empirically per monomer length.

Recommend (a) for MVP. The realised value is what matters for evaluation; the requested value is a target, not a guarantee.

### 11.2 Insertion of bases drawn from "local composition"

For SHIFT blocks and wobble insertions, the spec says "bases drawn from local composition". A simple implementation: sample uniformly from the last N bp of `output_sequence` (N=50). For an empty prefix, fall back to uniform random. The developer can refine if this produces detectable artefacts.

### 11.3 Coordinate update efficiency

Naive coordinate updates after every event are O(L) per event. For large arrays with many events this is acceptable (typical: <100 events, ~10 Mb array → <1 GB-op). If profiling shows it dominates, switch to interval-tree or rope representation. MVP: simple list-of-tuples is fine.

### 11.4 Multi-array configs

YAML can describe one array per file. To run a matrix of tests, use `batch` mode over a directory of YAMLs. An alternative — a top-level `arrays:` list in one YAML — is convenient but complicates the schema. Recommend the batch-mode approach for MVP.

### 11.5 Negative offsets in SHIFT

A negative `offset_bp` deletes bases. The behaviour at block boundaries (does it delete from the preceding HOR's last copy, fully consuming it if `|offset| > monomer_length`?) needs an explicit policy. Recommend: refuse `|offset_bp| > monomer_length / 2`; the use case for large negative shifts is undefined.

### 11.6 Reverse-complement orientation of inverted segments

Inversion reverse-complements DNA: complement each base AND reverse the order. Both steps required. Trivial to get wrong if one is forgotten — write a unit test that inverts `ACGTACGT` and checks the result is `ACGTACGT` (palindromic) or for `AAAA` produces `TTTT`.

## 12. Acceptance criteria for handoff

The simulator is considered ready when:

1. All 18 test configs in taxonomy §5.4 generate without error.
2. `truth.tsv` for each config matches expected ground truth (manually verified for at least T01, T05, T10 — simple TR, clean HOR, phase-shifted HOR).
3. Round-trip determinism: same seed → byte-identical FASTA.
4. Schema validation rejects malformed configs with informative error messages.
5. Documentation covers config writing (a "How to write a test case" appendix is welcome).

Once the simulator is at M7, the detector evaluation harness can run as a separate development track using the generated test set.

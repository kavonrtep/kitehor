tructural Taxonomy of Tandem Repeat Arrays

A purely structural catalogue of tandem-repeat-array architectures, with formal notation, simulator parameters, and line-width detection feasibility per category. Origin (centromere, B-chromosome, holocentromere, intronic VNTR, telomere…) is deliberately excluded — the same structures occur across many biological contexts.

## 1. Scope and assumptions

- **Input**: tandem-repeat arrays from TideCluster. Pre-annotated, oriented canonically, sufficiently long (typically ≥ tens of copies). No flanking non-tandem DNA. Wide range of monomer lengths (single-bp microsatellites to multi-kb megasatellites are out of the immediate target window; the typical range is ~50 bp to a few kb).
- **Out of scope**: array discovery, flank trimming, orientation search, telomere-specific G-skew handling. These are upstream concerns.
- **Goal**: extract structural properties + assign a structural class.

## 2. Notation

A compact grammar for describing arrays. Goals: human-readable, machine-parseable, sufficient to reconstruct any architecture in §3.

### 2.1 Atoms

- `M(L)` — a monomer template of length `L`. Generated randomly or supplied.
- `M_i` — the *i*-th slot of a HOR family (1-indexed). Each slot is itself a template; slots within one family are related variants (see §3.B).

### 2.2 Primitive structures

- `T(M, n)` — **simple tandem repeat**: `n` copies of `M`.
- `H([M_1, ..., M_k], n)` — **higher-order repeat**: `k`-monomer unit `(M_1, M_2, ..., M_k)` repeated `n` times. Slots `M_1, ..., M_k` are typically generated as related variants of a base monomer with controlled divergence `d` (see §2.3).

### 2.3 Generation parameters (modifiers on primitives)

Modifiers are written as keyword arguments inside the constructor or as suffix dots:

- `div=d` — mean pairwise divergence between HOR slots (intra-unit).
- `.mut(p)` — apply per-base point mutation rate `p` to all sequence positions in the structure.
- `.indel(p)` — per-base indel rate `p`.
- `.wobble(A, T?)` — apply continuous period drift of amplitude `A` (bp). Optional `T` = period of the wobble itself in rows; if absent, drift is a smoothed random walk.

### 2.4 Composition operators

Operators combine and disrupt structures.

- `X + Y` — concatenation: `Y` immediately follows `X`.
- `X | shift(δ) | Y` — concatenation with a discrete phase shift of `δ` bp (positive = `Y` starts `δ` bp later in its cycle than expected; negative = earlier).
- `INS(L, kind)` — insertion of `L` bp non-tandem sequence. `kind ∈ {random, AT_rich, GC_rich, retro_like, segdup_like}`.
- `INV(X)` — reverse complement of `X` as a block.
- `DEL(X, start, length)` — delete a sub-interval from `X`.
- `DUP(X, start, length)` — duplicate a sub-interval inside `X`.
- `HYB(M_i, M_j)` — a single hybrid (chimeric) monomer combining the 5′ half of slot `i` and the 3′ half of slot `j`. Used inside an `H(...)` to model single-copy chimeric variants.

### 2.5 Worked examples

```
# 1000-copy simple TR of 170 bp monomer with 2% mutation
T(M(170), 1000).mut(0.02)

# Canonical 12-slot HOR, 200 copies, 15% inter-slot divergence
H([M_1..M_12], 200, div=0.15).mut(0.02)

# Phase-shifted HOR: 100 HOR units, +85 bp phase jump, 100 more
H([M_1..M_12], 100, div=0.15) | shift(85) | H([M_1..M_12], 100, div=0.15)

# Wobbled simple TR
T(M(170), 1000).mut(0.02).wobble(2.0, T=500)

# HOR with a TE-like insertion in the middle
H([M_1..M_12], 50) + INS(5000, retro_like) + H([M_1..M_12], 50)

# Nested HOR: 3-slot inner unit, 4-slot outer
H_outer([H_inner([M_1, M_2, M_3], 1),
         H_inner([M_4, M_5, M_6], 1),
         H_inner([M_7, M_8, M_9], 1),
         H_inner([M_10, M_11, M_12], 1)], 100)

# Mixed families
T(M_A(170), 500) + T(M_B(220), 500)

# Local inversion of 10 HOR units inside a longer array
H([M_1..M_12], 100) + INV(H([M_1..M_12], 10)) + H([M_1..M_12], 100)

# HOR with a single hybrid monomer at position 547
H([M_1..M_12], 100).insert_hybrid(at=547, HYB(M_3, M_4))

# Two layered strata of the same period with different consensus
T(M_old(170), 300).mut(0.05) + T(M_young(170), 1000).mut(0.01)
```

## 3. Structural catalogue

For each category:
- **Definition** — structural only.
- **Schematic** — ASCII.
- **Syntax** — in §2 notation.
- **Simulator parameters** — what must be exposed.
- **Detection feasibility** — line-width method only. Rated:
  - **HIGH** — line-width features give a clean signature.
  - **MEDIUM** — detectable but with caveats or auxiliary computation.
  - **LOW** — line-width alone is weak; needs ancillary signal.
  - **FUNDAMENTAL** — provably indistinguishable from another category at the sequence-statistics level; needs external information (or is by definition the same thing).
- **Confusable with** — the most likely classification error and its disambiguator.

### A. Primitive structures

#### A1. Simple tandem repeat — `T`

**Definition.** A monomer repeated *n* times with mutation/indel noise. No internal periodicity at any non-trivial sub-multiple of the monomer.

```
M M M M M M M M ...
```

**Syntax.** `T(M(L), n).mut(p).indel(q)`

**Simulator parameters.** `L` (monomer length), `n` (copies), `p` (mutation), `q` (indel), monomer sequence (random or supplied), GC content if random.

**Detection feasibility — HIGH.**
- Best width *w* = *L*. Column conservation maximal. R(1) high. No preferred non-trivial *k* > 1.
- This is the calibration case.

**Confusable with.** Degenerated HOR at width `kL` (§D1, FUNDAMENTAL case). Distinguished only by ghost stripes at a sub-period — and not always.

#### A2. Higher-order repeat — `H`

**Definition.** A `k`-slot unit `(M_1, ..., M_k)` repeated `n` times, where slots are related variants with controlled divergence `d`. Same-slot copies across units are more similar to each other than different-slot copies within the same unit.

```
[M_1 M_2 ... M_k] [M_1 M_2 ... M_k] [M_1 M_2 ... M_k] ...
```

**Syntax.** `H([M_1..M_k], n, div=d).mut(p).indel(q)`

**Simulator parameters.** Monomer length `L`, multiplicity `k`, inter-slot divergence `d`, number of HOR units `n`, point mutation `p`, indel rate `q`.

**Detection feasibility — HIGH in moderate divergence regime (`d` ≈ 0.05 to 0.40).**
- Column IC at width *L*: moderate (columns are multimodal across `k` slot consensuses).
- Column IC at width `kL`: maximal.
- R(*L*): low to moderate; R(`kL`): strong. The ratio R(`kL`)/R(*L*) is the canonical diagnostic.
- Phase separation large.

**Confusable with.**
- *Simple TR at `kL`* when `d` is too large (regime collapse, see D1). Disambiguator: residual column structure at *L*.
- *Simple TR at `L`* when `d` is too small. Disambiguator: phase separation drops to zero.

#### A3. Nested HOR — `H` inside `H`

**Definition.** An HOR whose unit is itself a small HOR. Total multiplicity = `k_outer · k_inner`. Distinct from A2 with `k = k_outer · k_inner` because the slot relationships are hierarchical, not flat.

```
[(M_1 M_2 M_3)(M_4 M_5 M_6)(M_7 M_8 M_9)(M_10 M_11 M_12)] × n
   inner unit 1   inner unit 2  inner unit 3   inner unit 4
```

**Syntax.** `H([H([M_1..M_3], 1), H([M_4..M_6], 1), ...], n)`

**Simulator parameters.** `k_inner`, `k_outer`, inter-slot divergence at each level (potentially different — inner slots may diverge more than outer slot-group consensuses, or vice versa), monomer length, copies.

**Detection feasibility — HIGH.**
- R(k) shows commensurate peaks at `L`, `k_inner · L`, `k_outer · k_inner · L`.
- Column IC has a hierarchical signature: highest at outermost period.
- The hierarchical R(k) pattern is the diagnostic. A flat A2 HOR with `k = 12` shows one strong peak; a nested 3 × 4 shows peaks at 3 and 12 (and possibly 4 and 6 depending on slot relationships).

**Confusable with.**
- *Flat HOR with same total k* — disambiguator: R(k) at intermediate divisors (3, 4, 6) reveals hierarchical structure.

### B. Local modifications

These are applied to a primitive structure and produce a derived structure. They can compose.

#### B1. Mutations and indels — `.mut(p).indel(q)`

**Definition.** Per-base random substitutions and length-changing edits.

**Detection feasibility — N/A (baseline).** These are the ambient noise level. Every other detection feasibility rating assumes some baseline `p` and `q` are present. Effects:
- Increasing `p` raises column entropy uniformly across widths.
- Increasing `q` introduces drift-like artifacts at the indel scale.

The simulator must always allow `p` and `q` so calibration curves can be built across noise levels.

#### B2. Continuous wobble — `.wobble(A, T?)`

**Definition.** The period drifts smoothly along the array, either as a smoothed random walk (no fixed period) or sinusoidally with period `T`. Monomer length is effectively non-integer.

```
M M' M  M' M   M'  M  ...
^ same monomer, but separations drift smoothly
```

**Syntax.** `T(M, n).wobble(A=2.0, T=500)` or `H(...).wobble(A=1.5)`

**Simulator parameters.** Amplitude `A` (bp), optional period `T` (rows), drift model (random walk or sinusoid). Drift is realised by sampling a per-row shift `δ(r)` and inserting/deleting bases accordingly.

**Detection feasibility — HIGH.**
- The per-row best-shift signal `best_shift(r)` is non-zero and smoothly varying. Mean is zero (after width refinement) but **standard deviation = wobble amplitude**.
- 1D FFT of detrended `best_shift(r)` reveals the wobble period if present.
- Column IC moderately degraded compared to wobble-free case.

**Confusable with.**
- *Slightly wrong width* — disambiguator: drift from wrong width is constant; wobble varies along the array.
- *Phase shift* — disambiguator: wobble is smooth; phase shift is a step.

#### B3. Discrete phase shift — `| shift(δ) |`

**Definition.** At one or more positions in the array, the phase of the underlying repeat jumps discontinuously by `δ` bp. The structure on both sides of the shift is otherwise the same.

```
[M_1 M_2 M_3 M_1 M_2 M_3] | shift(+L) | [M_2 M_3 M_1 M_2 M_3 M_1]
                       phase advanced by one slot
```

**Syntax.** `X | shift(δ) | Y`

**Simulator parameters.** Positions of shifts (list), magnitudes `δ` of shifts.

**Detection feasibility — HIGH.**
- `best_shift(r)` shows a step function. Threshold on `|Δ best_shift|` detects the breakpoint cleanly.
- For an HOR with phase shift of one slot (`δ = L`), R(k) at the HOR width is degraded by a sawtooth pattern; segmenting at the breakpoint and recomputing R(k) per segment restores the clean signal.

**Confusable with.**
- *Wobble* — disambiguator: step vs smooth.
- *Local inversion boundary* — disambiguator: inversion produces breakpoints at both ends with reverse-complement match between them; phase shift does not.

#### B4. Hybrid monomer — `HYB`

**Definition.** Within an HOR, a small number of unit copies contain a chimeric monomer in one slot — for example slot 3 of copy 547 is a 5′-half of slot 3 + 3′-half of slot 4 (or any analogous fusion).

```
copy 546: M_1 M_2 M_3 M_4 ...
copy 547: M_1 M_2 [HYB(3,4)] M_4 ...   <-- one hybrid
copy 548: M_1 M_2 M_3 M_4 ...
```

**Syntax.** `H([M_1..M_k], n).insert_hybrid(at=547, HYB(M_3, M_4))`

**Simulator parameters.** Number of hybrids, positions (which HOR copy and which slot), source slots `i, j` of the chimera.

**Detection feasibility — MEDIUM.**
- Aggregate R(k) and column IC are barely affected by a few hybrids.
- Detection requires per-row outlier analysis: row similarity to expected slot-consensus drops below typical for that slot. Equivalently, detect rows whose k-mer embedding lies between two slot-consensus embeddings.
- For a few hybrids in a long array, this is local detection, not array-level classification. The array still classifies as HOR.

**Confusable with.**
- *Random row in degraded region* — disambiguator: a hybrid row sits between two slot consensuses, not in random space.

#### B5. Local internal duplication / deletion — `DUP, DEL`

**Definition.** A contiguous span of monomers (or HOR units) is duplicated or deleted internally.

**Syntax.** `H([M_1..M_k], n).dup(start_unit=50, length=5)` or `.del(start_unit=50, length=5)`

**Simulator parameters.** Position (unit index), length (units), type (dup/del).

**Detection feasibility — MEDIUM.**
- Sliding-window HOR-unit count changes locally.
- Phase signal may show a transient at the boundary, especially for deletions that are not multiples of the HOR width.
- A deletion exactly equal to one HOR unit is undetectable in column IC and R(k); it shows only as a copy-number anomaly relative to expected coordinates.

**Confusable with.**
- Nothing structural; this is a quantitative anomaly within an otherwise stable architecture. Often reported as a "structural variant within HOR" rather than a class change.

#### B6. Local inversion — `INV`

**Definition.** A contiguous segment of the array is reverse-complemented relative to its surroundings.

```
H(...,100) + INV(H(...,10)) + H(...,100)
       forward      reverse        forward
```

**Syntax.** `X + INV(Y) + Z`

**Simulator parameters.** Position, length of inverted segment.

**Detection feasibility — MEDIUM (LOW without strand awareness).**
- Without strand-aware comparison, the inverted segment appears as a discontinuity at both boundaries: column IC drops sharply at the entry and exit. The shift signal shows two breakpoints.
- With strand-aware comparison (compare each row against both forward and reverse-complement of the consensus), the inverted segment is identified as a reverse-complement-matching block.
- An inversion of one HOR unit is harder to detect than a longer one because the boundary signal dominates.
- **Deferred to v2 of the pipeline.** MVP detects inversions as two unexplained breakpoints; strand-aware testing is a planned extension once the core method is validated.

**Confusable with.**
- *Two phase shifts close together* — disambiguator: reverse-complement match between the two breakpoints.
- *TE insertion* — disambiguator: inversion preserves the satellite consensus (in RC); TE insertion has unrelated sequence.

### C. Compositional architectures (multiple atomic blocks)

#### C1. Non-tandem insertion — `INS`

**Definition.** A foreign sequence (random, retroelement-like, segmental duplication, etc.) interrupts the array.

**Syntax.** `X + INS(L, kind) + Y` where `X` and `Y` are typically the same architecture.

**Simulator parameters.** Position, length `L`, kind (controls composition — uniform, AT-rich, GC-rich, or a synthetic LTR-internal-LTR pattern for `retro_like`).

**Detection feasibility — HIGH.**
- Column IC at the satellite width collapses inside the insertion.
- R(k) at the satellite width breaks; reappears after the insertion.
- Sliding-window analysis flags the insertion as a gap. The shift signal does not propagate through the insertion.
- The insertion itself may have its own internal periodicity (retroelement LTR structure) detectable as a secondary best-*w*.

**Confusable with.**
- Very long deletion plus very different fill — disambiguator: actual insertion has detectable internal structure or random composition; pure deletion would not show a long uncharacterised block.

#### C2. Multiple coexisting periods — concatenated `H` blocks with different `k`

**Definition.** Two or more sub-regions of one array each have a clean HOR architecture, but with different multiplicities or different monomer families.

```
H([M_1..M_12], 100) + H([N_1..N_8], 100)
   region 1, k=12       region 2, k=8
```

**Syntax.** `H_A + H_B` (concatenation of distinct HOR blocks).

**Simulator parameters.** Number of blocks, parameters per block, transition position(s).

**Detection feasibility — HIGH.**
- Per-window best-*w* and best-*k* differ between blocks.
- Whole-array R(k) shows competing peaks; sliding-window R(k) resolves which is which.
- The transition itself produces a breakpoint in the shift signal.

**Confusable with.**
- *Nested HOR* — disambiguator: nested HOR's commensurate peaks coexist throughout; coexisting periods occupy different regions.

#### C3. Mixed unrelated families

**Definition.** Same as C2 but the monomers are entirely unrelated rather than being different HOR multiplicities of related material.

**Syntax.** `T(M_A, 500) + T(M_B, 500)` or any combination.

**Simulator parameters.** Per-block monomer source.

**Detection feasibility — HIGH.** Same as C2; even cleaner because the monomer consensuses share no information.

**Confusable with.** *Stratification* (C4) when families happen to have similar lengths — disambiguator: family monomers are unrelated by sequence; strata are sequence-related variants.

#### C4. Stratification (same period, divergent consensus)

**Definition.** Two or more sub-regions share the same period (e.g. monomer length 170 bp throughout) but the monomer consensus differs between regions. The regions are *related* in sequence (typically homologous) but have undergone independent drift, so columns conserve within a region but differ between regions.

```
T(M_old(170), 300).mut(0.05) + T(M_young(170), 1000).mut(0.01)
   diverged stratum             active stratum
```

**Syntax.** `T(M_A, n_1) + T(M_B, n_2)` where `M_A` and `M_B` are related variants of a shared ancestor.

**Simulator parameters.** Number of strata, per-stratum consensus divergence from a shared ancestor, per-stratum mutation rate.

**Detection feasibility — MEDIUM.**
- Per-window column IC is high within each stratum but drops at the boundary.
- Per-window monomer consensus changes; a global consensus is bimodal at many columns.
- Whole-array column IC underestimates the true within-stratum conservation.
- Detection requires sliding-window consensus comparison, not just sliding-window R(k).

**Confusable with.** *Mixed unrelated families* (C3) when divergence is high — disambiguator: stratum consensuses share appreciable similarity; unrelated families do not.

### D. Boundary regimes

#### D1. Degenerated HOR

**Definition.** A structure that was generated as `H([M_1..M_k], n, div=d)` with `d` so large that no statistical periodicity remains at width `L`. Mathematically, the smallest valid period is `kL` and the architecture is indistinguishable from `T(M(kL), n)` where the monomer is the entire former HOR unit.

```
[A_1 A_2 A_3 A_4][A_1' A_2' A_3' A_4'][A_1'' A_2'' A_3'' A_4'']
                  effectively becomes a simple TR at width 4L
```

**Syntax.** `H([M_1..M_k], n, div=d)` with `d` very large (typically > 0.5).

**Simulator parameters.** Same as A2 but the divergence is in the collapse range.

**Detection feasibility — FUNDAMENTAL LIMIT.**
- At width `L`: column IC collapses, R(*L*) low.
- At width `kL`: column IC high, R(1) high.
- The mathematically correct call is `simple_TR at width kL`. Any tool that calls this "HOR" is reading structure that is not statistically there.
- The line-width method correctly falls through to `simple_TR at kL` when the regime collapses (per spec §8).
- Faint residual stripes at `L` exist if `d` is at the boundary, but a clean call requires external evidence (e.g. comparison with sister haplotypes that retain the HOR signal).

**Confusable with.** *Pure simple TR with monomer = kL* (A1 at large `L`). Disambiguator: ghost stripes at sub-period `L` — but not always present. This is the irreducible ambiguity the field has acknowledged for decades.

#### D2. Decay zones at array boundaries

**Out of scope.** Assumed trimmed by TideCluster. If residual non-tandem flanks are present, they manifest as low column conservation at all widths over the first/last `O(monomer_length)` of the array. The pipeline should not refuse such cases; it should report the property and let the user decide whether to re-trim.

## 4. Cross-reference table

The structural taxonomy enumerates many architectures, but the pipeline output uses a small set of **class labels** (`simple_TR`, `HOR`, `irregular_HOR`, `mixed`, `ambiguous`) plus a **property vector** that carries the structural detail. The table below shows which structural variants change the class and which are reported only as properties.

**Reported as a class change** — A1, A2, A3 (these are the underlying repeat architectures); C2, C3 (different architectures in different blocks → `mixed`); D1 (regime collapse → class fallthrough). 

**Reported only as properties on top of an A1/A2/A3 base class** — B2 (`wobble_amplitude`, `wobble_periodicity`), B3 (`n_phase_shifts`, `phase_shift_positions`, `phase_shift_offsets`), B4 (`hybrid_fraction`), B5 (`local_cnv_events`), B6 (`n_inversions`, `inversion_positions` — deferred), C1 (`n_insertions`, `insertion_positions`, `insertion_lengths`).

**Reported as `irregularity_score` contribution** — C4 stratification (per-window column-IC variance), elevated B5 counts, and any case where the underlying architecture is preserved but realised noisily.

The rationale: a HOR with a phase shift halfway through is still a HOR (same monomer length, same `k`, same slot consensuses on both sides). The shift is a localised event captured by a property. Only when the segments on either side of a shift have *different* underlying architectures does the class change (to `mixed`).

| ID | Category | Notation | Detection | Key signature |
|---|---|---|---|---|
| A1 | Simple TR | `T(M, n)` | HIGH | one peak in R(k); high col IC at `L`; no preferred multiple |
| A2 | HOR | `H([M_1..M_k], n, div=d)` | HIGH (moderate `d`) | R(kL)/R(L) ≫ 1; phase separation |
| A3 | Nested HOR | `H([H([..], 1)..], n)` | HIGH | commensurate peaks at multiple lags |
| B2 | Wobble | `.wobble(A, T?)` | HIGH | std of `best_shift(r)` > 0; FFT peak if periodic |
| B3 | Phase shift | `| shift(δ) |` | HIGH | step in `best_shift(r)` |
| B4 | Hybrid monomer | `HYB(M_i, M_j)` | MEDIUM | per-row outlier between two slot consensuses |
| B5 | Local CNV | `.dup() / .del()` | MEDIUM | window-local copy count anomaly |
| B6 | Local inversion | `INV(X)` | MEDIUM | needs strand-aware compare; two breakpoints with RC match |
| C1 | Non-tandem insertion | `INS(L, kind)` | HIGH | column IC gap; shift signal interrupts |
| C2 | Coexisting periods | `H_A + H_B` (diff k) | HIGH | sliding-window best-w changes |
| C3 | Mixed families | `T_A + T_B` (unrelated) | HIGH | sliding-window best-w changes + unrelated consensuses |
| C4 | Stratification | `T(M_A, .) + T(M_B, .)` related | MEDIUM | sliding-window consensus changes, period stable |
| D1 | Degenerated HOR | `H(..., div=large)` | FUNDAMENTAL | mathematically = simple_TR at `kL` |
| D2 | Boundary decay | (out of scope) | — | upstream concern |

## 5. Simulator design

The simulator interprets a structured config that mirrors the §2 grammar. Three pieces: a YAML input schema, a truth output schema, and a concrete wobble-realisation algorithm.

### 5.1 YAML schema

Use a separate `templates` section. Each template is a named, fully-specified object; structural blocks reference templates by ID. This makes "two HOR blocks share the same slot consensuses" explicit and reproducible. A free-form tag like `monomer_seed` would be ambiguous (same random seed? same base monomer? same already-mutated sequence?) and will eventually cause silent bugs.

```yaml
schema_version: 1
seed: 42

global:
  mutation_rate: 0.02
  indel_rate: 0.005
  output: prefix

templates:
  alpha_A:
    type: HOR_slots
    monomer_length_bp: 171
    k: 12
    source: random            # random | sequence | file
    gc_content: 0.45
    inter_slot_divergence: 0.15

  simple_M170:
    type: monomer
    monomer_length_bp: 170
    source: random
    gc_content: 0.45

structure:
  - type: HOR
    template: alpha_A
    n_copies: 100

  - type: SHIFT
    offset_bp: 85

  - type: HOR
    template: alpha_A            # same slot consensuses → genuinely one phase-shifted HOR
    n_copies: 100

  - type: INSERTION
    length_bp: 5000
    kind: retro_like             # random | AT_rich | GC_rich | retro_like | segdup_like

  - type: HOR
    template: alpha_A
    n_copies: 50

modifiers:
  - target: all                  # all | block_index_N | range
    wobble:
      amplitude_bp: 1.5
      period_rows: 500           # absent or 0 = aperiodic random walk
      model: sinusoidal          # sinusoidal | random_walk
      realisation: integer_edits # see §5.3

post_generation:
  - type: HYBRID
    at_copy: 27
    slot: 4
    source_slots: [4, 5]

  - type: INVERSION
    start_copy: 60
    length_copies: 5
```

#### Field-name conventions

Consistent suffixes prevent parser bugs and make the truth file generation straightforward:

| Concept                       | Field                    |
|---|---|
| schema version                | `schema_version`         |
| output prefix                 | `global.output`          |
| block kind                    | `type`                   |
| monomer length                | `monomer_length_bp`      |
| HOR multiplicity              | `k`                      |
| HOR copies                    | `n_copies`               |
| inter-slot divergence         | `inter_slot_divergence`  |
| point mutation rate           | `mutation_rate`          |
| indel rate                    | `indel_rate`             |
| wobble amplitude              | `amplitude_bp`           |
| wobble period (rows)          | `period_rows`            |
| phase shift                   | `offset_bp`              |
| insertion length              | `length_bp`              |
| template reference            | `template`               |

Avoid mixing `block`/`type`, `length`/`length_bp`, `offset`/`offset_bp`, `monomer_seed`/`template`. Pick one and stay with it.

### 5.2 Truth file schema

Split into two files to keep generative truth distinct from derived measurements.

#### `truth.tsv` — true structural parameters only

```
array_id
length_bp
truth_class
base_width_bp
hor_k
hor_length_bp
n_complete_copies
wobble_amplitude_bp
wobble_periodicity_bp
n_phase_shifts
phase_shift_positions
phase_shift_offsets
n_segments
reason
structural_expression           # mandatory; grammar string from §2
schema_version
events_json                     # optional; see below
```

For fields that have no meaningful generative value (e.g. `column_conservation`, `phase_separation`, `confidence`, `irregularity_score` — these are detector measurements, not parameters of the simulator), they belong in detector output only.

#### `structural_expression` is mandatory

Class-level accuracy is too weak. A prediction that gets `class=HOR` right but misses `k=12` or the phase shift would score 100% on class accuracy and be useless. The structural expression makes finer-grained scoring possible. Example:

```
H([M_1..M_12],100,div=0.15) | shift(85) | H([M_1..M_12],100,div=0.15)
```

The test harness can compare predicted vs truth at multiple levels (class, base_width, k, phase shifts, etc.) without needing custom per-test scoring code.

#### `events_json` for structural variants not in the property table

Hybrids, insertions, inversions, and local CNV events are not first-class detector outputs in the MVP property table. Keep the TSV simple by emitting one JSON column with the event list:

```json
[
  {"type":"INSERTION","start_bp":1026000,"length_bp":5000,"kind":"retro_like"},
  {"type":"HYBRID","copy":27,"slot":4,"source_slots":[4,5]},
  {"type":"INVERSION","start_copy":60,"length_copies":5}
]
```

This preserves truth for events that may become first-class properties in v2 without bloating the MVP schema.

### 5.3 Wobble realisation — integer edits with residual accumulator

DNA is discrete. A "0.3 bp shift" is not a real base-level operation. Sub-base interpolation would create an artificial numerical signal rather than realistic sequence. Use integer insertions/deletions spaced along the array to realise non-integer average drift.

```text
For each row r in the array:
    desired_shift[r] ← target wobble curve at r
        sinusoidal:    A · sin(2π·r / T)
        random walk:   smoothed cumulative-sum noise, std ≈ A

    residual += desired_shift[r] - desired_shift[r-1]

    while residual ≥ 1:
        insert 1 base near this row boundary
        residual -= 1

    while residual ≤ -1:
        delete 1 base near this row boundary
        residual += 1
```

Bases inserted/deleted are drawn from the local sequence composition (a base copied from a nearby position is the simplest realistic model; pure random bases would produce a colour-distinct streak in the line-width raster).

This produces non-integer average drift while every realised edit is a real base-level operation. The detector sees what it would see in real degraded arrays: small insertions and deletions that accumulate into smooth shift along the array. Sub-base interpolation is not added to MVP; if ever needed, it would be a diagnostic stress test, not the default.

### 5.4 Test matrix to generate

Minimum coverage:

| Test | Structure |
|---|---|
| T1 | `T(M(170), 1000).mut(0.02)` |
| T2 | `T(M(170), 1000).mut(0.02).indel(0.01)` |
| T3 | `T(M(170), 1000).mut(0.02).wobble(2.0)` |
| T4 | `T(M(170), 1000).mut(0.02).wobble(2.0, T=500)` |
| T5 | `H([M_1..M_12], 200, div=0.15).mut(0.02)` |
| T6 | `H([M_1..M_4], 200, div=0.0).mut(0.02)` (regime A — should call simple_TR) |
| T7 | `H([M_1..M_4], 200, div=0.7).mut(0.02)` (regime C — should call simple_TR at 4L) |
| T8 | Divergence sweep at k=4: `d ∈ {0.0, 0.05, 0.15, 0.35, 0.55, 0.7}` |
| T9 | `H_outer([H_inner([M_1..M_3], 1) × 4], 100)` |
| T10 | `H([M_1..M_12], 100) \| shift(85) \| H([M_1..M_12], 100)` |
| T11 | `H([M_1..M_12], 50) + INS(5000, retro_like) + H([M_1..M_12], 50)` |
| T12 | `H([M_1..M_12], 100) + INV(H([M_1..M_12], 10)) + H([M_1..M_12], 100)` |
| T13 | `H([M_1..M_12], 100) + H([N_1..N_8], 100)` (coexisting k) |
| T14 | `T(M_A(170), 500) + T(M_B(220), 500)` (mixed families) |
| T15 | `T(M_old(170), 300).mut(0.05) + T(M_young(170), 1000).mut(0.01)` (stratification) |
| T16 | `H([M_1..M_12], 200).insert_hybrid(at=100, HYB(M_3, M_4))` |
| T17 | random sequence (negative control) |
| T18 | `T(M(170), 1000).mut(0.02)` with GC=0.2 (composition bias) |

Each runs at multiple noise levels.

## 6. Quick reference: detection by line-width feature

| Feature | Detects |
|---|---|
| Column IC at width *w* | base period validity; A1 vs A2 base-width call; D1 collapse |
| Column IC at width `kw` | A2 vs A1; A3 nested hierarchy |
| R(1) at *w* | A1 (high) vs A2 (lower) |
| R(k) at *w*, k≥2 | A2 multiplicity; A3 hierarchy |
| Phase separation | A2 vs A1 (must be high for A2) |
| Sliding-window column IC | C2, C3, C4 boundaries; C1 insertions; A3 nesting; B6 inversion entry/exit |
| Sliding-window best-*w* | C2, C3 |
| Sliding-window monomer consensus | C4 stratification |
| `best_shift(r)` mean | width refinement |
| `best_shift(r)` std | B2 wobble amplitude |
| FFT of `best_shift(r)` | B2 wobble period |
| Steps in `best_shift(r)` | B3 phase shifts; B6 inversion boundaries; C1 insertion boundaries |
| Reverse-complement column-IC | B6 inversion confirmation |
| Per-row similarity to slot consensus | B4 hybrid; B6 inversion content |
| Local block features | B5 CNV |

## 7. Things this taxonomy deliberately omits

- **Origin / biological context.** Same structure in a centromere, a B-chromosome, or an intronic VNTR has the same fingerprint.
- **Monomer length per se.** A 5 bp simple TR and a 500 bp simple TR are the same class A1. The pipeline reports the length as a property; classification doesn't depend on it.
- **Telomere-specific signatures (G-skew).** Out of scope for a generic line-width pipeline; telomeric arrays are class A1 with side information.
- **Functional annotation (CENH3 binding, transcribed rDNA, etc.).** External data, not derivable from sequence structure.
- **Higher-order spatial organisation across chromosomes** (e.g. pan-centromere comparisons). One array at a time.
- **rDNA-specific architecture** (intergenic spacer sub-repeats inside a transcribed unit) is class A3 nested-period structurally; the biological labelling is separate.

## 8. Open questions and decisions

**Resolved:**
1. **Hybrid monomer threshold** — calibrate on simulator output. Report `hybrid_fraction` as a property; the threshold above which the array-level class shifts from `HOR` to `irregular_HOR` is set during calibration.
2. **Strand-aware comparison for inversions** — deferred to v2 of the pipeline. MVP detects inversions as two unexplained breakpoints. Strand-aware row comparison (compare against forward and reverse-complement consensus, take max similarity) added in v2 at ~2× row-comparison cost.
3. **Nested HOR depth** — two levels only. Three-level nesting is not in the test set unless real data demands it.
4. **Phase shift vs wobble vs class** — both are properties, not class labels. Class is determined by the underlying repeat architecture; wobble and phase shifts describe how that architecture is realised along the array. A `HOR` with phase shifts is still a `HOR`. Only when the segments on either side of a shift have *different* underlying architectures (different `base_width` or different `k`) does the class become `mixed`.

**Remaining:**

5. **Stratification combined with phase shift.** If two strata share the same period but their phase differs at the boundary, this combines C4 with B3. The grammar handles it; the classifier needs an explicit check: when a phase shift is detected, compare segment monomer consensuses, not just `base_width` and `k`. If consensuses differ substantially → call `mixed` rather than reporting the shift as a benign property. Calibration on simulator output decides where this threshold sits.

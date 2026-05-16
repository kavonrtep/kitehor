ne-Width Characterization of Tandem Repeat Arrays — v2

A method to characterize tandem-repeat arrays by treating candidate periods as line-width hypotheses, extracting a set of quantitative properties (period structure, HOR multiplicity, drift, wobble, phase shifts, irregularity), and producing both a property vector and a summary class label.

---

## 1. Concept glossary

| Term | Definition |
|---|---|
| **Monomer** | Smallest repeated sequence unit. Its length is the **monomer length** (e.g. 171 bp for human alpha satellite). |
| **Period** | Any length at which the sequence repeats. The monomer length is a period; integer multiples of it are also periods. |
| **HOR (higher-order repeat)** | A repeat unit composed of an ordered series of **related but distinguishable** monomer variants derived from a common ancestor, e.g. `A₁ A₂ A₃ A₄ A₁ A₂ A₃ A₄`. The defining property: copies of the same slot across HOR units (A₁↔A₁) are more similar than monomers in different slots within the same unit (A₁↔A₂). |
| **HOR multiplicity (k)** | Number of monomer slots per HOR unit. Simple tandem repeats have `k = 1`. |
| **HOR unit length** | `monomer_length × k`. |
| **Inter-monomer divergence** | Mean sequence dissimilarity between slots within a HOR unit (A₁ vs A₂). Controls whether the HOR structure is mathematically detectable — see §3.1. |
| **Line width (w)** | Candidate period currently being tested as a wrap width. |
| **Row** | One line of the wrapped sequence, length `w`. |
| **Row autocorrelation R(k)** | Mean similarity between rows `i` and `i+k`. |
| **Phase** | Position within the HOR cycle. Two rows are in-phase if they correspond to the same monomer slot. |
| **Drift** | Constant offset between adjacent rows due to wrong line width. Slope of the diagonal pattern. |
| **Wobble** | Slow continuous variation of the true period along the array, observable as smooth horizontal meandering of vertical edges at the correct width. |
| **Phase shift** | Discrete, localised jump in monomer phase — array structure continues on both sides but with a different offset. |
| **Domain transition** | Change in `k` or in the underlying monomer between segments. |

## 2. Worked example

Input: a 4.2 Mb human alpha satellite array. Period generator returns `{171, 342, 2052}`.

At `w = 171`:
- Each row is one candidate monomer.
- Column conservation high; `R(1) = 0.38`; `R(12) = 0.86`; phase separation `0.41`.
- Vertical edges dominate, wobbling left–right by 1–2 bp across the array (biological wobble at correct width).
- `class_hint = HOR_base_width`, `k = 12`.

At `w = 2052`:
- Each row is one complete HOR unit.
- `R(1) = 0.91`, phase separation `~0`. No internal phase structure at this width.
- `class_hint = HOR_unit_width`.

At `w = 342`:
- Partial structure, edges slant diagonally — wrong width.
- `class_hint = unsupported_width`.

Final call:
- `base_width = 171`, `k = 12`, `hor_length = 2052`, `class = HOR`.
- `wobble_amplitude = 1.4 bp`, `n_phase_shifts = 0`, `confidence = 0.92`.

A second example: same array but with a single phase shift halfway through. At `w = 171`, the shift signal jumps from `~0` to `+85` at row 12,000 and stays there. Both segments still classify as HOR with the same `base_width = 171` and `k = 12`. The array class is `HOR`; phase shift detail is reported as `n_phase_shifts = 1`, `phase_shift_positions = [2,052,000 bp]`, `phase_shift_offsets = [+85]`.

A third example illustrates the regime boundary. Same period candidates, but the array is synthetic with `k = 4` slots that have diverged to ~50% pairwise identity. At `w = 171`: column IC is low (columns mix unrelated bases from the four slot consensuses), but `R(4) >> R(1)`. At `w = 684` (HOR unit length): column IC high, `R(1)` high. Because column IC at base_width fails the HOR threshold, the call falls through to `simple_TR` at width 684 with `reason = "regime C: HOR with k=4 base period has degenerated"`. This is the correct mathematical interpretation: the base period is no longer a valid statistical period of the sequence.

## 3. Problem statement

Given a tandem-repeat array sequence and a list of candidate periods from an upstream generator, extract a property vector describing the array's repeat structure and assign a summary class:

```
simple_TR | HOR | irregular_HOR | mixed | ambiguous
```

The method evaluates each candidate period as a line-width hypothesis, extracts features from the wrapped 2D representation, and combines width-level evidence into array-level properties.

**Phase shifts and wobble are properties, not classes.** A HOR with discrete phase shifts is structurally still a HOR (same `base_width`, same `k`, same slot consensuses on both sides of the shift). A simple TR with continuous wobble is structurally still a simple TR. The class captures the underlying repeat architecture; properties such as `n_phase_shifts`, `phase_shift_positions`, `wobble_amplitude`, and `wobble_periodicity` capture how that architecture is realised along the array.

**Assumptions about input.** Arrays are pre-annotated tandem repeats: oriented canonically, free of non-tandem flanks, and starting on or near a monomer boundary. No offset search is performed. Violations manifest as low column conservation at all widths.

### 3.1 The HOR detection regime

A HOR with multiplicity `k` is mathematically detectable only when the inter-monomer divergence within a HOR unit falls in a specific range. The structure collapses to simple tandem repeat at both extremes. The detector must produce the correct label across the whole spectrum, not force a HOR call wherever phase separation looks suggestive.

**Regime A — zero or near-zero divergence (A₁ ≈ A₂ ≈ ... ≈ A_k).**
All slots are essentially identical copies. The HOR period exists only as a harmonic of the base period. At `base_width`: column conservation high, `R(1) ≈ R(k)`, phase separation low. Correct label: `simple_TR` at `base_width`.

**Regime B — moderate divergence (slots distinguishable but homologous).**
Same-slot copies across HOR units are more similar than different-slot monomers within a unit. At `base_width`: column conservation moderate (columns are multimodal across the `k` slot consensus bases), `R(1) < R(k)`, phase separation high. At `HOR_unit_width`: column conservation also high, `R(1)` high. Correct label: `HOR`.

**Regime C — extreme divergence (slots have become unrelated sequences).**
The base period is no longer a statistical period of the sequence. At `base_width`: column conservation collapses (columns mix unrelated positions), `R(1)` low, `R(k)` high — but the base period is mathematically not valid anymore. At `HOR_unit_width`: column conservation high, `R(1)` high. Correct label: `simple_TR` at `HOR_unit_width`. The HOR has degenerated into a longer simple repeat.

**Implication for detection.** A HOR call requires:
- `base_width` shows adequate column conservation (slots remain recognizable as variants), AND
- `R(k) > R(1)` with adequate phase separation, AND
- `HOR_unit_width` also shows high column conservation (the HOR unit is itself conserved).

If `base_width` column conservation is low but `HOR_unit_width` has strong simple-repeat signal, the correct call is `simple_TR` at the HOR-unit width — the structure has degenerated. The detector must not call HOR purely on the row-autocorrelation signal.

**Inter-monomer identity is itself a useful property** and is reported when the call is HOR. It tells the user where on the spectrum the array sits and is informative about the evolutionary state of the array.

## 4. Properties extracted

The pipeline reports the following properties per array. Class is a derived summary; the property vector is the primary output.

| Property | Description |
|---|---|
| `base_width_bp` | Inferred monomer length. |
| `hor_k` | HOR multiplicity (1 = simple TR). |
| `hor_length_bp` | `base_width × hor_k`. |
| `n_complete_copies` | Estimated number of full HOR units (or monomers if `k=1`). |
| `column_conservation` | Mean information content per column at the chosen width, background-corrected. |
| `phase_separation` | `R(best_k) − background_R` (HOR phase signal strength). |
| `mean_shift_bp` | Global row-to-row shift; should be ~0 after width refinement. |
| `wobble_amplitude_bp` | Std. dev. of the local shift signal around its mean. |
| `wobble_periodicity_bp` | Period of wobble itself, if periodic (else NA). |
| `n_phase_shifts` | Count of discrete jumps in the shift signal. |
| `phase_shift_positions` | List of array coordinates where shifts occur. |
| `phase_shift_offsets` | Magnitude of each jump in bp. |
| `irregularity_score` | Variance of local-block features across the array. |
| `inter_monomer_identity` | Estimated mean sequence identity between different HOR slots (regime indicator). Reported only when class = HOR; NA otherwise. Derived from `R(1)` at `base_width` after accounting for the embedding's similarity floor. |
| `confidence` | Probabilistic score for the assigned class. |
| `class` | Summary label (see §3). |
| `n_segments` | Number of homogeneous segments (>1 only if phase shifts or domain transitions detected). |

Plus, as separate files:
- Consensus monomer (FASTA) and consensus HOR unit (FASTA, if `k>1`).
- Per-segment property table for arrays with `n_segments > 1`.
- Width-level feature table for diagnostics.

## 5. Inputs and outputs

### Inputs

```
arrays.fa                  # tandem-repeat array sequences, one per record
period_candidates.tsv      # array_id, period_bp, period_score, source
```

### Outputs

```
prefix.properties.tsv      # array-level property table (primary output)
prefix.segments.tsv        # per-segment properties where n_segments > 1
prefix.width_features.tsv  # per-width feature table (diagnostic)
prefix.consensus.fa        # consensus monomers and HOR units
prefix.diagnostics.json    # full structured diagnostics
prefix.*.png               # optional raster + shift-signal plots
```

## 6. Pipeline overview

```
FASTA + candidate periods
   │
   ▼
expand candidates → tested line widths (period ± neighborhood, plus divisors)
   │
   ▼
for each width w (cheap first, prune):
   │   column conservation       (background-corrected IC)
   │   ── if too low, mark unsupported and skip the rest ──
   │   row embeddings            (k-mer based)
   │   row autocorrelation R(k)
   │   edge-field statistics     (vertical/horizontal change rates)
   │   shift signal best_shift(r) along the array
   │       → mean_shift  (drift)
   │       → wobble_amplitude
   │       → breakpoints → n_phase_shifts
   │   local-block features → irregularity_score
   │   width-level class hint
   │
   ▼
combine width evidence → array-level properties
   │
   ▼
refine width using mean_shift (if non-zero, try w + mean_shift)
   │
   ▼
build consensus monomer and HOR unit
   │
   ▼
write TSV, JSON, FASTA, optional PNG
```

## 7. Feature definitions

### 7.1 Candidate width expansion

For each input period `p`: include `p ± N` (default `N=3`), plus divisors of `p` that fall within `[min_width, max_width]`, plus their neighborhoods. Deduplicate, cap to `max_widths_per_array`.

Divisors matter because the strongest period from the generator may be the HOR unit length, while the true base unit is a divisor.

### 7.2 Column conservation (background-corrected)

For each column at width `w`, count A/C/G/T over rows. Compute KL divergence against the array-wide base composition:

```
IC(col) = Σ_b p_b · log2(p_b / q_b)
```

where `p_b` is the column frequency and `q_b` is the array-wide frequency. This corrects for compositional bias (AT-rich satellites would otherwise show inflated IC against a uniform background).

Width-level features: `mean_column_IC`, `fraction_conserved_columns`.

### 7.3 Row embedding (k-mer)

For each row, compute canonical k-mer counts (default `k = 4`, 256 dimensions), then L2-normalise. Row similarity is the dot product.

Rationale: k-mer composition is position-invariant within a row, so the embedding is tolerant of small intra-row indels and of imperfect width registration. A positional bin embedding (the earlier proposal) is brittle to indels because a single insertion shifts every downstream bin.

If memory matters, feature-hash to 64–128 dimensions.

### 7.4 Row autocorrelation

```
R(k) = mean_{i} dot(emb[i], emb[i+k])   for k = 1 .. K
```

### 7.5 Edge-field statistics

Build two base-change indicator fields:

```
diff_x[r, c] = 1 if S[r·w + c]     ≠ S[r·w + c + 1]      else 0
diff_y[r, c] = 1 if S[r·w + c]     ≠ S[(r+1)·w + c]      else 0
```

Aggregate:

```
horizontal_edge_rate = mean(diff_x)             # within-row variability
vertical_edge_rate   = mean(diff_y)             # row-to-row variability at fixed col
```

At a correct width, vertical edges are sparse (columns conserve). At a wrong width or under chaotic structure, `vertical_edge_rate` rises.

Also compute the column-edge profile:

```
column_edge_rate[c] = mean_r diff_y[r, c]
```

Autocorrelate this profile along `c` (length `w`). A peak at lag `k` is an independent vote for HOR multiplicity `k`, derived from edge geometry rather than row embeddings. Two independent signals voting on the same `k` is a strong confidence boost.

### 7.6 Shift signal (unified drift / wobble / phase-shift detection)

For each adjacent row pair, compute similarity as a function of horizontal shift `s ∈ [-S, +S]` (default `S = 5`):

```
match(r, s) = fraction of c where S[r·w + c] == S[(r+1)·w + c + s]
best_shift(r) = argmax_s match(r, s)
```

`best_shift(r)` is a 1D signal along the array. Decompose it:

- **Global mean** `mean_shift = mean(best_shift)`. If non-zero and consistent, the tested width is off; refined width is `w + mean_shift`. Re-evaluate at the refined width.
- **Wobble amplitude** `wobble_amplitude = std(best_shift)` after removing breakpoint segments. Captures the smooth meandering seen in correct-width rasters.
- **Wobble periodicity**: 1D FFT (or autocorrelation) of `best_shift` after detrending. If a strong peak exists, report its period; else NA.
- **Breakpoints**: positions `r` where `|best_shift(r) − best_shift(r-1)| ≥ breakpoint_threshold` (default 3), or where adjacent-row similarity drops sharply for a single row then recovers. Count → `n_phase_shifts`. Positions in bp → `phase_shift_positions`. Magnitudes → `phase_shift_offsets`.

This one pass yields drift, wobble (amplitude and periodicity), and phase-shift segmentation.

### 7.7 Phase separation

```
best_k          = argmax_{k≥2} R(k)
same_phase      = R(best_k)
background      = median R(k') for k' near best_k but not a multiple
phase_separation = same_phase − background
```

HOR call requires `phase_separation ≥ threshold` (default 0.15) **in addition** to `R(best_k)` being high, to avoid calling homogeneous simple repeats as HORs.

### 7.8 Primitive multiplicity correction

If both `k = 12` and `k = 6` show high `R`, prefer the smaller:

```
for d in divisors(best_k):
    if R(d) ≥ R(best_k) − δ:    (default δ = 0.05)
        best_k ← d
```

### 7.9 Local irregularity

Split rows into blocks of size `B` (default `B = 100` or `L / 50`, whichever is larger). Per block compute: column IC, best lag, phase separation, mean shift. Aggregate as block-to-block variance:

```
irregularity_score = weighted sum of normalised block variances
```

High irregularity with strong HOR signal → `irregular_HOR`. High irregularity with weak HOR signal → `ambiguous`.

### 7.10 Segmentation

If `n_phase_shifts > 0`, segment the array at shift positions and re-run §7.2–7.9 per segment. If all segments share `base_width` and `k`, the array inherits the segment class (`HOR` or `simple_TR`) and phase shifts are reported as properties (`n_phase_shifts`, `phase_shift_positions`, `phase_shift_offsets`) with per-segment properties available in the diagnostics. If segments differ in `k` or `base_width`, report as `mixed`.

## 8. Classification logic

Decision rules, applied at array level after width refinement. Two central principles:

1. **Phase shifts are a property, not a class.** When the shift signal shows discrete jumps, the array is segmented and each segment classified independently. If all segments share the same structure (same `base_width`, same `k`), the array inherits that class and the phase shifts are reported as properties. Only if segments differ in `base_width` or `k` does the array become `mixed`.

2. **HOR call requires both periods to be statistically valid.** A HOR is called only when *both* the base period and the HOR-unit period are statistically valid. If only the longer period is valid, the call is `simple_TR` at the HOR-unit width (regime C — the HOR has degenerated).

```
if no width has column_IC ≥ threshold_min:
    class = ambiguous

elif n_phase_shifts > 0:
    # Segment at shift positions, classify each segment by the rules below,
    # then combine.
    if all segments share base_width and k:
        class = (class of any segment)         # HOR, simple_TR, or irregular_HOR
        # phase shifts reported as properties; do not modify the class label
    else:
        class = mixed

elif HOR_candidate_exists:
    # Candidate: some width w_b has best_k ≥ 2, phase_sep ≥ threshold,
    #            primitive-corrected, ≥ min_hor_copies
    if column_IC(w_b) ≥ threshold_hor_base               # base period valid
       AND column_IC(w_b · k) ≥ threshold_hor_unit       # HOR unit period valid
       AND irregularity_score low:
        class = HOR
    elif column_IC(w_b) ≥ threshold_hor_base
       AND column_IC(w_b · k) ≥ threshold_hor_unit
       AND irregularity_score elevated:
        class = irregular_HOR
    elif column_IC(w_b) < threshold_hor_base
       AND column_IC(w_b · k) ≥ threshold_simple_tr:
        # Regime C: base period has degenerated, only the HOR-unit period remains
        class = simple_TR  (at base_width = w_b · k, with note in `reason`)
    else:
        class = ambiguous

elif simple_TR width found (high column_IC, high R(1), low phase_sep):
    class = simple_TR

else:
    class = ambiguous
```

Threshold ordering: `threshold_hor_base < threshold_simple_tr ≈ threshold_hor_unit`. The base-width threshold is more permissive than the simple-TR or HOR-unit thresholds because at moderate inter-monomer divergence the columns are multimodal, which is correct for HOR but would not pass a simple_TR threshold.

The `reason` field in the output should record which regime was identified, so a downstream user can see when a `simple_TR` call came from regime C (degenerated HOR) versus regime A (true simple repeat), and whether phase shifts contributed to a `mixed` call.

## 9. Confidence score

A single sigmoid combining the main evidence components:

```
logit = α · phase_separation
      + β · (R(best_k) − R(unrelated_k))
      + γ · log10(n_complete_copies + 1)
      + δ · mean_column_IC
      − ε · irregularity_score
      − ζ · wobble_amplitude / w
      − η · |mean_shift| / w
confidence = sigmoid(logit)
```

Weights `α..η` calibrated on the synthetic test set so:
- Clean HOR cases land at `≥ 0.9`
- Clear ambiguous cases land near `0.5`
- Negative controls land at `≤ 0.2`

For `simple_TR` cases the `phase_separation` and `R(best_k)` terms are replaced by `R(1)` and `(R(1) − R(2))`.

State the formula and weights in the docs. This is the most-questioned number downstream.

## 10. Simulator

A standalone subcommand `simulate` that generates labelled FASTA test cases. The ground-truth labels are written alongside as a TSV matching the property table schema, so the test harness can compute property-level errors, not just class-level.

The simulator is specified in detail in the **structural taxonomy document, §5**:

- **§5.1 YAML schema** — structured input with a `templates` section (named, fully-specified slot-set objects) referenced by structural blocks. Replaces the earlier flat CLI parameter set. Key consequence: two HOR blocks referencing the same template share slot consensuses, making "one phase-shifted HOR" structurally distinct from "two HORs with the same k", which is otherwise ambiguous.

- **§5.2 Truth file schema** — two outputs. `truth.tsv` mirrors the detector's `properties.tsv` for properties that have a true generative value (omitting detector-derived measurements such as `column_conservation` or `confidence`), plus a mandatory `structural_expression` field carrying the grammar string. An optional `events_json` column carries structural variants (hybrids, insertions, inversions, local CNVs) that are not first-class detector properties in MVP.

- **§5.3 Wobble realisation** — integer insertions/deletions placed via a residual accumulator against a desired smooth shift curve. Produces non-integer average drift while every realised edit is a real base-level operation. Sub-base interpolation is not used.

- **§5.4 Test matrix** — minimum coverage matrix expressed in the §2 grammar, including a divergence sweep at fixed `k` to test the full HOR detection regime (A → B → C).

A summary CLI wrapper is acceptable for ad hoc single-array runs:

```
hordetect-image simulate --config test_case_T5.yaml --out runs/T5
```

but YAML is the primary input. CLI flags should map one-to-one onto YAML fields where they exist, to keep the truth-file generation logic single-sourced.

## 11. Testing plan

All tests run on simulator output with known ground truth. Property-level errors are computed in addition to class accuracy.

| # | Scenario | Expected properties |
|---|---|---|
| 1 | Simple TR, 171 bp × 1000 copies, 2% mut | `class=simple_TR`, `base_width=171`, `k=1` |
| 2 | Clean HOR, 171 bp base × k=12 × 200 copies, 2% mut, hor-divergence=0.1 | `class=HOR`, `base_width=171`, `k=12`, `hor_length=2052` |
| 3 | Same as #2, periods include 171 and 2052 | base width correctly identified as 171, not 2052 |
| 4 | True k=6, candidates include 6 and 12 | `k=6` after primitive correction |
| 5 | True width 171, candidate width 170 (one off) | `mean_shift ≠ 0`, width refined to 171, then clean call |
| 6 | Clean HOR with single phase shift mid-array | `class=HOR`, `n_phase_shifts=1`, position correct ±100 bp, `phase_shift_offsets` reported |
| 7 | Clean HOR with wobble amplitude 1.5 bp | `wobble_amplitude` reported, base call unaffected |
| 8 | HOR with periodic wobble (period = 500 rows) | `wobble_periodicity ≈ 500` |
| 9 | HOR with local deletion of one HOR unit | `class=irregular_HOR`, `irregularity_score` elevated |
| 10 | Mixed: half is k=12 HOR, half is simple TR same monomer | `class=mixed`, two segments reported |
| 11 | Random sequence, no repeat structure | `class=ambiguous`, low confidence |
| 12 | AT-rich (GC=0.2) simple TR | column IC and confidence not inflated by composition |
| 13 | Divergence sweep at k=4: hor-divergence ∈ {0.0, 0.05, 0.15, 0.35, 0.55} | Regime A (`simple_TR` at base), regime A boundary, regime B (`HOR`), regime B boundary, regime C (`simple_TR` at HOR-unit width). `inter_monomer_identity` tracks `1 - divergence` when HOR is called. |
| 14 | Regime-C check: k=4, hor-divergence=0.6, k=12, divergence=0.7 | Both correctly fall through to `simple_TR` at HOR-unit width, not `ambiguous`, and `reason` records "regime C / degenerated HOR". |

## 12. Output schemas

### `properties.tsv` (primary output)

```
array_id  length_bp  class  base_width_bp  hor_k  hor_length_bp  n_complete_copies
column_conservation  phase_separation  mean_shift_bp  wobble_amplitude_bp
wobble_periodicity_bp  n_phase_shifts  phase_shift_positions  phase_shift_offsets
irregularity_score  confidence  n_segments  reason
```

`phase_shift_positions` and `phase_shift_offsets` are comma-separated lists (or NA).

### `segments.tsv` (when `n_segments > 1`)

```
array_id  segment_id  start_bp  end_bp  class  base_width_bp  hor_k
column_conservation  phase_separation  wobble_amplitude_bp  irregularity_score
```

### `width_features.tsv` (diagnostic)

```
array_id  width_bp  rows  column_IC  fraction_conserved_columns
row_lag1_similarity  best_lag  best_lag_score  phase_separation
vertical_edge_rate  column_edge_autocorr_k  column_edge_autocorr_score
mean_shift_bp  wobble_amplitude_bp  n_phase_shifts  irregularity_score
class_hint
```

### `consensus.fa`

```
>array1.monomer  length=171
ACGT...
>array1.hor_unit  length=2052  k=12
ACGT...
```

### `diagnostics.json`

Full structured form: array-level properties, per-width features, per-segment features, breakpoint signals, R(k) curves, shift signal samples. Intended for programmatic consumption and plotting.

## 13. Implementation plan

### Milestone 1 — width expansion + column conservation
Read FASTA + periods; expand widths; compute background-corrected IC per width; write `width_features.tsv` with column-conservation columns only.

### Milestone 2 — row embeddings + autocorrelation
K-mer embeddings; `R(k)`; `best_lag`, `phase_separation`. Add to width feature table.

### Milestone 3 — edge field + shift signal
`diff_x`, `diff_y`, `column_edge_rate`, autocorrelation of `column_edge_rate`. Compute `best_shift(r)` along array; derive `mean_shift`, `wobble_amplitude`, `wobble_periodicity`, breakpoints. Width refinement when `mean_shift ≠ 0`.

### Milestone 4 — array-level classification + segmentation
Combine width evidence; primitive correction; segmentation when phase shifts detected; per-segment recomputation; confidence formula; emit `properties.tsv` and `segments.tsv`.

### Milestone 5 — consensus + diagnostics
Build consensus monomer and HOR unit by column-vote; emit `consensus.fa`; full `diagnostics.json`; optional PNG raster and shift-signal plots.

### Milestone 6 — simulator
Implement `simulate` subcommand with the parameters in §10. Generate all test cases in §11 and check property-level accuracy.

## 14. Dependencies (Rust)

```toml
clap         = "4"   # CLI
needletail   = "0.7" # FASTA parsing
rayon        = "1"   # parallelism over arrays and widths
serde        = "1"
serde_json   = "1"
csv          = "1"
ahash        = "0.8"
image        = "0.25" # optional PNG output
anyhow       = "1"
thiserror    = "2"
log          = "0.4"
env_logger   = "0.11"
rand         = "0.8"  # simulator
rand_distr   = "0.4"  # simulator
```

Deferred: `rustfft` (only if direct lag computation becomes a bottleneck), `ndarray` (if 2D numeric pipelines grow). No `opencv`, no `imageproc` for MVP.

## 15. Out of scope for MVP

- Hilbert curve representation
- CNN / learned classifier
- Dense `L²` dotplot
- Whole-genome tandem-array discovery (input is pre-annotated)
- Offset search / non-tandem flank trimming (input is pre-annotated)
- Reverse-complement orientation search (input is canonical)
- Variant HOR reconstruction (which monomer occupies which slot)
- Species-specific centromere annotation
- 2D FFT, Hough/Radon transforms

These can be added in v2 if line-width features prove predictive on real data.

## 16. Complexity

Per array, with `L` = length, `W` = tested widths, `K` = max HOR multiplicity, `d` = embedding dim, `S` = shift range:

```
column stats:        O(W · L)
edge field:          O(W · L)
row embeddings:      O(W · L)
row autocorrelation: O(Σ_i (L/w_i) · K · d)
shift signal:        O(Σ_i (L/w_i) · S · w_i) = O(W · L · S)
```

Dominant terms scale with `W · L`, not `L²`. With `W ≈ 20`–`50` from the period generator, this is tractable for multi-megabase arrays.

## 17. One-paragraph summary

The tool characterizes tandem-repeat arrays by evaluating candidate periods as line widths. For each width, it extracts column conservation (with background correction), k-mer-based row autocorrelation, edge-field statistics, and a per-row shift signal that simultaneously yields drift, wobble amplitude, wobble periodicity, and discrete phase-shift breakpoints. Width evidence is combined into a property vector covering base width, HOR multiplicity, copy number, wobble, phase-shift positions, and irregularity; a class label and confidence score are derived. A simulator generates labelled test data with controllable mutations, indels, wobble, phase shifts, and structural disruptions for property-level evaluation. The method is O(W · L), uses pre-annotated tandem arrays as input, and avoids dense dotplots, CNNs, and Hilbert curves.

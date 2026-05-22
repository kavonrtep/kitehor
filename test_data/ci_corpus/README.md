# `ci_corpus/` — diverse CI fixture for the rule-proto pipeline

Curated 13-record corpus exercising five of the eight `combined_class`
values: `hor, pure_ssr, tr, tr_with_nested_tr, tr_with_ssr`.

The corpus is small (~few hundred KB) by design so end-to-end
`kitehor analyze` against it runs in well under a second and can be
gated in CI.

Provenance is documented in `manifest.tsv`. The synthetic-fixture
records (rDNA-like + sub-TR + various edge cases) come verbatim
from `tools/rule_proto/subrepeat/synthetic.fasta` (regenerable via
`tools/rule_proto/subrepeat/make_fixtures.py`). The `sim_hor_*`
records are reproducible via `kitehor simulate` with the parameters
listed in `manifest.tsv`.

## Regenerating

```bash
KITE=./target/release/kitehor
TMP=$(mktemp -d)
$KITE simulate --monomer-size 100 --multiplicity 3 --copies 80 \
    --sub-rate-intra 0.04 --sub-rate-inter 0.02 \
    --case-id sim_hor_k3 --out $TMP/sim_hor_k3.fa
$KITE simulate --monomer-size 150 --multiplicity 5 --copies 60 \
    --sub-rate-intra 0.05 --sub-rate-inter 0.03 \
    --case-id sim_hor_k5 --out $TMP/sim_hor_k5.fa
$KITE simulate --monomer-size 200 --multiplicity 7 --copies 40 \
    --sub-rate-intra 0.06 --sub-rate-inter 0.04 \
    --case-id sim_hor_k7 --out $TMP/sim_hor_k7.fa
cat tools/rule_proto/subrepeat/synthetic.fasta \
    $TMP/sim_hor_k3.fa $TMP/sim_hor_k5.fa $TMP/sim_hor_k7.fa \
    > test_data/ci_corpus/sequences.fasta
```

## Missing classes

The current corpus does not yet cover `hor_with_ssr`,
`tr_with_subrepeat`, or `unresolved`. Adding them requires more
targeted synthesis (e.g., a HOR whose monomer contains a SSR for
`hor_with_ssr`; a TR with an internal 2-3-copy duplication of an
F-length segment for `tr_with_subrepeat`). When that synthesis lands,
extend `manifest.tsv` and `tests/analyze_ci_corpus.rs`.

# Release runbook

Tag-driven pipeline. One annotated tag (`v<MAJOR>.<MINOR>.<PATCH>`)
on `main` fires both the GitHub release and the conda upload.

## Workflows

| File | Trigger | Jobs |
|---|---|---|
| [`.github/workflows/ci.yml`](../.github/workflows/ci.yml) | push to `main`, PR | `fmt`, `clippy`, `test` (linux-x86_64) |
| [`.github/workflows/release.yml`](../.github/workflows/release.yml) | tag `v*.*.*` (also `workflow_dispatch`) | `check-tag` → `build` → `publish` → `conda` |

The release pipeline is modelled on
[`kavonrtep/dottir`'s setup](https://github.com/kavonrtep/dottir/tree/main/.github/workflows).

## Cutting a release

```bash
# 1. Bump version in Cargo.toml (single source of truth).
#    The release workflow's `check-tag` step asserts the tag and the
#    Cargo.toml version match — fast-fail if they don't.
$EDITOR Cargo.toml
cargo build --release        # refresh Cargo.lock
git commit -am "release: v0.9.X"

# 2. Tag + push.
git tag v0.9.X
git push origin main
git push origin v0.9.X
```

The push of the tag triggers `release.yml`:

1. **`check-tag`** — parses `${GITHUB_REF_NAME}`, asserts
   `cargo_version == tag_minus_v_prefix`. Hard fail otherwise.
2. **`build`** — `cargo build --release --locked` on
   `ubuntu-22.04`, target `x86_64-unknown-linux-gnu`. Strips the
   binary, packs `kitehor + README + LICENSE-*` into a
   `.tar.gz`, emits a `.sha256` sidecar.
3. **`publish`** — `gh release create v0.9.X --generate-notes
   --verify-tag` with the tarball + a combined `SHA256SUMS`.
4. **`conda`** — `conda build conda/kitehor/ --output-folder
   conda-out` (with `KITEHOR_VERSION=0.9.X`), then
   `anaconda upload --user petrnovak --label main --force` the
   resulting `linux-64/kitehor-*.conda` (or `.tar.bz2`). The build
   step runs inside conda-forge's
   `quay.io/condaforge/linux-anvil-cos7-x86_64` docker image
   (CentOS 7, glibc 2.17) so the resulting binary stays portable
   across all current LTS distros — Ubuntu 20.04, CentOS / RHEL 8,
   Debian 11. Running `cargo install` directly on the `ubuntu-22.04`
   runner pulls in glibc 2.35 symbols (`GLIBC_2.32` / `2.33` / `2.34`
   refs) that break on anything older — that's what bit v0.9.2; see
   [`docs/kitehor_upstream_issues.md`](kitehor_upstream_issues.md).

## Secrets required

| Secret | Where | Used by |
|---|---|---|
| `ANACONDA_API_TOKEN` | Settings → Secrets → Actions | `release.yml → conda` step |
| `GITHUB_TOKEN` | Auto-injected by GitHub | `release.yml → publish` step |

The user's setup already has `ANACONDA_API_TOKEN` provisioned.

## Re-running against an existing tag

`workflow_dispatch` is wired with a single `tag` input, so a release
can be re-attempted (e.g. after the conda upload fails) without
deleting and re-creating the tag:

GitHub → Actions → Release → "Run workflow" → enter `v0.9.X`.

## Local pre-flight

The repo ships tracked git hooks under [`.githooks/`](../.githooks/)
that mirror the CI gates. **One-time install per clone**:

```bash
git config core.hooksPath .githooks
```

After that:

- **`pre-commit`** runs `cargo fmt --check` + `cargo clippy -- -D warnings`
  on every commit that touches `.rs` / `.toml` / `.lock`. ~5 s
  incremental. Bypass with `git commit --no-verify`.
- **`pre-push`** runs `cargo test --release --locked` on every push.
  Bypass with `git push --no-verify`.

The same gates also runnable manually:

```bash
cargo fmt --all --check
cargo clippy --release --all-targets --locked --no-deps -- -D warnings
cargo test --release --locked
```

All tests use small datasets shipped in `test_data/` or
`tools/rule_proto/fixtures/` — no external downloads, no
large-corpus dependencies.

## Conda recipe

[`conda/kitehor/meta.yaml`](../conda/kitehor/meta.yaml) — the
version is templated from the workflow's `KITEHOR_VERSION` env var
(default `0.0.0.dev` when built locally without the env var). The
build script is a single `cargo install --locked --no-track --bin
kitehor --root $PREFIX --path .`. The `test:` block runs `--version`,
`--help`, and an end-to-end `analyze` on a `kitehor simulate`-built
synthetic FASTA.

## First release: v0.9.2

v0.9.0 and v0.9.1 attempts both failed at the conda job:

- **v0.9.0**: multi-line `{# … #}` jinja comment in `meta.yaml`
  that conda-build's parser rejected. Fixed in `fdf78c2`.
- **v0.9.1**: `{{ stdlib('c') }}` + `{{ compiler('rust') }}` macros
  expanded to bare placeholders (`c_linux-64`) because conda-forge's
  variant config wasn't being applied to the build env. Fixed by
  mirroring the dottir pattern: drop `stdlib('c')`, drop the rust
  compiler macro, depend on `rust >=1.85` as a plain package.

v0.9.2 is the first release that ships a published conda package.

### v0.9.2 — glibc-too-new follow-up

The v0.9.2 conda package builds and uploads cleanly, but the resulting
binary inherits the runner's glibc 2.35 symbols and won't run on any
host with glibc < 2.34 (Ubuntu 20.04, CentOS 8, Debian 11). Reported in
[`docs/kitehor_upstream_issues.md`](kitehor_upstream_issues.md). The
v0.9.3 release ships the fix.

## v0.12.0 — three new subcommands + cascade gates + full doc coverage

Non-breaking release. Adds three independent observation subcommands
(`rescore`, `report`, `irregularity`), tightens the rule-classify and
summary-merge cascade with two new threshold gates, and ships
per-column documentation for every TSV the tool emits.

What it ships:

- **feat(rescore)**: new `kitehor rescore` subcommand — banded
  semi-global edit-distance scoring of sampled adjacent tile pairs
  per kite peak. Appends **9 columns** to the input peaks TSV:
  identity stats (`identity_med`, `identity_iqr`, `identity_p25`,
  `identity_n`), shift diagnostics (`shift_med`,
  `shift_consistency`), two derived flags (`phantom`, `subrepeat`),
  and the standalone `coverage_frac`. Phantom catches sub-period
  harmonics via shifted alignments; subrepeat catches localised
  short motifs inside the founder via bimodal identity
  distributions, founder-gated to keep FP rate ≈ 0.4 % on
  ground_truth_v2. Period-relative auto-band (`max(20, 2·slop,
  ⌈0.02·P⌉)`) keeps long monomers from saturating. Additive: no
  downstream stage changes; cascade still decides on `score2_norm`.
  Calibration: 0 phantom/subrepeat FPs on true HOR-unit periods;
  HOR-unit `identity_med` beats monomer `identity_med` on 100 % of
  clean HORs and 97.5 % overall (1300-case ground_truth_v2). 62
  unit tests + 2 integration tests.
- **feat(report)**: new `kitehor report` subcommand — observation-
  only whole-array TSV. 20 columns: kite peaks (raw + clustered),
  SSR coverage, irregularity metrics (indel events + drift +
  dropout rate), and array metadata. No `combined_class`, no
  rule-classify verdicts. Designed as a sibling to `analyze` for
  consumers who want raw numbers to filter on themselves rather
  than the cascade's calls.
- **feat(irregularity)**: new `kitehor irregularity` subcommand —
  Rust port of the v2 distance-residual + phase-bin indel-event
  prototype. 14-column TSV split into the Approach-6 indel /
  dropout signals (`indel_event_count`, `indel_burden_pct`,
  `indel_max_shift_bp`, `indel_drift_bp_per_kb` versus
  `dropout_rate_per_pair`). Available as a standalone subcommand
  and as the irregularity backend powering `report` columns 13–20.
- **feat(rule_classify) — new gate**: HOR verdicts now require at
  least `--min-tile-copies` (default 6) complete tile copies in the
  array. Records with too few copies fall back to `simple_tr` /
  `unresolved`. Surfaces the constraint that an HOR call needs
  enough tile copies for the clustering bump to be meaningful.
- **feat(summary) — density gate**: `tr_with_subrepeat` cascade
  decision additionally requires `tv_density ≥
  --subrepeat-density-min` (default 0.4, tuned down from 0.7).
  Surfaces a new `subrepeat_coverage_pct` column in `summary.tsv`
  (column 32) so downstream filters can read the underlying
  density. Detail: [`docs/irregularity_and_subrepeat_v0_12.md`](irregularity_and_subrepeat_v0_12.md) §1.
- **feat(simulator)**: legacy `simulate` gains discrete indel
  events and compound subrepeat monomers; useful for stress-
  testing the new irregularity backend on known event tracks.
- **docs**: every subcommand now has a reference page with
  per-column descriptions for every TSV it emits:
  - [`docs/rule_proto.md`](rule_proto.md) — 7-stage column tables
    (kite, peaks, verdicts, tandem_validate, ssr, ssr.regions,
    summary)
  - [`docs/rescore.md`](rescore.md) + [`docs/report.md`](report.md)
    (existing; refreshed)
  - [`docs/irregularity.md`](irregularity.md) (new)
  - [`docs/detect.md`](detect.md) (new)
  - [`docs/synth.md`](synth.md) (new)
  - [`docs/simulate.md`](simulate.md) (new)
  - [`docs/onboarding_pipelines.md`](onboarding_pipelines.md) —
    rescore + report orientation map
  - `README.md` subcommand table updated to link every reference
- **schema(summary.tsv)**: gains `subrepeat_coverage_pct` (column
  32). Column count: 32 → 33. Existing column indices unchanged.

Non-breaking notes for operators:

- `combined_class` values are unchanged from v0.11.0.
- The new subcommands are independent of `analyze` — `analyze`
  still emits its 7 per-stage TSVs and nothing else. Run
  `rescore`, `report`, `irregularity` separately when you want the
  extra signal.
- Rule-classify's `min_tile_copies` gate may flip a small number
  of borderline HOR calls to `simple_tr` / `unresolved` on arrays
  with very few tile copies (below 6 by default). Override via
  `--min-tile-copies 0` to restore v0.11 behaviour.
- Summary's `subrepeat_density_min` gate (0.4) may suppress a
  small number of borderline `tr_with_subrepeat` calls. Override
  via `--subrepeat-density-min 0` to restore v0.11 behaviour.

Pre-flight passed:
- 403 unit + integration tests pass, 1 ignored
  (`tandem_validate_python_parity`, opt-in;
  `cargo test --release --locked`).
- `cargo clippy --release --all-targets --locked --no-deps -- -D warnings` clean.
- `cargo fmt --all --check` clean.
- Local `kitehor --version` reports `0.12.0`.

To ship:

```bash
git tag v0.12.0     # already created by the /release skill
git push origin main
git push origin v0.12.0
```

## v0.11.0 — SSR cascade fix + 2 new classes (breaking)

**Breaking release.** Fixes a `pure_ssr` over-call bug from v0.10
*and* expands the combined_class enumeration from 7 → 9 so every
verdict category has a parallel `_with_ssr` partner.

What it ships:

- **fix(summary)**: SSR cascade no longer fires `pure_ssr` for
  arrays as low as ~4% SSR coverage. Under v0.10 the `ssr-scan`
  `consensus_single` path (used when the kite top period contains
  a known SSR motif) overwrote `dominant_motif_coverage_pct` with
  the candidate monomer's *self*-coverage on a synthetic dimer
  (~100% by construction), and the cascade's `dom_pct ≥ 80` check
  then over-triggered. The cascade now reads
  `ssr_raw_total_coverage_pct` (the array-scale total from the raw
  scanner) exclusively for all SSR decisions. The per-record
  `ssr_flag` column is recomputed from the same raw total so it
  stops contradicting the call. The consensus path's
  `dominant_motif*` fields remain in `summary.tsv` as informational
  labels.
- **feat(summary)**: two new combined_class values so the cascade
  covers all 4 verdict×SSR combinations: `unresolved_with_ssr`
  (unresolved + SSR ≥ 30%) and `tr_with_subrepeat_with_ssr`
  (tandem_validate `localized_subrepeat` + SSR ≥ 30%). The
  has-SSR threshold is exposed as the new
  `--ssr-has-pct-threshold` flag on `analyze` + `summary-merge`;
  default 30% matches `--ssr-flag-threshold-pct` so the
  recomputed flag and the cascade agree by default.
- **schema(summary.tsv)**: `ssr_raw_total_coverage_pct` promoted
  from "optional SSR diagnostic" to always-emitted since it now
  drives the cascade.

Breaking changes operators need to know:

- `combined_class` may differ from v0.10 for any record that
  previously fired `pure_ssr` on a consensus_single result. On the
  IPIP200579 pangenome corpus (3024 records), 53 records moved out
  of `pure_ssr` (99 → 46); they redistributed into `tr` (+7),
  `tr_with_ssr` (+10), `unresolved_with_ssr` (new, 78 records),
  and `tr_with_subrepeat_with_ssr` (new, 2 records).
- Library `summary::combined_class()` signature changed from
  `(verdict, ssr_flag, ssr_dom_pct, tv_decision, &cfg)` to
  `(verdict, ssr_raw_total_pct, tv_decision, &cfg)`.
- `summary::Config` gains `ssr_has_pct_threshold: f64` (default 30.0).
- Two new combined_class enum values to handle in downstream code.

Pre-flight passed:
- 379 unit + integration tests pass, 3 ignored (`cargo test --release --locked`).
- `cargo clippy --release --all-targets --locked --no-deps -- -D warnings` clean.
- `cargo fmt --all --check` clean.
- Local `kitehor --version` reports `0.11.0`.

To ship:

```bash
git tag v0.11.0
git push origin main
git push origin v0.11.0
```

## v0.10.0 — unified tandem_validate stage (breaking)

**Breaking release.** The CLI surface lost two subcommands and one
`combined_class` value; the `summary.tsv` schema lost nine columns and
gained ten. Pipelines pinning v0.9.x output formats need to be updated
before adopting v0.10.0.

What it ships:

- **feat(tandem_validate)**: unified spatial-localization subrepeat
  detector — port of `tools/rule_proto/tandem_validate.py` (spec v5).
  One density + spatial- + phase-contrast check replaces the prior
  two-stage `subrepeat-scan` + `hor-validate` pair, testing both
  within-tile and array-scale heterogeneity at one geometric scale
  (window capped at host, candidates capped at `host/3`). New
  `kitehor tandem-validate` subcommand exposes all 13 thresholds as
  flags. See `docs/new/tandem_validate_spec.md` for the algorithm
  and `docs/new/tandem_validate_port_plan.md` for the rollout plan.
- **refactor(analyze)**: pipeline collapses from 5 → 4 post-classify
  stages and the cascade from 8 → 7 classes (`tr_with_nested_tr`
  retired; both prior subrepeat-style triggers merge into
  `tr_with_subrepeat`). The `analyze` orchestrator now writes 7
  per-stage TSVs instead of 9. Net diff: −1402 LOC after deleting
  `src/subrepeat/` and `src/hor_validate/`.
- **test(tandem_validate)**: ignored-by-default Python-parity test
  (`tests/tandem_validate_python_parity.rs`) — runs both the Rust
  and Python implementations on a 6-record synthetic fixture and
  asserts `decision_hint` matches. Run before tagging via
  `cargo test --release --test tandem_validate_python_parity -- --ignored`.
- **docs**: README, CLAUDE.md, and `docs/rule_proto.md` updated for
  the 7-stage / 7-class layout. `docs/new/rule_proto_impl_plan.md`
  retains the original port intent but carries a "Historical
  document — partially superseded by v0.10" header at the top.
- **fix(ci)**: cos7 conda build now runs from a writable workdir and
  publishes idempotently (no-op on a re-tag that ships the same
  `.conda` artifact). Carried over from post-v0.9.3 work.

Breaking changes operators need to know:

- `kitehor subrepeat-scan` and `kitehor hor-validate` subcommands
  removed.
- `kitehor summary-merge` — `--subrepeat` and `--within-tile` flags
  replaced by a single `--tandem-validate`.
- `tr_with_nested_tr` is no longer a possible `combined_class`
  value (semantically equivalent records now fire
  `tr_with_subrepeat`; some former `tr_with_nested_tr` records
  legitimately fall through to `tr` or `unresolved` per the v5
  out-of-scope list in the spec).
- `summary.tsv`: drops `length_bp`, all `subrepeat_*` columns, and
  all `density_hint` / `founder_density` / `phase_contrast` /
  `density_n_windows` columns. Gains 10 `tv_*` columns from the new
  detector. Join on `record_id` against `<prefix>.kite.tsv` if you
  need `length_bp`.

Pre-flight passed:
- 374 unit + integration tests pass, 3 ignored (`cargo test --release --locked`).
- `cargo clippy --release --all-targets --locked --no-deps -- -D warnings` clean.
- `cargo fmt --all --check` clean.
- Local `kitehor --version` reports `0.10.0`.

To ship:

```bash
git tag v0.10.0
git push origin main
git push origin v0.10.0
```

## v0.9.3 — first portable conda binary

Behaviour-equivalent to v0.9.2 plus:

- **fix(conda)**: the `conda` job now runs `conda build` inside
  `quay.io/condaforge/linux-anvil-cos7-x86_64` (CentOS 7, glibc 2.17),
  the canonical conda-forge build env. No recipe changes — only the
  build host. Resulting binary is portable across every current LTS
  distro (Ubuntu 20.04 / Debian 11 / CentOS 8 / RHEL 8 and newer).
- **feat(kite)**: `kite-periodicity` and `analyze` gain
  `--periodogram <PATH>` — a FASTA-like bundle of the per-record
  neighbour-distance histogram + smoothed background for plotting
  (mirrors TideCluster's in-memory `profile_list`). The previous
  `--dump-profile <DIR>` flag (per-record sparse TSVs) is removed.

Pre-flight passed:
- 288 unit + 60+ integration tests pass (`cargo test --release --locked`).
- `cargo clippy -- -D warnings` clean.
- `cargo fmt --check` clean.
- Local `kitehor --version` reports `0.9.3`.

To ship:

```bash
git tag v0.9.3
git push origin main
git push origin v0.9.3
```

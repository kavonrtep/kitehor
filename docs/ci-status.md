# CI / release pipeline — status & plan

Last updated: 2026-05-15.

This document tracks the GitHub Actions + conda packaging work for
`kitehor`. It is both a **state report** (what is done, what is left)
and the **canonical reference** for the original implementation plan.

## Phase tracker

| Phase | Subject | State | Commit |
|---|---|---|---|
| **1** | Code prep: bake models, clear clippy, format, pin toolchain, polish Cargo.toml | **Done** | `0b93097` |
| **2** | `ci.yml` — fmt + clippy + test + smoke on push & PR | **Done** | `f8c91bb` |
| **3** | `release.yml` + `conda-release.yml` + `conda/kitehor/meta.yaml` | **Done** | `600eb5d` |
| **4** | First release tag (`0.1.0`): one-time secrets setup, push tag, verify both pipelines | **Pending** | — |

Local validation gates at the time of writing:

| | result |
|---|---|
| `cargo build --release` | 26 MB binary (22 MB of model JSON baked in) |
| `cargo test --release --lib` | 47/47 |
| `cargo test --release --test integration_smoke -- --ignored` | 1/1 |
| `cargo clippy --release --all-targets --no-deps -- -D warnings` | clean |
| `cargo fmt --check` | clean |
| Local smoke verdicts on `test_data/smoke/` | tandem_pure→tandem · hor_k3→hor (k=3, founder=100, tile=300) · hor_k5→hor (k=5, founder=150, tile=750) |

## Locked decisions

| | Choice | Rationale |
|---|---|---|
| Versioning | Single source: `Cargo.toml` `version`. Tag format `[0-9]+.[0-9]+.[0-9]+*` (bare numeric, no `v` prefix). | Matches TideCluster pattern; one tag fires both workflows. |
| Models | Baked into binary via `include_bytes!`. ~26 MB binary, fully self-contained. | Zero runtime path-resolution code; conda recipe collapses to one `cargo install` line; eliminates version skew. |
| Platforms | linux-64 (gnu) + linux-64 (musl, static). | gnu covers the default conda env; musl gives a portable static binary. macOS deferred. |
| Conda source | `source: path: ../..` (build from checked-out tree). | Matches TideCluster pattern; no tarball/sha256 management. `conda-release.yml` runs independently of `release.yml`. |
| CI gates | `cargo test`, `cargo fmt --check`, `cargo clippy -- -D warnings`. | Catches regressions, style drift, and lint regressions. `cargo audit` deferred (advisory-only later). |

## Implementation plan (with completed sections marked)

### Phase 1 — Code changes ✓ done in `0b93097`

- **Bake models with `include_bytes!`**
  - `src/classifier.rs`: added `BAKED_HOR_MODEL` and `BAKED_K_MODEL` constants;
    new `RandomForest::load_json_bytes(&[u8])` shares parser with the file-path
    variant. Removed the now-unused `ModelsCfg` struct from `ClassifierConfig`.
  - `src/main.rs`: default loader uses baked bytes; `--hor-model` / `--k-model`
    flags still override with a path.
  - `config/classifier.toml`: dropped the `[models]` table.
- **Cleared 23 clippy errors** under `-D warnings`: abs_diff, byte-string
  literals, `iter::enumerate`, `unwrap` after `is_some`, useless casts,
  doc-list indentation, unused lifetimes.
- **`cargo fmt`** over the tree (10 files).
- **`rust-toolchain.toml`** pins `channel = "1.95"` + rustfmt + clippy components.
- **`Cargo.toml`** polish: `repository`, `homepage`, `documentation`,
  `keywords`, `categories`, `exclude` block (keeps test_data, ground_truth,
  tools, docs, tests out of any future `cargo package` artefact).

### Phase 2 — `ci.yml` ✓ done in `f8c91bb`

`.github/workflows/ci.yml`. Single `ubuntu-latest` job on push/PR to `main`:

1. checkout
2. `rustup show` (installs the 1.95 toolchain pinned in `rust-toolchain.toml`)
3. `Swatinem/rust-cache@v2` (cargo registry + `target/`)
4. `cargo fmt --all -- --check`
5. `cargo clippy --release --all-targets --no-deps -- -D warnings`
6. `cargo test --release --lib`
7. `cargo build --release`
8. `cargo test --release --test integration_smoke -- --ignored`

Concurrency group `ci-${{ github.ref }}` cancels stale runs on the same ref.
30-minute timeout. Cold-cache run estimated at ~5 min; warm runs ~1.5 min.

### Phase 3 — release pipeline + conda recipe ✓ done in `600eb5d`

Three new files, all triggered by `push: tags: ['[0-9]+.[0-9]+.[0-9]+*']`.

**`.github/workflows/release.yml`** — GitHub Release with binary tarballs.

| Job | Steps |
|---|---|
| `version-check` | Asserts tag equals `Cargo.toml` `[package].version`. Fast-fail. Outputs `tag` for downstream jobs. |
| `build` (matrix: `x86_64-unknown-linux-gnu` + `x86_64-unknown-linux-musl`) | `rustup target add`; `apt install musl-tools` for the musl branch; cache; `cargo build --release --locked --target …`; strip; tarball `kitehor-<tag>-<target>.tar.gz` containing binary + LICENSE-MIT + LICENSE-APACHE + README.md; emit `.sha256` sibling. |
| `publish` | Downloads matrix artefacts, builds combined `SHA256SUMS`, `gh release create --generate-notes --verify-tag` with both tarballs + `SHA256SUMS` attached. |

`workflow_dispatch(tag)` supported for re-running a release without re-tagging.

**`.github/workflows/conda-release.yml`** — TideCluster-pattern conda upload.

| Step | Purpose |
|---|---|
| Tag = Cargo.toml version check | Fast-fail |
| `ANACONDA_API_TOKEN` set check | Fast-fail (don't burn 30 min then fail at upload) |
| `conda-incubator/setup-miniconda@v3` (miniforge, mamba, strict channels: `conda-forge,kavonrtep`) | Build env |
| `mamba install -y conda-build anaconda-client "setuptools<81"` | `setuptools<81` pin: `anaconda-client`'s `binstar_client` still imports `pkg_resources` which setuptools 81 removed (TideCluster playbook §7 pitfall #4). |
| `cp Cargo.toml conda/kitehor/Cargo.toml` | `meta.yaml`'s `load_file_regex(from_recipe_dir=True)` reads the version from a Cargo.toml that lives inside the recipe dir — conda-build copies the recipe to a temp workdir before rendering. |
| `conda-build -c conda-forge --output-folder build_out conda/kitehor` | Produces `build_out/linux-64/kitehor-<version>-h….conda` |
| `anaconda upload --user kavonrtep --label main` | Publish to `anaconda.org/kavonrtep/main` |
| `actions/upload-artifact` (`if: always()`) | Debug-only artefact upload of the built package |

**`conda/kitehor/meta.yaml`** — Rust adaptation of the conda-forge example.

- Version: parsed from Cargo.toml via `load_file_regex` with the anchored
  pattern `\nversion\s*=\s*"([^"]+)"` so `rust-version = "..."` later in the
  file does not match the package version.
- `source: path: ../..`
- `build.script: cargo install --locked --no-track --bin kitehor --root $PREFIX --path .`
- `requirements.build`: `stdlib('c')`, `compiler('c')`, `compiler('rust')`.
- `test.commands`:
  - `kitehor --version | grep -F "<version>"` (catches stale lockfile / wrong checkout)
  - `kitehor --help | grep -qi 'kitehor'`
  - `kitehor kite-periodicity --help | grep -qi 'classify'`
  - End-to-end: `simulate --monomer-size 100 --multiplicity 3 --copies 80` then
    `kite-periodicity --classify` and assert the output TSV is non-empty
    with the expected `sim_0000` row. Exercises the baked-in classifier in
    the installed package.

### Phase 4 — First release (pending)

See **§ Release runbook** below.

## Release runbook

The first release also doubles as the validation of the full pipeline. Once
green, subsequent releases reduce to "bump Cargo.toml, tag, push tag".

### One-time setup (do once per repo)

1. **Move the staged repo out of the sandbox** (we currently live at
   `kite2/kitehor/`; needs to be a sibling once the container restriction
   lifts):
   ```bash
   cp -a kite2/kitehor/. ~/PycharmProjects/kitehor/
   cd ~/PycharmProjects/kitehor
   ```
2. **Add the GitHub remote and push `main`**:
   ```bash
   git remote add origin git@github.com:kavonrtep/kitehor.git
   git push -u origin main
   ```
   First push triggers `ci.yml` — confirms the workflow runs cleanly against
   hosted runners.
3. **Create the anaconda.org API token**: anaconda.org → Settings → Access →
   API tokens → create token with **Write to API** scope. Copy the token.
4. **Add it as a GitHub Actions secret**: github.com/kavonrtep/kitehor →
   Settings → Secrets and variables → Actions → New repository secret.
   Name: `ANACONDA_API_TOKEN`. Paste the token.

### Cut the release

```bash
git tag 0.1.0 -m "Initial release"
git push origin 0.1.0
```

Both `release` and `conda-release` workflows fire in parallel.

### Expected outcomes

- **GitHub Release `0.1.0`** with three attached files:
  - `kitehor-0.1.0-x86_64-unknown-linux-gnu.tar.gz` (~7.5 MB)
  - `kitehor-0.1.0-x86_64-unknown-linux-musl.tar.gz` (similar)
  - `SHA256SUMS`
- **`anaconda.org/kavonrtep/kitehor`** package `0.1.0` for `linux-64`.
- A fresh env should resolve:
  ```bash
  mamba create -n kitehor-smoke -c kavonrtep -c conda-forge kitehor=0.1.0
  mamba activate kitehor-smoke
  kitehor --version    # → "kitehor 0.1.0"
  ```

### If something fails

```bash
git tag -d 0.1.0
git push --delete origin 0.1.0
# fix the issue on main; retag.
```

Don't use the same tag twice without an intentional reason — anaconda.org
will reject a duplicate upload at the same `--label main`.

## Engineering notes / pitfalls

- **`serde_json` float parsing** — the `float_roundtrip` cargo feature is
  required (already in `Cargo.toml`). Default parser is fast but not
  correctly-rounded; without it, RF split values parse 1 ULP off and tree
  traversals occasionally flip leaf direction. This was caught during the
  port and recorded in user-level memory.
- **`setuptools<81` pin** — needed for `anaconda-client`'s `binstar_client`
  which still imports `pkg_resources`. Already pinned in
  `conda-release.yml`.
- **`conda-build` vs `conda mambabuild`** — we use plain `conda-build`.
  `mambabuild` requires `boa`'s conda plugin to register correctly and is
  flaky on hosted runners. Same builder underneath; no plugin layer.
- **Version regex** — `\nversion\s*=` anchors to start-of-line so
  `rust-version = "..."` (later in `Cargo.toml`) does not falsely match.
  First newline-anchored `version =` is the `[package]` one.
- **Sandbox quirk** — `kitehor/` currently lives at `kite2/kitehor/`. When
  pushing to GitHub, work from a moved-out copy so the parent `kite2`
  repo's `.git` does not interfere.
- **First tag drift risk** — if `ANACONDA_API_TOKEN` is missing on the
  first tag, `conda-release.yml` fails fast at the token-check step (we
  added an explicit check). `release.yml` is unaffected and still produces
  the GitHub Release.

## Document layout

`docs/` is the canonical home for all project documentation that is not
the public-facing `README.md` (root) or the local dev guide `CLAUDE.md`
(root). Historical / pre-trim design notes live under `docs/archive/`
which is gitignored.

```
docs/
├── ci-status.md     ← this file
└── archive/         (gitignored: historical design docs)
```

Add new docs as siblings of this file (`docs/<topic>.md`).

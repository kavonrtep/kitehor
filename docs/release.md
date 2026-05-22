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
   resulting `linux-64/kitehor-*.conda` (or `.tar.bz2`).

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

## First release: v0.9.1

v0.9.0 was attempted but the conda job failed at the recipe-render
step (multi-line `{# … #}` jinja comment that conda-build's parser
rejected — fixed in `fdf78c2`). v0.9.1 is the first release that
ships a published conda package.

Pre-flight passed:
- 352 tests pass, 0 fail (`cargo test --release --locked`).
- `cargo clippy -- -D warnings` clean.
- `cargo fmt --check` clean.
- Local `kitehor --version` reports `0.9.1`.

To ship:

```bash
git tag v0.9.1
git push origin v0.9.1
```

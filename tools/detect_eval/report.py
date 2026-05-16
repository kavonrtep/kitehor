#!/usr/bin/env python3
"""Render an interactive HTML dashboard for a kite → detect run.

Inputs are the artefacts produced by:

  kitehor kite-periodicity <fastas...> -o <kite.tsv> --classify \
      --emit-periods <periods.tsv>
  kitehor detect-batch --fasta-dir <fasta_flat> \
      --periods-dir <by_stem> --out-dir <det_out>
  python3 tools/detect_eval/eval.py ...

Output is a folder containing one `index.html` and an `assets/` dir
holding per-case line-width PNGs. Open `index.html` in any browser.

The dashboard has three linked panels:

  1. Aggregate (per-category accuracy bars + confusion matrix).
  2. Filterable per-case table.
  3. Case detail (periodogram + line-width raster + side-by-side
     kite / detect / truth rows), populated on row-click.

Design contract: docs/reports/kite_v2_dashboard/ is the canonical
output path for the v2 corpus; pass any --out-dir for other runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

# Plotly + jinja2 + PIL are required; pandas / matplotlib intentionally
# not used so the tool runs in minimal envs.
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Template
from PIL import Image

# Same category → expected_class mapping as eval.py — keep in sync.
CATEGORY_TO_EXPECTED = {
    "simple_tr":             "simple_TR",
    "hor_clean":             "HOR",
    "hor_wobble":            "HOR",
    "hor_shift":             "HOR",
    "hor_insertion":         "HOR",
    "hor_event_hybrid":      "HOR",
    "hor_event_inversion":   "HOR",
    "hor_event_duplication": "HOR",
    "hor_event_deletion":    "HOR",
    "mixed":                 "mixed",
    "random":                "ambiguous",
    "gc_bias":               "simple_TR",
}

CLASS_LABELS = ["simple_TR", "HOR", "irregular_HOR", "mixed", "ambiguous"]

# Base palette matches `src/detect/viz.rs::write_raster_png`.
BASE_COLOURS = {
    b"A": (0x2c, 0xa0, 0x2c),
    b"C": (0x1f, 0x77, 0xb4),
    b"G": (0xd6, 0x27, 0x28),
    b"T": (0x94, 0x67, 0xbd),
}
N_COLOUR = (0x99, 0x99, 0x99)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_kite_predictions(path: Path) -> dict[str, dict[str, str]]:
    """`kite-periodicity -o <out>.tsv` → dict by case_id."""
    out: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        out[row["case_id"]] = row
    return out


def load_periods_by_case(path: Path) -> dict[str, list[dict[str, str]]]:
    """`--emit-periods` output: rows grouped by `array_id`."""
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(path):
        out[row["array_id"]].append(row)
    return out


def load_detect_properties(properties_dir: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for p in properties_dir.rglob("*.properties.tsv"):
        for row in read_tsv(p):
            out[row["array_id"]] = row
    return out


def load_truth(corpus_root: Path) -> dict[str, dict[str, str]]:
    """Collects per-case `*.truth.tsv` rows by array_id."""
    out: dict[str, dict[str, str]] = {}
    for p in corpus_root.rglob("*.truth.tsv"):
        for row in read_tsv(p):
            out[row["array_id"]] = row
    return out


def load_fasta_record(fasta_path: Path) -> tuple[str, bytes]:
    """Read a single-record FASTA; return (id, normalized bytes)."""
    with fasta_path.open("rb") as f:
        lines = f.read().splitlines()
    header = b""
    seq_parts: list[bytes] = []
    for line in lines:
        if line.startswith(b">"):
            header = line[1:].split()[0]
        else:
            seq_parts.append(line.strip().upper())
    seq = b"".join(seq_parts)
    # Normalize non-ACGT to N.
    norm = bytes(b if b in b"ACGT" else ord("N") for b in seq)
    return header.decode("utf-8"), norm


# ---------------------------------------------------------------------------
# Line-width raster
# ---------------------------------------------------------------------------

def pick_raster_width(case: dict[str, Any]) -> Optional[int]:
    """Choose a width to wrap at for the line-width thumbnail.

    Priority:
      1. detector's chosen `base_width_bp` (if resolved class).
      2. kite founder (rule = HOR).
      3. kite tile (rule = HOR or Tandem).
      4. truth.tsv `base_width_bp` (fallback for ambiguous detection
         cases where we still want a visualization).
    Returns None if no usable width was found.
    """
    candidates = [
        case.get("detect_base_width"),
        case.get("kite_founder"),
        case.get("kite_tile"),
        case.get("truth_base_width"),
    ]
    for w in candidates:
        try:
            wi = int(w)
            if wi > 0:
                return wi
        except (TypeError, ValueError):
            continue
    return None


def write_line_width_png(
    seq: bytes,
    width: int,
    out_path: Path,
    max_rows: int = 256,
    max_cols: int = 600,
) -> None:
    """Render a wrapped-sequence colour heatmap. Downsampled if huge."""
    n_rows = len(seq) // width
    if n_rows == 0 or width == 0:
        return
    # Subsample rows + columns to keep the PNG small.
    row_stride = max(1, (n_rows + max_rows - 1) // max_rows)
    col_stride = max(1, (width + max_cols - 1) // max_cols)
    img_rows = (n_rows + row_stride - 1) // row_stride
    img_cols = (width + col_stride - 1) // col_stride
    img = Image.new("RGB", (img_cols, img_rows))
    px = img.load()
    for r_dst in range(img_rows):
        src_r = r_dst * row_stride
        base = src_r * width
        for c_dst in range(img_cols):
            src_c = c_dst * col_stride
            b = seq[base + src_c : base + src_c + 1]
            px[c_dst, r_dst] = BASE_COLOURS.get(b, N_COLOUR)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_path, optimize=True)


# ---------------------------------------------------------------------------
# Aggregate building blocks
# ---------------------------------------------------------------------------

def build_cases(
    manifest: list[dict[str, str]],
    kite_pred: dict[str, dict[str, str]],
    periods: dict[str, list[dict[str, str]]],
    detect_props: dict[str, dict[str, str]],
    truth: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """Build one case dict per manifest row. Missing inputs → None values."""
    cases: list[dict[str, Any]] = []
    for m in manifest:
        cid = m["case_id"]
        cat = m["category"]
        expected = CATEGORY_TO_EXPECTED.get(cat)
        kp = kite_pred.get(cid, {})
        dp = detect_props.get(cid, {})
        tr = truth.get(cid, {})
        per = periods.get(cid, [])
        c = {
            "case_id": cid,
            "category": cat,
            "expected": expected or "unknown",
            "got": dp.get("class", "missing"),
            "correct": (dp.get("class") == expected) if expected else False,
            # detect summary
            "detect_base_width": _opt_str(dp.get("base_width_bp")),
            "detect_hor_k": _opt_str(dp.get("hor_k")),
            "detect_conf": _opt_float(dp.get("confidence")),
            "detect_reason": dp.get("reason", ""),
            "detect_phase_sep": _opt_float(dp.get("phase_separation")),
            "detect_column_ic": _opt_float(dp.get("column_conservation")),
            "detect_irregularity": _opt_float(dp.get("irregularity_score")),
            "detect_n_phase_shifts": _opt_str(dp.get("n_phase_shifts")),
            # kite summary
            "kite_verdict": kp.get("verdict", "—"),
            "kite_founder": _opt_str(kp.get("founder")),
            "kite_multiplicity": _opt_str(kp.get("multiplicity")),
            "kite_tile": _opt_str(kp.get("tile")),
            "kite_share": _opt_float(kp.get("share")),
            "kite_n_peaks": _opt_str(kp.get("n_peaks_kept")),
            # truth
            "truth_class": tr.get("truth_class", "—"),
            "truth_base_width": _opt_str(tr.get("base_width_bp")),
            "truth_hor_k": _opt_str(tr.get("hor_k")),
            # period rows (for periodogram inline)
            "periods": [
                {
                    "period_bp": int(p["period_bp"]),
                    "period_score": float(p["period_score"]),
                    "source": p["source"],
                }
                for p in per
            ],
            "line_width_png": f"assets/line_width/{cid}.png",
        }
        cases.append(c)
    return cases


def _opt_str(v: Optional[str]) -> Optional[str]:
    if v is None or v == "" or v == "NA":
        return None
    return v


def _opt_float(v: Optional[str]) -> Optional[float]:
    if v is None or v == "" or v == "NA":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def aggregate_per_category(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_cat: dict[str, dict[str, int]] = defaultdict(lambda: {"correct": 0, "n": 0})
    for c in cases:
        by_cat[c["category"]]["n"] += 1
        if c["correct"]:
            by_cat[c["category"]]["correct"] += 1
    out = []
    for cat in sorted(by_cat):
        v = by_cat[cat]
        out.append({
            "category": cat,
            "correct": v["correct"],
            "n": v["n"],
            "pct": 100.0 * v["correct"] / v["n"] if v["n"] else 0.0,
        })
    return out


def aggregate_confusion(cases: list[dict[str, Any]]) -> list[list[int]]:
    """Confusion matrix indexed by CLASS_LABELS (rows=expected, cols=got).
    Last row holds "other expected" (e.g., random→ambiguous).
    """
    cm: dict[tuple[str, str], int] = defaultdict(int)
    for c in cases:
        cm[(c["expected"], c["got"])] += 1
    rows: list[list[int]] = []
    for exp in CLASS_LABELS:
        row = [cm.get((exp, got), 0) for got in CLASS_LABELS]
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Plotly figures
# ---------------------------------------------------------------------------

def fig_per_category(per_cat: list[dict[str, Any]]) -> str:
    cats = [r["category"] for r in per_cat]
    pcts = [r["pct"] for r in per_cat]
    text = [f"{r['correct']}/{r['n']}" for r in per_cat]
    fig = go.Figure(
        go.Bar(
            x=pcts, y=cats, orientation="h", text=text, textposition="inside",
            marker={"color": ["#2c8a3a" if p >= 90 else "#d97b00" if p >= 70 else "#c0392b" for p in pcts]},
            hovertemplate="%{y}<br>%{text} = %{x:.1f}%<extra></extra>",
        )
    )
    fig.update_layout(
        title="Per-category accuracy",
        xaxis={"range": [0, 100], "title": "% correct"},
        yaxis={"autorange": "reversed"},
        margin={"t": 40, "l": 130, "r": 20, "b": 40},
        height=380,
    )
    return pio.to_html(fig, include_plotlyjs=False, full_html=False, div_id="fig-per-cat")


def fig_confusion(cm: list[list[int]]) -> str:
    fig = go.Figure(
        go.Heatmap(
            z=cm, x=CLASS_LABELS, y=CLASS_LABELS,
            text=cm, texttemplate="%{text}",
            colorscale="Blues",
            hovertemplate="expected=%{y}<br>got=%{x}<br>n=%{z}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Confusion matrix (expected × got)",
        xaxis={"title": "detector's class", "side": "bottom"},
        yaxis={"title": "expected (oracle)", "autorange": "reversed"},
        margin={"t": 40, "l": 100, "r": 20, "b": 60},
        height=380,
    )
    return pio.to_html(fig, include_plotlyjs=False, full_html=False, div_id="fig-confusion")


def fig_kite_verdict_stack(cases: list[dict[str, Any]]) -> str:
    """Stacked bar: rule-classifier verdict distribution per category."""
    cats_set = sorted({c["category"] for c in cases})
    verdicts = ["hor", "tandem", "unresolved", "no_signal", "—"]
    counts: dict[tuple[str, str], int] = defaultdict(int)
    for c in cases:
        v = c.get("kite_verdict") or "—"
        counts[(c["category"], v)] += 1
    fig = go.Figure()
    for v in verdicts:
        ys = [counts.get((cat, v), 0) for cat in cats_set]
        fig.add_bar(name=v, x=cats_set, y=ys)
    fig.update_layout(
        barmode="stack",
        title="Kite rule verdict per category",
        xaxis={"title": "category"},
        yaxis={"title": "n cases"},
        margin={"t": 40, "l": 60, "r": 20, "b": 80},
        height=380,
        legend={"orientation": "h", "y": -0.2},
    )
    return pio.to_html(fig, include_plotlyjs=False, full_html=False, div_id="fig-kite-verdict")


# ---------------------------------------------------------------------------
# HTML template (inline)
# ---------------------------------------------------------------------------

TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         margin: 0; padding: 0 24px 60px; color: #222; }
  h1 { font-weight: 600; }
  .headline { display: flex; gap: 24px; align-items: baseline; padding: 16px 0; }
  .headline .big { font-size: 28px; font-weight: 600; }
  .headline .muted { color: #888; }
  .figs { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  .filters { display: flex; gap: 12px; padding: 12px 0; align-items: center; flex-wrap: wrap; }
  .filters label { font-size: 13px; }
  .filters select { padding: 4px 6px; font-size: 13px; }
  .filters button { padding: 4px 10px; }
  table.cases { border-collapse: collapse; width: 100%; font-size: 12px; }
  table.cases th, table.cases td { border-bottom: 1px solid #eee; padding: 4px 6px; text-align: left; }
  table.cases th { cursor: pointer; background: #fafafa; position: sticky; top: 0; }
  table.cases tr.row:hover { background: #f0f6ff; cursor: pointer; }
  table.cases tr.row.selected { background: #cfe3ff; }
  table.cases td.ok { color: #2c8a3a; }
  table.cases td.bad { color: #c0392b; font-weight: 600; }
  .pager { padding: 8px 0; font-size: 13px; }
  .pager button { margin: 0 4px; padding: 2px 8px; }
  .detail { border: 1px solid #ddd; margin-top: 16px; padding: 16px; background: #fafafa; }
  .detail h2 { margin-top: 0; }
  .detail .row { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  .detail img { max-width: 100%; height: auto; border: 1px solid #ccc; background: #fff; }
  .detail .raster-wrap { max-height: 380px; overflow: auto; }
  .detail .reason { background: #f4f4f4; padding: 8px; font-family: ui-monospace, monospace; font-size: 12px; white-space: pre-wrap; }
  .detail dl { display: grid; grid-template-columns: 140px 1fr; gap: 4px 12px; font-size: 13px; margin: 0; }
  .detail dt { color: #666; }
  .detail dd { margin: 0; font-family: ui-monospace, monospace; }
  .pill { display: inline-block; padding: 1px 8px; border-radius: 10px; font-size: 11px; color: #fff; }
  .pill.HOR { background: #2c8a3a; }
  .pill.simple_TR { background: #1f77b4; }
  .pill.irregular_HOR { background: #d97b00; }
  .pill.mixed { background: #6c3483; }
  .pill.ambiguous { background: #888; }
  .pill.missing { background: #c0392b; }
</style>
</head>
<body>

<h1>{{ title }}</h1>
<div class="headline">
  <div><span class="big">{{ overall_pct }}%</span>
       <span class="muted">overall ({{ overall_correct }}/{{ overall_n }})</span></div>
  {% if oracle_pct is not none %}
  <div class="muted">oracle baseline: {{ oracle_pct }}%</div>
  {% endif %}
  <div class="muted">commit {{ commit }} • {{ generated_at }}</div>
</div>

<div class="figs">
  {{ fig_per_cat | safe }}
  {{ fig_confusion | safe }}
</div>
<div>
  {{ fig_kite_verdict | safe }}
</div>

<h2>Cases ({{ n_cases }} total)</h2>
<div class="filters">
  <label>category: <select id="f-cat"><option value="">all</option>
    {% for cat in categories %}<option>{{ cat }}</option>{% endfor %}
  </select></label>
  <label>expected: <select id="f-exp"><option value="">all</option>
    {% for c in CLASS_LABELS %}<option>{{ c }}</option>{% endfor %}
  </select></label>
  <label>detector got: <select id="f-got"><option value="">all</option>
    {% for c in CLASS_LABELS %}<option>{{ c }}</option>{% endfor %}
  </select></label>
  <label>kite verdict: <select id="f-kv"><option value="">all</option>
    <option>hor</option><option>tandem</option>
    <option>unresolved</option><option>no_signal</option>
  </select></label>
  <label><input type="checkbox" id="f-wrong"> wrong only</label>
  <button onclick="resetFilters()">reset</button>
  <span id="row-count" class="muted"></span>
</div>

<div class="pager">
  <button onclick="page(-1)">◀ prev</button>
  <span id="page-label">—</span>
  <button onclick="page(1)">next ▶</button>
  <span class="muted">(page size {{ page_size }})</span>
</div>

<table class="cases" id="cases-table">
  <thead><tr>
    <th data-key="case_id">case_id</th>
    <th data-key="category">category</th>
    <th data-key="expected">expected</th>
    <th data-key="got">got</th>
    <th data-key="correct">✓</th>
    <th data-key="kite_verdict">kite verdict</th>
    <th data-key="kite_founder">founder</th>
    <th data-key="kite_tile">tile</th>
    <th data-key="detect_base_width">base_w</th>
    <th data-key="detect_hor_k">k</th>
    <th data-key="detect_conf">conf</th>
  </tr></thead>
  <tbody></tbody>
</table>

<div class="detail" id="detail" style="display:none">
  <h2 id="d-case-id"></h2>
  <div class="row">
    <div>
      <dl>
        <dt>category</dt><dd id="d-category"></dd>
        <dt>expected</dt><dd><span id="d-expected" class="pill"></span></dd>
        <dt>detect got</dt><dd><span id="d-got" class="pill"></span> <span id="d-correct"></span></dd>
        <dt>truth class</dt><dd id="d-truth-class"></dd>
        <dt>truth (k, base)</dt><dd id="d-truth-k-base"></dd>
        <dt>detect (k, base)</dt><dd id="d-det-k-base"></dd>
        <dt>confidence</dt><dd id="d-conf"></dd>
        <dt>kite verdict</dt><dd id="d-kv"></dd>
        <dt>kite (founder, k, tile)</dt><dd id="d-kfk"></dd>
        <dt>kite share / n_peaks</dt><dd id="d-kshare"></dd>
        <dt>phase_sep / col_IC</dt><dd id="d-phase-ic"></dd>
        <dt>irreg / n_shifts</dt><dd id="d-irreg-shifts"></dd>
      </dl>
      <h3>Detector reason</h3>
      <div class="reason" id="d-reason"></div>
    </div>
    <div>
      <h3>Periodogram (kite peaks emitted to detector)</h3>
      <div id="d-periodogram" style="height:300px"></div>
      <h3>Line-width raster <span class="muted" id="d-rast-w" style="font-size:12px"></span></h3>
      <div class="raster-wrap"><img id="d-rast-img" alt="line-width"></div>
    </div>
  </div>
</div>

<script>
const CASES = {{ cases_json | safe }};
const PAGE_SIZE = {{ page_size }};
let filtered = CASES.slice();
let pageIdx = 0;
let sortKey = "case_id"; let sortAsc = true;
let selectedId = null;

function applyFilters() {
  const cat = document.getElementById("f-cat").value;
  const exp = document.getElementById("f-exp").value;
  const got = document.getElementById("f-got").value;
  const kv  = document.getElementById("f-kv").value;
  const wrong = document.getElementById("f-wrong").checked;
  filtered = CASES.filter(c =>
    (!cat || c.category === cat) &&
    (!exp || c.expected === exp) &&
    (!got || c.got === got) &&
    (!kv  || c.kite_verdict === kv) &&
    (!wrong || !c.correct)
  );
  sortRows();
  pageIdx = 0;
  render();
}

function sortRows() {
  filtered.sort((a, b) => {
    let av = a[sortKey], bv = b[sortKey];
    if (av === null || av === undefined) av = "";
    if (bv === null || bv === undefined) bv = "";
    if (typeof av === "number" || typeof bv === "number") {
      av = Number(av) || -Infinity; bv = Number(bv) || -Infinity;
      return sortAsc ? av - bv : bv - av;
    }
    av = String(av); bv = String(bv);
    return sortAsc ? av.localeCompare(bv) : bv.localeCompare(av);
  });
}

function render() {
  const tbody = document.querySelector("#cases-table tbody");
  tbody.innerHTML = "";
  const start = pageIdx * PAGE_SIZE;
  const end = Math.min(start + PAGE_SIZE, filtered.length);
  for (let i = start; i < end; i++) {
    const c = filtered[i];
    const tr = document.createElement("tr");
    tr.className = "row" + (c.case_id === selectedId ? " selected" : "");
    tr.dataset.id = c.case_id;
    const correctCls = c.correct ? "ok" : "bad";
    tr.innerHTML = `
      <td>${c.case_id}</td>
      <td>${c.category}</td>
      <td>${c.expected}</td>
      <td>${c.got}</td>
      <td class="${correctCls}">${c.correct ? "✓" : "✗"}</td>
      <td>${c.kite_verdict}</td>
      <td>${c.kite_founder ?? ""}</td>
      <td>${c.kite_tile ?? ""}</td>
      <td>${c.detect_base_width ?? ""}</td>
      <td>${c.detect_hor_k ?? ""}</td>
      <td>${c.detect_conf !== null && c.detect_conf !== undefined ? c.detect_conf.toFixed(3) : ""}</td>
    `;
    tr.addEventListener("click", () => selectCase(c.case_id));
    tbody.appendChild(tr);
  }
  document.getElementById("row-count").textContent =
    `${filtered.length} matching cases`;
  document.getElementById("page-label").textContent =
    `page ${pageIdx + 1} / ${Math.max(1, Math.ceil(filtered.length / PAGE_SIZE))}`;
}

function page(d) {
  const max = Math.ceil(filtered.length / PAGE_SIZE) - 1;
  pageIdx = Math.max(0, Math.min(max, pageIdx + d));
  render();
}

function resetFilters() {
  ["f-cat","f-exp","f-got","f-kv"].forEach(id => document.getElementById(id).value = "");
  document.getElementById("f-wrong").checked = false;
  applyFilters();
}

function selectCase(id) {
  selectedId = id;
  const c = CASES.find(x => x.case_id === id);
  if (!c) return;
  const set = (k, v) => { document.getElementById(k).textContent = (v ?? ""); };
  set("d-case-id", c.case_id);
  set("d-category", c.category);
  const exp = document.getElementById("d-expected");
  exp.className = "pill " + c.expected; exp.textContent = c.expected;
  const got = document.getElementById("d-got");
  got.className = "pill " + c.got; got.textContent = c.got;
  set("d-correct", c.correct ? "✓ correct" : "✗ wrong");
  set("d-truth-class", c.truth_class);
  set("d-truth-k-base", `k=${c.truth_hor_k ?? "—"}  base=${c.truth_base_width ?? "—"}`);
  set("d-det-k-base", `k=${c.detect_hor_k ?? "—"}  base=${c.detect_base_width ?? "—"}`);
  set("d-conf", c.detect_conf !== null && c.detect_conf !== undefined ? c.detect_conf.toFixed(4) : "—");
  set("d-kv", c.kite_verdict);
  set("d-kfk", `founder=${c.kite_founder ?? "—"}  k=${c.kite_multiplicity ?? "—"}  tile=${c.kite_tile ?? "—"}`);
  set("d-kshare", `share=${c.kite_share !== null && c.kite_share !== undefined ? c.kite_share.toFixed(3) : "—"}  n_peaks=${c.kite_n_peaks ?? "—"}`);
  set("d-phase-ic", `phase_sep=${c.detect_phase_sep !== null ? (c.detect_phase_sep ?? "").toFixed?.(4) ?? "—" : "—"}  IC=${c.detect_column_ic !== null ? (c.detect_column_ic ?? "").toFixed?.(4) ?? "—" : "—"}`);
  set("d-irreg-shifts", `irreg=${c.detect_irregularity !== null ? (c.detect_irregularity ?? "").toFixed?.(3) ?? "—" : "—"}  n_shifts=${c.detect_n_phase_shifts ?? "—"}`);
  set("d-reason", c.detect_reason);
  // Periodogram plot.
  const trace = {
    x: c.periods.map(p => p.period_bp),
    y: c.periods.map(p => p.period_score),
    text: c.periods.map(p => p.source),
    mode: "markers",
    type: "scatter",
    marker: {
      size: 12,
      color: c.periods.map(p => SOURCE_COLOUR[p.source] ?? "#888")
    },
    hovertemplate: "period=%{x} bp<br>score=%{y:.3f}<br>%{text}<extra></extra>"
  };
  Plotly.newPlot("d-periodogram", [trace], {
    margin: {t: 10, l: 50, r: 10, b: 50},
    xaxis: {title: "period (bp)"},
    yaxis: {title: "period_score", range: [0, 1]},
    shapes: [{
      type: "line", x0: 0, x1: 1, xref: "paper", y0: 0.85, y1: 0.85,
      line: {color: "#c0392b", width: 1, dash: "dash"}
    }],
    annotations: [{
      x: 0.99, xref: "paper", y: 0.87, yref: "y", xanchor: "right",
      text: "strong_period_score = 0.85", showarrow: false,
      font: {color: "#c0392b", size: 10}
    }]
  }, {displayModeBar: false});
  // Raster image.
  document.getElementById("d-rast-img").src = c.line_width_png;
  document.getElementById("d-rast-w").textContent = "(rendered at width=" + (c.detect_base_width ?? c.kite_founder ?? c.kite_tile ?? c.truth_base_width ?? "?") + " bp)";
  document.getElementById("detail").style.display = "";
  // Highlight the clicked row.
  document.querySelectorAll("#cases-table tbody tr.row").forEach(tr => {
    tr.classList.toggle("selected", tr.dataset.id === id);
  });
  document.getElementById("detail").scrollIntoView({behavior: "smooth", block: "start"});
}

const SOURCE_COLOUR = {
  kite_founder: "#2c8a3a",
  kite_tile: "#1f77b4",
  kite_monomer: "#2c8a3a",
  kite_secondary: "#888",
  kite_peak: "#d97b00"
};

document.querySelectorAll("#cases-table th").forEach(th => {
  th.addEventListener("click", () => {
    const k = th.dataset.key;
    if (sortKey === k) sortAsc = !sortAsc;
    else { sortKey = k; sortAsc = true; }
    sortRows();
    pageIdx = 0;
    render();
  });
});
["f-cat","f-exp","f-got","f-kv","f-wrong"].forEach(id =>
  document.getElementById(id).addEventListener("change", applyFilters)
);

applyFilters();
</script>

</body>
</html>
"""


def render_dashboard(
    out_dir: Path,
    title: str,
    cases: list[dict[str, Any]],
    oracle_pct: Optional[float],
    commit: str,
    generated_at: str,
    page_size: int = 200,
) -> None:
    per_cat = aggregate_per_category(cases)
    cm = aggregate_confusion(cases)
    overall_correct = sum(1 for c in cases if c["correct"])
    overall_n = len(cases)
    overall_pct = round(100.0 * overall_correct / overall_n, 1) if overall_n else 0.0
    figs = {
        "fig_per_cat": fig_per_category(per_cat),
        "fig_confusion": fig_confusion(cm),
        "fig_kite_verdict": fig_kite_verdict_stack(cases),
    }
    html = Template(TEMPLATE).render(
        title=title,
        cases_json=json.dumps(cases, default=str),
        categories=sorted({c["category"] for c in cases}),
        CLASS_LABELS=CLASS_LABELS,
        overall_correct=overall_correct,
        overall_n=overall_n,
        overall_pct=overall_pct,
        oracle_pct=oracle_pct,
        commit=commit,
        generated_at=generated_at,
        n_cases=overall_n,
        page_size=page_size,
        **figs,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "index.html").write_text(html)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    import datetime
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--kite", type=Path, required=True,
                    help="Predictions TSV from `kite-periodicity -o ...`.")
    ap.add_argument("--periods", type=Path, required=True,
                    help="`--emit-periods` output (v2 detector schema).")
    ap.add_argument("--properties-dir", type=Path, required=True,
                    help="Detector batch output directory.")
    ap.add_argument("--truth-root", type=Path, required=True,
                    help="Corpus root containing per-case `*.truth.tsv`.")
    ap.add_argument("--fasta-dir", type=Path, required=True,
                    help="Flat dir of per-case FASTAs (symlinks fine).")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--title", type=str, default="kite → detect — v2 corpus")
    ap.add_argument("--commit", type=str, default="HEAD",
                    help="Commit hash to embed in the report header.")
    ap.add_argument("--oracle-pct", type=float, default=None,
                    help="Oracle baseline (e.g., 94.4) shown in the header.")
    ap.add_argument("--skip-png", action="store_true",
                    help="Skip line-width PNG rendering (fast iteration).")
    ap.add_argument("--page-size", type=int, default=200)
    args = ap.parse_args()

    print("loading inputs...", file=sys.stderr)
    manifest = read_tsv(args.manifest)
    kite_pred = load_kite_predictions(args.kite)
    periods = load_periods_by_case(args.periods)
    detect_props = load_detect_properties(args.properties_dir)
    truth = load_truth(args.truth_root)
    cases = build_cases(manifest, kite_pred, periods, detect_props, truth)
    print(f"  manifest: {len(manifest)}", file=sys.stderr)
    print(f"  kite preds: {len(kite_pred)}", file=sys.stderr)
    print(f"  detect props: {len(detect_props)}", file=sys.stderr)
    print(f"  truth rows: {len(truth)}", file=sys.stderr)
    print(f"  cases assembled: {len(cases)}", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    asset_root = args.out_dir / "assets" / "line_width"
    asset_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_png:
        print(f"rendering line-width PNGs into {asset_root}...", file=sys.stderr)
        for i, c in enumerate(cases):
            cid = c["case_id"]
            fa = args.fasta_dir / f"{cid}.fa"
            if not fa.exists():
                continue
            w = pick_raster_width(c)
            if not w:
                continue
            _hdr, seq = load_fasta_record(fa)
            write_line_width_png(seq, w, asset_root / f"{cid}.png")
            if (i + 1) % 200 == 0:
                print(f"  {i + 1}/{len(cases)}", file=sys.stderr)

    print(f"rendering index.html into {args.out_dir}...", file=sys.stderr)
    render_dashboard(
        out_dir=args.out_dir,
        title=args.title,
        cases=cases,
        oracle_pct=args.oracle_pct,
        commit=args.commit,
        generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        page_size=args.page_size,
    )
    print(f"done → open {args.out_dir / 'index.html'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Pick putative-subrepeat candidates from a rescore TSV and render
dottir dotplots for each.

Selection criteria (rows must satisfy ALL):
  * founder_period known (non-NA)
  * period in [--period-min, --period-max]
  * period <= founder_period * --ratio-max  (nested-period gate)
  * scan_occupancy_frac in [--occ-min, --occ-max]
  * scan_n_intervals >= --min-intervals
  * array_length >= 6 * founder_period (so we can skip 2 + plot N more)

For each surviving (record, period, founder) triple:
  1. Extract a sub-region of the array from offset
     `--skip-founders * founder` for `--plot-founders * founder` bp.
  2. Write a single-record FASTA.
  3. Invoke `dottir batch` on it (self-comparison) with the requested
     window / pixel-fac / width.
  4. Save the PNG (and dottir's .params.toml sidecar) under <out-dir>.

Output filename pattern:
    P{period}_F{founder}_occ{occ4d}_{safe_recordid}.png

Diversity cap: at most `--per-family` candidates per TRC family.

Usage:
    python3 tools/subrepeat_scan/render_candidates.py \\
        --fasta test_data/IPIP200579_2026-04-14_combined.fasta \\
        --rescore /tmp/rescore_v012/ipip_scan.peaks.tsv \\
        --dottir /tmp/dottir-src/target/release/dottir \\
        --out-dir /tmp/dotplots_subrepeats
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def family(record_id: str) -> str:
    # TRC_104:chr3_411509737_411516552 → TRC_104
    return record_id.split(":")[0]


def pick_candidates(
    rescore_tsv: Path,
    *,
    period_min: int,
    period_max: int,
    ratio_max: float,
    occ_min: float,
    occ_max: float,
    min_intervals: int,
    skip_founders: int,
    plot_founders: int,
    per_family: int,
) -> list[dict]:
    """Return a list of candidate dicts in selection order."""
    candidates: list[dict] = []
    with open(rescore_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                period = int(row["period"])
                founder = row.get("founder_period", "NA")
                if founder == "NA":
                    continue
                founder = int(founder)
                array_length = int(row["array_length"])
                occ = row.get("scan_occupancy_frac", "NA")
                n_int = row.get("scan_n_intervals", "NA")
                if occ == "NA" or n_int == "NA":
                    continue
                occ = float(occ)
                n_int = int(n_int)
            except (KeyError, ValueError):
                continue
            if not (period_min <= period <= period_max):
                continue
            if period > founder * ratio_max:
                continue
            if not (occ_min <= occ <= occ_max):
                continue
            if n_int < min_intervals:
                continue
            if array_length < (skip_founders + plot_founders) * founder:
                continue
            candidates.append(
                {
                    "record_id": row["case_id"],
                    "period": period,
                    "founder": founder,
                    "occ": occ,
                    "n_int": n_int,
                    "array_length": array_length,
                }
            )

    # Sort by occupancy desc, then by n_intervals desc.
    candidates.sort(key=lambda c: (-c["occ"], -c["n_int"]))

    # Diversity cap: at most N per TRC family.
    seen: defaultdict[str, int] = defaultdict(int)
    kept: list[dict] = []
    for c in candidates:
        fam = family(c["record_id"])
        if seen[fam] >= per_family:
            continue
        seen[fam] += 1
        kept.append(c)
    return kept


def iter_fasta(path: Path):
    with open(path) as f:
        header = None
        buf: list[bytes] = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, b"".join(buf)
                header = line[1:].split()[0]
                buf = []
            elif header is not None:
                buf.append(line.encode("ascii"))
        if header is not None:
            yield header, b"".join(buf)


def index_fasta(path: Path) -> dict[str, bytes]:
    return {rid: seq for rid, seq in iter_fasta(path)}


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


def write_subregion(
    seq: bytes,
    record_id: str,
    start: int,
    end: int,
    out_path: Path,
) -> None:
    sub = seq[start:end]
    with open(out_path, "w") as f:
        f.write(f">{record_id}|{start}-{end}\n")
        for i in range(0, len(sub), 80):
            f.write(sub[i : i + 80].decode("ascii") + "\n")


def render_one(
    dottir: Path,
    fa_path: Path,
    out_png: Path,
    *,
    window: int,
    pixel_fac: float,
    width: int,
) -> tuple[bool, str]:
    cmd = [
        str(dottir),
        "batch",
        str(fa_path),
        "-o",
        str(out_png),
        "-W",
        str(window),
        "--pixel-fac",
        str(pixel_fac),
        "--width",
        str(width),
        "--no-sidecar",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return False, "timeout"
    if r.returncode != 0:
        return False, r.stderr.strip().splitlines()[-1] if r.stderr else "nonzero exit"
    return True, ""


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--rescore", type=Path, required=True,
                   help="rescore peaks TSV (with scan_occupancy_frac column)")
    p.add_argument("--dottir", type=Path, required=True,
                   help="path to the dottir CLI binary")
    p.add_argument("--out-dir", type=Path, required=True)

    # Selection criteria
    p.add_argument("--period-min", type=int, default=20)
    p.add_argument("--period-max", type=int, default=200)
    p.add_argument("--ratio-max", type=float, default=0.25,
                   help="max period/founder_period")
    p.add_argument("--occ-min", type=float, default=0.10)
    p.add_argument("--occ-max", type=float, default=0.80)
    p.add_argument("--min-intervals", type=int, default=5)
    p.add_argument("--per-family", type=int, default=2,
                   help="max candidates per TRC family")
    p.add_argument("--top", type=int, default=25,
                   help="cap total number of candidates rendered")

    # Sub-region geometry
    p.add_argument("--skip-founders", type=int, default=2,
                   help="skip the first N founder copies (boundary irregularity)")
    p.add_argument("--plot-founders", type=int, default=5,
                   help="render this many founder copies")

    # dottir rendering knobs
    p.add_argument("--window", type=int, default=15,
                   help="dottir window size -W")
    p.add_argument("--pixel-fac", type=int, default=0,
                   help="dottir --pixel-fac integer (0 = auto)")
    p.add_argument("--png-width", type=int, default=1500,
                   help="dottir --width in pixels")

    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning {args.rescore} for candidates...", file=sys.stderr)
    candidates = pick_candidates(
        args.rescore,
        period_min=args.period_min,
        period_max=args.period_max,
        ratio_max=args.ratio_max,
        occ_min=args.occ_min,
        occ_max=args.occ_max,
        min_intervals=args.min_intervals,
        skip_founders=args.skip_founders,
        plot_founders=args.plot_founders,
        per_family=args.per_family,
    )
    if args.top > 0:
        candidates = candidates[: args.top]
    print(f"Selected {len(candidates)} candidate(s)", file=sys.stderr)
    if not candidates:
        return

    print(f"Indexing FASTA {args.fasta} ...", file=sys.stderr)
    needed = {c["record_id"] for c in candidates}
    seqs: dict[str, bytes] = {}
    for rid, seq in iter_fasta(args.fasta):
        if rid in needed:
            seqs[rid] = seq
            if len(seqs) == len(needed):
                break

    missing = needed - seqs.keys()
    if missing:
        print(f"warning: {len(missing)} record(s) missing from FASTA", file=sys.stderr)

    manifest_rows = []
    n_ok = 0
    n_fail = 0
    for c in candidates:
        rid = c["record_id"]
        if rid not in seqs:
            continue
        seq = seqs[rid]
        founder = c["founder"]
        start = args.skip_founders * founder
        end = start + args.plot_founders * founder
        end = min(end, len(seq))
        if end - start < 3 * founder:
            continue

        safe_rid = safe_name(rid)
        fa_path = args.out_dir / f"P{c['period']}_F{founder}_{safe_rid}.fa"
        png_path = args.out_dir / (
            f"P{c['period']:03d}_F{founder:04d}_occ{int(c['occ'] * 10000):04d}_"
            f"{safe_rid}.png"
        )
        write_subregion(seq, rid, start, end, fa_path)
        ok, err = render_one(
            args.dottir,
            fa_path,
            png_path,
            window=args.window,
            pixel_fac=args.pixel_fac,
            width=args.png_width,
        )
        # Drop the intermediate sub-FASTA after successful render.
        if ok:
            n_ok += 1
            fa_path.unlink(missing_ok=True)
        else:
            n_fail += 1
            print(f"  FAIL {png_path.name}: {err}", file=sys.stderr)

        manifest_rows.append(
            {
                "record_id": rid,
                "period": c["period"],
                "founder": founder,
                "occupancy_frac": c["occ"],
                "n_intervals": c["n_int"],
                "array_length": c["array_length"],
                "region_start_bp": start,
                "region_end_bp": end,
                "png": png_path.name if ok else "(failed)",
            }
        )

    manifest_path = args.out_dir / "manifest.tsv"
    with open(manifest_path, "w") as f:
        if manifest_rows:
            headers = list(manifest_rows[0].keys())
            f.write("\t".join(headers) + "\n")
            for r in manifest_rows:
                f.write("\t".join(str(r[h]) for h in headers) + "\n")

    print(f"rendered {n_ok}, failed {n_fail}; manifest at {manifest_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Self-self k-mer dotplot for one FASTA record, rendered as a PNG.

Useful for visually inspecting nested tandem repeats — diagonals at
short offsets (subrepeats) are easy to see against the founder-scale
diagonal pattern.

The plot is a scatter of `(i, j)` coordinates where `seq[i:i+k] ==
seq[j:j+k]` and `j > i`. The main diagonal `(i, i)` is omitted by
construction.

Usage:
    python3 tools/subrepeat_scan/dotplot.py \\
        --fasta <fasta> --record <record_id> \\
        --out <out.png> \\
        [--k 12] [--start 0] [--end 10000] [--max-lag 2000]

Notes:
- `--k 12` is a sensible default for centromeric arrays (long enough
  to be specific, short enough to fire on real repetition).
- Very large windows (e.g. 100 kb × 100 kb) produce many points; use
  `--start` / `--end` to focus on a sub-region.
- `--max-lag` limits plotting to short-range repeats (helps see the
  subrepeat scale without losing the picture to the founder diagonal).
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_record(fasta_path: Path, record_id: str) -> bytes:
    """Return the sequence of `record_id` from `fasta_path`."""
    with open(fasta_path) as f:
        target = False
        buf = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                rid = line[1:].split()[0]
                if rid == record_id:
                    target = True
                elif target:
                    # finished collecting target
                    break
                continue
            if target:
                buf.append(line.encode("ascii"))
    if not buf:
        raise SystemExit(f"record {record_id!r} not found in {fasta_path}")
    return b"".join(buf)


def build_kmer_dots(seq: bytes, k: int, max_lag: int | None):
    """Return arrays (xs, ys) of dot coordinates."""
    table = {b"A": 0, b"C": 1, b"G": 2, b"T": 3,
             b"a": 0, b"c": 1, b"g": 2, b"t": 3}
    L = len(seq)
    positions = defaultdict(list)
    h = 0
    valid = 0
    mask = (1 << (2 * k)) - 1
    for i, bb in enumerate(seq):
        code = table.get(bytes([bb]))
        if code is None:
            h = 0
            valid = 0
            continue
        h = ((h << 2) | code) & mask
        valid += 1
        if valid >= k:
            positions[h].append(i + 1 - k)

    xs: list[int] = []
    ys: list[int] = []
    for poses in positions.values():
        if len(poses) < 2:
            continue
        for a_idx, a in enumerate(poses):
            for b in poses[a_idx + 1 :]:
                lag = b - a
                if max_lag is not None and lag > max_lag:
                    break  # poses is sorted ascending
                xs.append(a)
                ys.append(b)
    return np.array(xs, dtype=np.int64), np.array(ys, dtype=np.int64)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--record", required=True, help="record id to plot")
    p.add_argument("--out", type=Path, required=True,
                   help="output PNG path")
    p.add_argument("--k", type=int, default=12,
                   help="k-mer length (default 12)")
    p.add_argument("--start", type=int, default=0,
                   help="sub-region start (bp); default 0")
    p.add_argument("--end", type=int, default=0,
                   help="sub-region end (bp); 0 = whole record")
    p.add_argument("--max-lag", type=int, default=0,
                   help="cap maximum dot lag (j − i) in bp; 0 = no cap")
    p.add_argument("--canvas", type=int, default=1500,
                   help="output canvas size in pixels (square); the "
                        "region is binned into canvas×canvas cells")
    args = p.parse_args()

    seq = read_record(args.fasta, args.record)
    L = len(seq)
    print(f"loaded {args.record}: {L} bp", file=sys.stderr)
    s = max(0, args.start)
    e = args.end if args.end > 0 else L
    e = min(e, L)
    if e - s < args.k:
        sys.exit(f"sub-region too short ({e - s} bp) for k={args.k}")
    sub = seq[s:e]
    max_lag = args.max_lag if args.max_lag > 0 else None

    xs, ys = build_kmer_dots(sub, args.k, max_lag)
    print(f"sub-region {s}..{e} ({e - s} bp), k={args.k}, "
          f"max_lag={max_lag}, n_dots={len(xs)}", file=sys.stderr)

    # Render via PIL — bin dots into a fixed-size canvas (typical
    # regions are 5–100 kb, far more than reasonable pixel resolution).
    from PIL import Image, ImageDraw, ImageFont

    canvas = args.canvas
    if canvas <= 0:
        canvas = 1500
    region_len = e - s
    counts = np.zeros((canvas, canvas), dtype=np.int32)
    if len(xs) > 0:
        # Map each (x, y) to a bin
        bx = ((xs.astype(np.int64)) * canvas // region_len)
        by = ((ys.astype(np.int64)) * canvas // region_len)
        bx = np.clip(bx, 0, canvas - 1)
        by = np.clip(by, 0, canvas - 1)
        np.add.at(counts, (by, bx), 1)

    # Symmetrise so the plot looks like a true square dotplot
    counts = np.maximum(counts, counts.T)
    # Render as inverted black-on-white. Any bin with a hit at all is
    # dark; intensity scales with sqrt(count). The sqrt instead of
    # log keeps single-hit bins clearly visible.
    if counts.max() > 0:
        vis = np.sqrt(counts.astype(np.float32))
        vis = 255.0 - 220.0 * (vis / vis.max())  # 35 = blackest
    else:
        vis = np.full(counts.shape, 255.0, dtype=np.float32)
    vis_img = vis.clip(0, 255).astype(np.uint8)

    # Draw the main diagonal in light gray as a hint of orientation
    # (do it post-binning so we don't blow out the histogram)
    diag_mask = np.eye(canvas, dtype=bool)
    vis_img[diag_mask] = np.minimum(vis_img[diag_mask], 220)

    img = Image.fromarray(vis_img, mode="L").convert("RGB")
    draw = ImageDraw.Draw(img)
    title = (
        f"{args.record}  k={args.k}  "
        f"region={s}-{e} bp ({region_len} bp)"
        + (f"  max_lag={max_lag}" if max_lag else "")
    )
    try:
        font = ImageFont.load_default()
        draw.text((6, 4), title, fill=(180, 0, 0), font=font)
        draw.text(
            (6, canvas - 16),
            f"x = position (bp); axis bin size ≈ {region_len / canvas:.1f} bp",
            fill=(180, 0, 0),
            font=font,
        )
    except Exception:
        pass

    args.out.parent.mkdir(parents=True, exist_ok=True)
    img.save(args.out)
    print(f"wrote {args.out}  ({canvas}x{canvas} canvas)", file=sys.stderr)


if __name__ == "__main__":
    main()

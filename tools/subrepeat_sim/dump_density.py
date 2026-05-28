#!/usr/bin/env python3
"""Dump the per-window density profile of k-mer pairs at distance P
for one record in a FASTA. Used to inspect why kmer_autocorr_founder
underperforms on TRC_104 vs simulation.

Usage:
    python3 tools/subrepeat_sim/dump_density.py \\
        <fasta> <record_id> <period> <founder> <out.tsv>

Output TSV columns:
    window_start, window_end, midpoints_in_window, density_per_bp,
    phase_within_founder
"""
import sys
from pathlib import Path
from collections import defaultdict


def read_record(fasta_path: Path, record_id: str) -> bytes:
    with open(fasta_path) as f:
        target = None
        buf = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if target is not None:
                    return b"".join(buf)
                if line[1:].split()[0] == record_id:
                    target = record_id
                continue
            if target is not None:
                buf.append(line.encode("ascii"))
    if target is None:
        raise SystemExit(f"record {record_id!r} not found")
    return b"".join(buf)


def kmer_positions(seq: bytes, k: int) -> dict:
    table = {b"A": 0, b"C": 1, b"G": 2, b"T": 3,
             b"a": 0, b"c": 1, b"g": 2, b"t": 3}
    out = defaultdict(list)
    h = 0
    valid = 0
    mask = (1 << (2 * k)) - 1
    for i, b in enumerate(seq):
        bb = bytes([b])
        code = table.get(bb)
        if code is None:
            h = 0
            valid = 0
            continue
        h = ((h << 2) | code) & mask
        valid += 1
        if valid >= k:
            out[h].append(i + 1 - k)
    return out


def main():
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    fasta = Path(sys.argv[1])
    record_id = sys.argv[2]
    period = int(sys.argv[3])
    founder = int(sys.argv[4])
    out_path = Path(sys.argv[5])

    seq = read_record(fasta, record_id)
    print(f"Loaded {record_id}: {len(seq)} bp", file=sys.stderr)

    k = 6
    distance_tol = 3
    positions = kmer_positions(seq, k)
    print(f"Unique k-mers: {len(positions)}", file=sys.stderr)

    # Collect midpoints of consecutive-pair distances in [P-tol, P+tol]
    midpoints = []
    for poses in positions.values():
        if len(poses) < 2:
            continue
        for i in range(1, len(poses)):
            d = poses[i] - poses[i - 1]
            if abs(d - period) <= distance_tol:
                midpoints.append((poses[i - 1] + poses[i]) // 2)

    print(f"Matching pair midpoints at P={period}±{distance_tol}: "
          f"{len(midpoints)}", file=sys.stderr)

    # Bin into sliding windows of width = founder // 6, step = win // 2
    # (matches rescore's autocorr metric default).
    win = max(20, period) // 2  # same as Rust default
    step = max(1, win // 2)
    n_windows = (len(seq) - win) // step + 1
    counts = [0] * n_windows
    for m in midpoints:
        idx = max(0, (m - win // 2) // step)
        idx = min(n_windows - 1, idx)
        counts[idx] += 1

    print(f"Windows: n={n_windows} win={win} step={step}", file=sys.stderr)

    with open(out_path, "w") as f:
        f.write("window_idx\twindow_start\twindow_end\tcount\tphase_mod_founder\n")
        for i, c in enumerate(counts):
            start = i * step
            end = start + win
            center = (start + end) // 2
            phase = center % founder
            f.write(f"{i}\t{start}\t{end}\t{c}\t{phase}\n")

    # Also dump a "phase-folded" summary: for each phase bin (founder/12),
    # accumulate counts across all founders.
    phase_bins = 12
    phase_hist = [0] * phase_bins
    for m in midpoints:
        phase = m % founder
        b = min(phase_bins - 1, phase * phase_bins // founder)
        phase_hist[b] += 1

    summary_path = out_path.with_suffix(".phase_hist.tsv")
    with open(summary_path, "w") as f:
        f.write("phase_bin\tphase_start_bp\tphase_end_bp\tcount\tfraction\n")
        total = sum(phase_hist)
        for b, c in enumerate(phase_hist):
            pstart = b * founder // phase_bins
            pend = (b + 1) * founder // phase_bins
            frac = c / total if total else 0
            f.write(f"{b}\t{pstart}\t{pend}\t{c}\t{frac:.4f}\n")
    print(f"Wrote {out_path}", file=sys.stderr)
    print(f"Wrote {summary_path}", file=sys.stderr)
    print(f"Total midpoints folded into phase histogram: {total}", file=sys.stderr)


if __name__ == "__main__":
    main()

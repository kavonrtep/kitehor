#!/usr/bin/env python3
"""Shifted self-alignment scan for nested short tandem repeats.

For each long-monomer (founder) along a tandem-repeat array, find
contiguous regions where the sequence is similar to itself shifted by
some short period `p`. Each such region is a candidate nested
subrepeat — a short tandem motif tiling inside the founder.

The algorithm is a thin per-base autocorrelation:

    1. For each candidate period p in [period_min, period_max]:
         match_p(i) = 1 if seq[i] == seq[i+p] else 0
    2. Smooth match_p over a window of width max(3*p, min_window) bp.
    3. Threshold the smoothed signal at `--id-threshold`.
    4. Take contiguous runs ≥ min_copies * p as candidate intervals.

Period aliasing (a 30 bp repeat also lights up at 60, 90, …) is
resolved by sorting candidate intervals by period ascending and
greedily dropping any longer-period interval that overlaps a
shorter-period one by ≥ `--overlap-frac`.

Per-founder occupancy is computed by intersecting each surviving
interval with `[k·F, (k+1)·F)` for k = 0..n_founders-1 and
summing the intersection lengths.

Usage:
    python3 tools/subrepeat_scan/scan.py \\
        --fasta <fasta> \\
        --rescore <rescore.peaks.tsv> \\
        --out-prefix <prefix>

Outputs:
    <prefix>.intervals.tsv  — per-interval calls
    <prefix>.summary.tsv    — per-array occupancy summary
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------- FASTA I/O ------------------------------------------------------


def iter_fasta(path: Path):
    """Yield (record_id, sequence_bytes) per record. Whitespace-trimmed."""
    with open(path) as f:
        header = None
        buf = []
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


def read_rescore_index(
    rescore_tsv: Path, period_min: int, period_max: int
) -> tuple[dict[str, int], dict[str, list[int]]]:
    """Read a rescore peaks TSV and return:
      - record_id → founder_period (int)
      - record_id → list of kite-reported `period` values that fall in
        `[period_min, period_max]` (deduplicated, sorted ascending)
    """
    founder: dict[str, int] = {}
    periods: dict[str, set[int]] = {}
    with open(rescore_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rid = row.get("case_id") or row.get("record_id")
            if not rid:
                continue
            if rid not in founder:
                fp = row.get("founder_period", "NA")
                if fp and fp != "NA":
                    try:
                        founder[rid] = int(fp)
                    except ValueError:
                        pass
            p = row.get("period", "NA")
            if p and p != "NA":
                try:
                    pi = int(p)
                    if period_min <= pi <= period_max:
                        periods.setdefault(rid, set()).add(pi)
                except ValueError:
                    pass
    return founder, {k: sorted(v) for k, v in periods.items()}


# ---------- Core scan ------------------------------------------------------


@dataclass(frozen=True)
class Interval:
    """A candidate nested-subrepeat interval inside an array."""

    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    period: int
    score: float  # mean smoothed identity over the interval
    copies: float

    @property
    def length(self) -> int:
        return self.end - self.start


def shifted_identity(seq: np.ndarray, p: int) -> np.ndarray:
    """Per-base match indicator at shift p: 1 where seq[i] == seq[i+p].

    Returns a length-(L-p) uint8 array.
    """
    if p <= 0 or p >= len(seq):
        return np.zeros(0, dtype=np.uint8)
    return (seq[:-p] == seq[p:]).astype(np.uint8)


def windowed_match_rate(match: np.ndarray, p: int) -> np.ndarray:
    """For each i, return mean(match[i:i+p]) — the match rate inside
    one period-wide forward window starting at i.

    Returns a length-(len(match) - p + 1) float64 array.
    """
    if len(match) < p:
        return np.zeros(0, dtype=np.float64)
    csum = np.concatenate(([0], np.cumsum(match, dtype=np.int64)))
    n = len(match) - p + 1
    sums = csum[p : p + n] - csum[:n]
    return sums.astype(np.float64) / p


def find_runs(mask: np.ndarray, min_length: int) -> list[tuple[int, int]]:
    """Return list of (start, end) for contiguous True runs in `mask`
    with length ≥ min_length. End is exclusive."""
    if len(mask) == 0:
        return []
    # find rising / falling edges
    padded = np.concatenate(([False], mask, [False]))
    diff = np.diff(padded.astype(np.int8))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    runs = [(s, e) for s, e in zip(starts, ends) if (e - s) >= min_length]
    return runs


def scan_one_record(
    seq: np.ndarray,
    periods_to_test: list[int],
    id_threshold: float,
    min_copies: int,
) -> list[Interval]:
    """Scan one array for nested-subrepeat candidate intervals.

    Tests each period in `periods_to_test`. For each candidate
    period p, compute the windowed-match-rate profile (mean of
    `match` over windows of size p, where `match[i] = 1 iff
    seq[i] == seq[i+p]`). A real tandem of N copies of a motif of
    length p produces a contiguous run of length `(N − 1)·p`
    windows where the rate is ≥ id_threshold. We require at least
    `(min_copies − 1)·p` consecutive qualifying windows; the
    resulting sequence-region covers `(N · p)` bp end-to-end.
    """
    candidates: list[Interval] = []
    for p in periods_to_test:
        match = shifted_identity(seq, p)
        if len(match) < p:
            continue
        rate = windowed_match_rate(match, p)
        mask = rate >= id_threshold
        min_run_windows = (min_copies - 1) * p
        if min_run_windows < 1:
            min_run_windows = 1
        runs = find_runs(mask, min_length=min_run_windows)
        for win_start, win_end in runs:
            # win_end is exclusive over window-start indices. The
            # last window starts at win_end-1 and covers p bp, so
            # the sequence region runs to win_end-1 + p = win_end+p-1.
            seq_start = int(win_start)
            seq_end = int(win_end + p - 1)
            seq_end = min(seq_end, len(seq))
            score = float(np.mean(rate[win_start:win_end]))
            copies = (seq_end - seq_start) / p
            candidates.append(
                Interval(start=seq_start, end=seq_end, period=p,
                         score=score, copies=copies)
            )
    return candidates


def resolve_aliasing(
    candidates: list[Interval], overlap_frac: float
) -> list[Interval]:
    """Greedy: prefer shorter periods. An interval is dropped if it
    overlaps an already-kept interval by ≥ overlap_frac of its own
    length AND has a strictly larger period.
    """
    sorted_cands = sorted(candidates, key=lambda c: (c.period, -c.length))
    kept: list[Interval] = []
    for c in sorted_cands:
        drop = False
        for k in kept:
            if k.period >= c.period:
                continue
            ov_start = max(c.start, k.start)
            ov_end = min(c.end, k.end)
            ov_len = max(0, ov_end - ov_start)
            if ov_len >= overlap_frac * c.length:
                drop = True
                break
        if not drop:
            kept.append(c)
    return kept


def merge_overlapping_same_period(intervals: list[Interval]) -> list[Interval]:
    """Merge intervals at the same period that overlap or abut. Used
    after the per-period scan to coalesce broken runs."""
    if not intervals:
        return []
    by_period: dict[int, list[Interval]] = {}
    for c in intervals:
        by_period.setdefault(c.period, []).append(c)
    out: list[Interval] = []
    for p, cands in by_period.items():
        cands.sort(key=lambda c: c.start)
        cur = cands[0]
        for nxt in cands[1:]:
            if nxt.start <= cur.end:
                end = max(cur.end, nxt.end)
                cur = Interval(
                    start=cur.start, end=end, period=p,
                    score=max(cur.score, nxt.score),
                    copies=(end - cur.start) / p,
                )
            else:
                out.append(cur)
                cur = nxt
        out.append(cur)
    return out


def union_length(intervals: list[Interval]) -> int:
    """Total bp covered by the union of intervals."""
    if not intervals:
        return 0
    spans = sorted(((c.start, c.end) for c in intervals))
    total = 0
    cur_s, cur_e = spans[0]
    for s, e in spans[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return total


def per_founder_occupancy(
    intervals: list[Interval], founder_period: int, array_length: int
) -> list[tuple[int, int, int, int]]:
    """For each founder window, compute occupied bp.

    Returns list of (founder_idx, start_bp, end_bp, occupied_bp).
    """
    if founder_period <= 0:
        return []
    n_founders = max(1, array_length // founder_period)
    out = []
    for k in range(n_founders):
        win_s = k * founder_period
        win_e = min(array_length, (k + 1) * founder_period)
        # intersect every interval with this founder window
        clipped: list[Interval] = []
        for c in intervals:
            s = max(c.start, win_s)
            e = min(c.end, win_e)
            if e > s:
                clipped.append(
                    Interval(start=s, end=e, period=c.period,
                             score=c.score, copies=(e - s) / c.period)
                )
        occupied = union_length(clipped)
        out.append((k, win_s, win_e, occupied))
    return out


# ---------- Main -----------------------------------------------------------


def encode_seq(seq_bytes: bytes) -> np.ndarray:
    """Encode ACGT to small ints; non-ACGT bases become a sentinel
    (255) that never equals itself in the match function (so N-N
    pairs don't count as matches).
    """
    arr = np.frombuffer(seq_bytes.upper(), dtype=np.uint8).copy()
    table = np.full(256, 255, dtype=np.uint8)
    table[ord("A")] = 0
    table[ord("C")] = 1
    table[ord("G")] = 2
    table[ord("T")] = 3
    encoded = table[arr]
    # Replace Ns with random sentinel that breaks self-match: use 4
    # for forward, 5 for backward (we want N != N to keep them out
    # of high-identity runs).
    # Simpler: use position-dependent sentinel — guarantee mismatch.
    is_n = encoded == 255
    encoded[is_n] = (
        np.arange(len(encoded), dtype=np.uint8)[is_n] | 0x80
    )  # high bit set, position-varying
    return encoded


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, required=True,
                   help="input FASTA")
    p.add_argument("--rescore", type=Path,
                   help="rescore peaks TSV — used to read founder_period per record")
    p.add_argument("--founder-period", type=int, default=None,
                   help="apply this founder_period to every record (overrides --rescore)")
    p.add_argument("--out-prefix", type=Path, required=True,
                   help="output prefix; writes <prefix>.intervals.tsv and <prefix>.summary.tsv")
    p.add_argument("--period-min", type=int, default=20,
                   help="lower bound for periods to consider")
    p.add_argument("--period-max", type=int, default=45,
                   help="upper bound for periods to consider")
    p.add_argument("--all-periods", action="store_true",
                   help="test every integer period in [period_min, period_max] "
                        "for every record (default: restrict to kite-reported "
                        "candidates from --rescore, intersected with the bounds)")
    p.add_argument("--id-threshold", type=float, default=0.70,
                   help="min match rate inside a p-wide window for it "
                        "to qualify; tune for the expected per-copy divergence "
                        "(0.70 ~ 12 pct divergence, 0.80 ~ 7 pct)")
    p.add_argument("--min-copies", type=int, default=3,
                   help="minimum number of tandem copies — translates to "
                        "(min_copies − 1) · p consecutive qualifying windows")
    p.add_argument("--overlap-frac", type=float, default=0.5,
                   help="period-aliasing overlap threshold")
    p.add_argument("--records", default=None,
                   help="comma-separated record IDs to restrict scanning to")
    args = p.parse_args()

    if args.founder_period is None and args.rescore is None:
        sys.exit("error: pass --rescore or --founder-period")

    founder: dict[str, int] = {}
    kite_periods: dict[str, list[int]] = {}
    if args.rescore:
        founder, kite_periods = read_rescore_index(
            args.rescore, args.period_min, args.period_max
        )
        print(
            f"loaded {len(founder)} founder period(s) and "
            f"{sum(len(v) for v in kite_periods.values())} kite-reported "
            f"candidate period(s) in [{args.period_min}, {args.period_max}] "
            f"from {args.rescore}",
            file=sys.stderr,
        )

    only_records = None
    if args.records:
        only_records = set(args.records.split(","))

    intervals_out = args.out_prefix.with_suffix(".intervals.tsv")
    summary_out = args.out_prefix.with_suffix(".summary.tsv")
    intervals_out.parent.mkdir(parents=True, exist_ok=True)

    n_records_scanned = 0
    n_intervals_total = 0
    with open(intervals_out, "w") as fi, open(summary_out, "w") as fs:
        fi.write(
            "record_id\tinterval_idx\tstart_bp\tend_bp\tlength_bp\t"
            "period_bp\tcopies\tsmoothed_id_mean\n"
        )
        fs.write(
            "record_id\tarray_length\tfounder_period\tn_founders\t"
            "n_intervals\ttotal_occupied_bp\toccupancy_frac\t"
            "median_per_founder_occupancy_frac\n"
        )

        for rid, seq_bytes in iter_fasta(args.fasta):
            if only_records is not None and rid not in only_records:
                continue
            L = len(seq_bytes)
            if L < 3 * args.period_max:
                continue
            seq = encode_seq(seq_bytes)
            if args.all_periods or not args.rescore:
                periods_to_test = list(range(args.period_min,
                                             args.period_max + 1))
            else:
                periods_to_test = kite_periods.get(rid, [])
            if not periods_to_test:
                # Still write a summary row with zero intervals.
                cands: list[Interval] = []
            else:
                cands = scan_one_record(
                    seq,
                    periods_to_test=periods_to_test,
                    id_threshold=args.id_threshold,
                    min_copies=args.min_copies,
                )
            cands = merge_overlapping_same_period(cands)
            selected = resolve_aliasing(cands, overlap_frac=args.overlap_frac)
            selected.sort(key=lambda c: c.start)

            # Per-interval rows
            for i, c in enumerate(selected):
                fi.write(
                    f"{rid}\t{i}\t{c.start}\t{c.end}\t{c.length}\t"
                    f"{c.period}\t{c.copies:.2f}\t{c.score:.4f}\n"
                )

            # Summary row
            F = args.founder_period or founder.get(rid, 0)
            total_occ = union_length(selected)
            occ_frac = total_occ / L if L else 0.0
            if F > 0:
                per_f = per_founder_occupancy(selected, F, L)
                n_founders = len(per_f)
                fracs = [
                    occ / (end - start) if end > start else 0.0
                    for _, start, end, occ in per_f
                ]
                median_frac = float(np.median(fracs)) if fracs else 0.0
            else:
                n_founders = 0
                median_frac = 0.0

            fs.write(
                f"{rid}\t{L}\t{F if F else 'NA'}\t{n_founders}\t"
                f"{len(selected)}\t{total_occ}\t{occ_frac:.4f}\t"
                f"{median_frac:.4f}\n"
            )

            n_records_scanned += 1
            n_intervals_total += len(selected)

    print(
        f"scanned {n_records_scanned} record(s), "
        f"emitted {n_intervals_total} interval(s)",
        file=sys.stderr,
    )
    print(f"wrote {intervals_out}", file=sys.stderr)
    print(f"wrote {summary_out}", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the training feature matrix per record.

Inputs (per seed dir):
  - sequences.fasta
  - truth.tsv
  - kite/top3.tsv                  (from `hordetect kite-periodicity`)
  - kite/top3.tsv.peaks.tsv        (long-format)

Outputs:
  - <out>/features_seed<seed>.tsv   — one row per record, ~13 features + target

Features (final list):
  Kite (8):
    s1, s2, s3, s2_over_s1, s3_over_s1,
    family_size_best, tile_founder_ratio, tile_jitter
  Period structure (2):
    d1, log_d1_over_L
  Sequence diversity (3):
    distinct_kmers_per_bp, kmer_entropy, singletons_ratio

Target:
  hor_signal (continuous, from truth.tsv). NA replaced by 0.0 (null arrays).
  hor_binary (= truth.hor_order > 1).
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional


KITE_K = 6
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_fasta(path: Path):
    cur, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    yield cur, "".join(chunks)
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if cur is not None:
            yield cur, "".join(chunks)


def diversity_features(seq: str, k: int = KITE_K) -> dict:
    """Distinct k-mer count, Shannon entropy, singletons ratio. Skip
    k-mers containing N (matches kite.R / our Rust kite path)."""
    L = len(seq)
    if L < k + 1:
        return {"distinct_kmers_per_bp": float("nan"),
                "kmer_entropy": float("nan"),
                "singletons_ratio": float("nan")}
    counts: dict[str, int] = defaultdict(int)
    total = 0
    for i in range(L - k + 1):
        s = seq[i:i + k]
        if "N" in s:
            continue
        counts[s] += 1
        total += 1
    if total == 0:
        return {"distinct_kmers_per_bp": float("nan"),
                "kmer_entropy": float("nan"),
                "singletons_ratio": float("nan")}
    distinct = len(counts)
    singletons = sum(1 for v in counts.values() if v == 1)
    # Shannon entropy of k-mer frequency distribution, in bits.
    entropy = 0.0
    for c in counts.values():
        p = c / total
        entropy -= p * math.log2(p)
    return {
        "distinct_kmers_per_bp": distinct / L,
        "kmer_entropy": entropy,
        "singletons_ratio": singletons / total,
    }


def find_best_family(peaks: list[dict], qmax: int = 30,
                     tol_bp: int = 5, tol_rel: float = 0.02,
                     lo_period: int = 15, min_size: int = 2,
                     min_founder_top1_share: float = 0.5,
                     require_top_k: int = 3,
                     ) -> dict:
    """Mirror hor_call.rs::find_best_family + family_metrics for the
    feature extraction. Returns family_size_best, tile_founder_ratio,
    tile_period, tile_jitter (peaks within ±15 % of tile)."""
    if not peaks:
        return {"family_size_best": 0, "tile_founder_ratio": 0.0,
                "tile_period": float("nan"), "tile_jitter": 0}
    peaks = sorted(peaks, key=lambda p: -p["score"])
    top1_score = peaks[0]["score"]
    top_k_periods = [p["period"] for p in peaks[:require_top_k]]
    def matches(p, m):
        if m == 0:
            return None
        k = round(p / m)
        if k < 1 or k > qmax:
            return None
        expected = k * m
        diff = abs(p - expected)
        tol = max(tol_bp, int(tol_rel * expected))
        return k if diff <= tol else None
    for cand in peaks:
        m_f = cand["period"]
        if m_f < lo_period:
            continue
        if top1_score > 0 and cand["score"] < min_founder_top1_share * top1_score:
            continue
        family = []
        best_tile_k, best_tile_period, best_tile_score = 1, m_f, -1.0
        for p in peaks:
            k = matches(p["period"], m_f)
            if k is None:
                continue
            family.append((p["period"], p["score"], k))
            if k >= 2 and p["score"] > best_tile_score:
                best_tile_score = p["score"]
                best_tile_k = k
                best_tile_period = p["period"]
        if len(family) < min_size:
            continue
        if not all(matches(tp, m_f) is not None for tp in top_k_periods):
            continue
        tile_jitter = sum(
            1 for p in peaks
            if abs(p["period"] - best_tile_period) <= 0.15 * best_tile_period
        )
        ratio = (best_tile_score / cand["score"]
                 if cand["score"] > 0 and best_tile_score > 0 else 0.0)
        return {
            "family_size_best": len(family),
            "tile_founder_ratio": ratio,
            "tile_period": best_tile_period,
            "tile_jitter": tile_jitter,
        }
    return {"family_size_best": 0, "tile_founder_ratio": 0.0,
            "tile_period": float("nan"), "tile_jitter": 0}


def load_peaks_by_case(path: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            r = dict(zip(header, cols))
            out[r["case_id"]].append({
                "rank": int(r["rank"]),
                "period": int(r["period"]),
                "peak_height": float(r["peak_height"]),
                "score": float(r["score"]),
                "score2": float(r["score2"]),
                "background": float(r["background"]),
            })
    return out


def load_truth(path: Path) -> dict[str, dict]:
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            r = dict(zip(header, cols))
            out[r["case_id"]] = r
    return out


def parse_float_or_zero(s: str) -> float:
    try:
        v = float(s)
        if v != v:  # NaN
            return 0.0
        return v
    except Exception:
        return 0.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed-dir", type=Path, required=True,
                    help="e.g. eval/training_data/sim_seed101")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--kite-k", type=int, default=KITE_K)
    args = ap.parse_args(argv)

    fasta = args.seed_dir / "sequences.fasta"
    truth_path = args.seed_dir / "truth.tsv"
    peaks_path = args.seed_dir / "kite" / "top3.tsv.peaks.tsv"
    for p in [fasta, truth_path, peaks_path]:
        if not p.exists():
            print(f"error: missing {p}", file=sys.stderr)
            return 1
    truth = load_truth(truth_path)
    peaks = load_peaks_by_case(peaks_path)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "case_id", "stratum", "array_length",
        "s1", "s2", "s3", "s2_over_s1", "s3_over_s1",
        "family_size_best", "tile_founder_ratio", "tile_jitter",
        "d1", "d2", "d3", "log_d1_over_L",
        "d2_over_d1", "d3_over_d1",
        "max_d_top3_over_min_d_top3",
        # Extended kite peaks (positions 4..10 by score); each is the
        # period of the k-th kite peak after score2_norm + peak>bg
        # filtering. Zero if fewer than k peaks exist for this record.
        "d4", "d5", "d6", "d7", "d8", "d9", "d10",
        "family_founder_d", "family_tile_d",   # for downstream homology probing
        "distinct_kmers_per_bp", "kmer_entropy", "singletons_ratio",
        "hor_signal", "hor_binary",
    ]
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for case_id, seq in parse_fasta(fasta):
            L = len(seq)
            t = truth.get(case_id, {})
            stratum = case_id.rsplit("_", 1)[0]
            hor_order = int(t.get("hor_order", "1"))
            hor_binary = 1 if hor_order > 1 else 0
            hor_signal_raw = t.get("hor_signal", "NA")
            hor_signal = parse_float_or_zero(hor_signal_raw)
            div = diversity_features(seq, args.kite_k)
            pks = peaks.get(case_id, [])
            pks_sorted = sorted(pks, key=lambda p: p["rank"])
            s = [p["score"] for p in pks_sorted]
            d = [p["period"] for p in pks_sorted]
            def pad(arr, n, fill=0.0):
                return list(arr) + [fill] * (n - len(arr))
            s1, s2, s3 = pad(s, 3, 0.0)[:3]
            d1, d2, d3 = pad(d, 3, 0)[:3]
            d4, d5, d6, d7, d8, d9, d10 = pad(d, 10, 0)[3:10]
            s2os1 = (s2 / s1) if s1 > 0 else 0.0
            s3os1 = (s3 / s1) if s1 > 0 else 0.0
            log_d1 = (math.log(d1 / L) if d1 > 0 and L > 0 else 0.0)
            d2od1 = (d2 / d1) if d1 > 0 and d2 > 0 else 0.0
            d3od1 = (d3 / d1) if d1 > 0 and d3 > 0 else 0.0
            d_nonzero = [x for x in (d1, d2, d3) if x > 0]
            max_d_over_min_d = (max(d_nonzero) / min(d_nonzero)
                                if len(d_nonzero) >= 2 else 1.0)
            fam = find_best_family(pks)
            # Founder period: when a family was found, fam already has the
            # founder (the candidate m_f). Reconstruct from the family's
            # smallest k=1 member (== founder period of the best m_f).
            # As a robust shortcut: when family was found via the family
            # search, the founder equals the m_f the search chose. We
            # recompute it here by re-running the search and taking the
            # first qualifying candidate.
            family_founder_d = 0
            family_tile_d = 0
            if fam["family_size_best"] > 0 and not (fam["tile_period"] != fam["tile_period"]):
                # find_best_family wrote tile_period; founder is whichever
                # candidate m_f produced this family. We re-scan to extract
                # it (cheap given top peaks <= 30).
                peaks_sorted = sorted(pks, key=lambda p: -p["score"])
                top1 = peaks_sorted[0]["score"] if peaks_sorted else 0.0
                tol_bp, tol_rel, qmax = 5, 0.02, 30
                require_top_k = 3
                top_k_periods = [p["period"] for p in peaks_sorted[:require_top_k]]
                def fits(p, m):
                    if m == 0:
                        return False
                    k = round(p / m)
                    if k < 1 or k > qmax:
                        return False
                    e = k * m
                    tol = max(tol_bp, int(tol_rel * e))
                    return abs(p - e) <= tol
                for cand in peaks_sorted:
                    m_f = cand["period"]
                    if m_f < 15:
                        continue
                    if top1 > 0 and cand["score"] < 0.5 * top1:
                        continue
                    fam_count = sum(1 for p in peaks_sorted if fits(p["period"], m_f))
                    if fam_count < 2:
                        continue
                    if not all(fits(tp, m_f) for tp in top_k_periods):
                        continue
                    family_founder_d = m_f
                    family_tile_d = int(fam["tile_period"])
                    break
            row = [
                case_id, stratum, str(L),
                f"{s1:.6f}", f"{s2:.6f}", f"{s3:.6f}",
                f"{s2os1:.6f}", f"{s3os1:.6f}",
                str(fam["family_size_best"]),
                f"{fam['tile_founder_ratio']:.6f}",
                str(fam["tile_jitter"]),
                str(d1), str(d2), str(d3),
                f"{log_d1:.6f}",
                f"{d2od1:.6f}", f"{d3od1:.6f}",
                f"{max_d_over_min_d:.6f}",
                str(d4), str(d5), str(d6), str(d7), str(d8), str(d9), str(d10),
                str(family_founder_d), str(family_tile_d),
                f"{div['distinct_kmers_per_bp']:.6f}",
                f"{div['kmer_entropy']:.6f}",
                f"{div['singletons_ratio']:.6f}",
                f"{hor_signal:.6f}",
                str(hor_binary),
            ]
            out.write("\t".join(row) + "\n")
    print(f"wrote features for {len(truth)} records → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

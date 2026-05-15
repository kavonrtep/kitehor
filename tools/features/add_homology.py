#!/usr/bin/env python3
"""Add `h_d1` and `h_founder` columns to a features TSV by probing
block-mean homology via `hordetect probe-periods` (the Path B bridge).

For each record:
  - probe at d1 (kite top-1 period; always probed)
  - probe at family_founder_d when the kite family has a real founder
    (≠ d1, ≠ 0); else h_founder := h_d1

These two features should help the model distinguish:
  - real HOR (founder homology ~0.5-0.7, tile homology ~0.85+)
  - pure tandem (single-period homology ~0.85+ at the monomer)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--features", type=Path, required=True,
                    help="Existing features.tsv produced by extract_features.py")
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--hordetect", type=Path,
                    default=Path("hordetect/target/release/hordetect"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not args.hordetect.exists():
        debug = Path("hordetect/target/debug/hordetect")
        if debug.exists():
            args.hordetect = debug
        else:
            print(f"error: hordetect binary not found", file=sys.stderr)
            return 1

    # Read features.
    with open(args.features) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh]
    idx = {c: i for i, c in enumerate(header)}
    for needed in ("case_id", "d1", "family_founder_d"):
        if needed not in idx:
            print(f"error: column {needed} missing", file=sys.stderr)
            return 1

    # Build the probe TSV. Use a set per (case_id, period) to avoid
    # duplicate work.
    probes: set[tuple[str, int]] = set()
    for r in rows:
        case = r["case_id"]
        try:
            d1 = int(r["d1"])
            if d1 > 0:
                probes.add((case, d1))
            fd = int(r["family_founder_d"])
            if fd > 0:
                probes.add((case, fd))
        except Exception:
            pass
    print(f"unique (case, period) probes: {len(probes)}", file=sys.stderr)

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        probe_in = td_path / "probes.tsv"
        probe_out = td_path / "probes.out.tsv"
        with open(probe_in, "w") as fh:
            fh.write("case_id\tperiod\n")
            for case, p in sorted(probes):
                fh.write(f"{case}\t{p}\n")
        print("running hordetect probe-periods ...", file=sys.stderr)
        subprocess.run([
            str(args.hordetect), "probe-periods",
            str(args.fasta),
            "--periods", str(probe_in),
            "-o", str(probe_out),
        ], check=True)

        # Read homology results.
        h: dict[tuple[str, int], float] = {}
        with open(probe_out) as fh:
            next(fh)  # header
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                case, period, hh = cols[0], int(cols[1]), cols[2]
                if hh == "NA":
                    continue
                h[(case, period)] = float(hh)

    # Write extended TSV.
    new_cols = header + ["h_d1", "h_founder"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(new_cols) + "\n")
        for r in rows:
            case = r["case_id"]
            try:
                d1 = int(r["d1"])
                fd = int(r["family_founder_d"])
            except Exception:
                d1, fd = 0, 0
            h_d1 = h.get((case, d1), float("nan")) if d1 > 0 else float("nan")
            h_f  = h.get((case, fd), h_d1) if fd > 0 else h_d1
            def fmt(v):
                return "NA" if v != v else f"{v:.4f}"
            row = list(r.values()) + [fmt(h_d1), fmt(h_f)]
            fh.write("\t".join(row) + "\n")
    print(f"wrote {len(rows)} rows → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

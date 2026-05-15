#!/usr/bin/env python3
"""Generate synthetic tandem-repeat arrays for testing KITE / HOR scoring.

For each row of an input parameter TSV, builds one tandem-repeat array with a
controlled higher-order-repeat (HOR) structure and writes:

  sequences.fasta  one record per case, header ``>{case_id}``
  truth.tsv        construction parameters + measured signal metrics per case
  monomers.tsv     per-monomer lattice position and coordinates in the array
  events.tsv       log of conversion events applied to each array

Per-case simulation pipeline
----------------------------
1. Generate a random monomer M of length L.
2. Derive N founder monomers from M using intra-block rates -> HOR block B
   (for N=1 the founder equals M; intra rates are ignored).
3. Replicate B into K block copies, mutating every monomer copy independently
   with inter-block rates.
4. Apply post-construction conversion events to model concerted evolution:
   - block-level: copy one HOR block onto another (strengthens HOR signal)
   - monomer-level: copy one monomer onto another, regardless of block
     (degrades HOR signal; for N=1 produces a homogenized satellite)

Parameter TSV columns
---------------------
case_id              free-form identifier, used as FASTA record name
monomer_len          base monomer length L (int)
hor_order            HOR period N (int, >=1)
n_blocks             number of HOR block copies K (int, >=1)
sub_rate_intra       per-base substitution rate when deriving founders
sub_rate_inter       per-base substitution rate per block copy
indel_rate_intra     per-position indel rate (50/50 ins/del) for founders
indel_rate_inter     per-position indel rate per block copy
block_conversions    number of block-level conversion events (int, >=0)
monomer_conversions  number of monomer-level conversion events (int, >=0)
seed                 optional per-case integer seed; if blank, derived from
                     --seed and case_id

Identity in truth.tsv is computed as ``1 - edit_distance / max(len_a, len_b)``
(Levenshtein) over a random subsample of monomer pairs in each category. The
subsample keeps runtime bounded; the metric itself is alignment-aware and not
biased by indels.
"""

import argparse
import csv
import hashlib
import random
import sys
from pathlib import Path

BASES = "ACGT"
ALT = {b: tuple(c for c in BASES if c != b) for b in BASES}

EXAMPLE_PARAMS = """\
case_id\tmonomer_len\thor_order\tn_blocks\tsub_rate_intra\tsub_rate_inter\tindel_rate_intra\tindel_rate_inter\tblock_conversions\tmonomer_conversions\tseed
no_hor_clean\t170\t1\t40\t0\t0.02\t0\t0.005\t0\t0\t
no_hor_homogenized\t170\t1\t40\t0\t0.05\t0\t0.01\t0\t30\t
hor4_strong\t170\t4\t12\t0.10\t0.02\t0.01\t0.005\t8\t0\t
hor4_decaying\t170\t4\t12\t0.10\t0.05\t0.01\t0.01\t4\t12\t
hor3_no_conversion\t170\t3\t15\t0.08\t0.03\t0.005\t0.005\t0\t0\t
"""


def random_dna(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def mutate(seq, sub_rate, indel_rate, rng):
    """Apply substitutions and indels to a DNA string.

    Each position independently: with probability ``indel_rate`` an indel
    occurs (50/50 insertion of a random base before the position, or
    deletion); otherwise, with probability ``sub_rate`` the base is replaced
    by one of the other three bases uniformly.
    """
    out = []
    for base in seq:
        if rng.random() < indel_rate:
            if rng.random() < 0.5:
                out.append(rng.choice(BASES))
                out.append(base)
        else:
            if rng.random() < sub_rate:
                out.append(rng.choice(ALT[base]))
            else:
                out.append(base)
    return "".join(out)


def simulate_array(monomer_len, hor_order, n_blocks,
                   sub_intra, sub_inter, indel_intra, indel_inter,
                   block_conv, monomer_conv, rng,
                   submono_k=1):
    """Return (monomers, events) for one simulated array.

    monomers: list of dicts ``{block_idx, founder_idx, seq}`` in array order.
              ``block_idx`` and ``founder_idx`` are the *original* lattice
              coordinates; sequence content may have been overwritten by a
              later conversion event but the lattice labels are preserved.
    events:   list of dicts ``{event_order, scope, source_idx, target_idx}``.
              For ``scope == 'block'``, indices are block indices; for
              ``scope == 'monomer'``, indices are flat monomer indices into
              ``monomers``.

    submono_k > 1 builds the base monomer by tiling a smaller random motif
    of length ``monomer_len // submono_k`` ``submono_k`` times. This injects
    an internal sub-periodicity inside each monomer; useful for generating
    arrays where the autocorrelation spectrum has a sub-monomer peak (e.g.
    at L/4) without there being any biological HOR. Sub-monomer arrays
    must NOT be called HOR by the detector.
    """
    if submono_k >= 2:
        sub_len = max(1, monomer_len // submono_k)
        sub_motif = random_dna(sub_len, rng)
        base = (sub_motif * submono_k)[:monomer_len]
    else:
        base = random_dna(monomer_len, rng)

    if hor_order == 1:
        founders = [base]
    else:
        founders = [mutate(base, sub_intra, indel_intra, rng)
                    for _ in range(hor_order)]

    monomers = []
    for b in range(n_blocks):
        for f in range(hor_order):
            seq = mutate(founders[f], sub_inter, indel_inter, rng)
            monomers.append({"block_idx": b, "founder_idx": f, "seq": seq})

    events = []

    for _ in range(block_conv):
        if n_blocks < 2:
            break
        s, t = rng.sample(range(n_blocks), 2)
        for f in range(hor_order):
            monomers[t * hor_order + f]["seq"] = monomers[s * hor_order + f]["seq"]
        events.append({"event_order": len(events), "scope": "block",
                       "source_idx": s, "target_idx": t})

    total = len(monomers)
    for _ in range(monomer_conv):
        if total < 2:
            break
        s, t = rng.sample(range(total), 2)
        monomers[t]["seq"] = monomers[s]["seq"]
        events.append({"event_order": len(events), "scope": "monomer",
                       "source_idx": s, "target_idx": t})

    return monomers, events


def edit_distance(a, b):
    """Levenshtein distance between two strings (Wagner-Fischer, O(len_a*len_b))."""
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr.append(min(curr[j - 1] + 1, prev[j] + 1, prev[j - 1] + cost))
        prev = curr
    return prev[-1]


def aligned_identity(a, b):
    """Identity = 1 - edit_distance / max(len_a, len_b). NaN when both empty."""
    m = max(len(a), len(b))
    if m == 0:
        return float("nan")
    return 1.0 - edit_distance(a, b) / m


def diagnostic_metrics(monomers, hor_order, rng, max_pairs_per_category=40):
    """Mean alignment-based identity in three pair categories.

    Pairs are categorised by lattice labels (block_idx, founder_idx), then
    each category is independently sub-sampled to ``max_pairs_per_category``
    pairs to bound runtime.
    """
    intra_block, homologous, cross_position = [], [], []
    n = len(monomers)
    for i in range(n):
        for j in range(i + 1, n):
            mi, mj = monomers[i], monomers[j]
            if mi["block_idx"] == mj["block_idx"]:
                if hor_order > 1:
                    intra_block.append((i, j))
            elif mi["founder_idx"] == mj["founder_idx"]:
                homologous.append((i, j))
            else:
                cross_position.append((i, j))

    def mean_id(pairs):
        if not pairs:
            return float("nan")
        if len(pairs) > max_pairs_per_category:
            pairs = rng.sample(pairs, max_pairs_per_category)
        vals = [aligned_identity(monomers[i]["seq"], monomers[j]["seq"])
                for i, j in pairs]
        vals = [v for v in vals if v == v]
        return sum(vals) / len(vals) if vals else float("nan")

    intra = mean_id(intra_block)
    homol = mean_id(homologous)
    cross = mean_id(cross_position)
    signal = homol - cross if (homol == homol and cross == cross) else float("nan")

    return {
        "mean_intra_block_id": intra,
        "mean_homologous_id": homol,
        "mean_cross_position_id": cross,
        "hor_signal": signal,
    }


def fmt(v):
    if isinstance(v, float):
        return "NA" if v != v else f"{v:.4f}"
    return str(v)


def derive_seed(master_seed, case_id):
    h = hashlib.sha256(f"{master_seed}:{case_id}".encode()).digest()
    return int.from_bytes(h[:8], "little")


PARAM_FIELDS = (
    "case_id", "monomer_len", "hor_order", "n_blocks",
    "sub_rate_intra", "sub_rate_inter",
    "indel_rate_intra", "indel_rate_inter",
    "block_conversions", "monomer_conversions",
    "submono_k", "seed",
)
TRUTH_FIELDS = PARAM_FIELDS + (
    "array_length", "n_monomers",
    "mean_intra_block_id", "mean_homologous_id",
    "mean_cross_position_id", "hor_signal",
)


def parse_case(row, master_seed):
    case_id = row["case_id"].strip()
    seed_raw = (row.get("seed") or "").strip()
    seed = int(seed_raw) if seed_raw else derive_seed(master_seed, case_id)
    submono_raw = (row.get("submono_k") or "").strip()
    return {
        "case_id": case_id,
        "monomer_len": int(row["monomer_len"]),
        "hor_order": int(row["hor_order"]),
        "n_blocks": int(row["n_blocks"]),
        "sub_rate_intra": float(row["sub_rate_intra"]),
        "sub_rate_inter": float(row["sub_rate_inter"]),
        "indel_rate_intra": float(row["indel_rate_intra"]),
        "indel_rate_inter": float(row["indel_rate_inter"]),
        "block_conversions": int(row.get("block_conversions") or 0),
        "monomer_conversions": int(row.get("monomer_conversions") or 0),
        "submono_k": int(submono_raw) if submono_raw else 1,
        "seed": seed,
    }


def write_fasta_record(fh, name, seq, wrap=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), wrap):
        fh.write(seq[i:i + wrap] + "\n")


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("-p", "--params", help="parameter grid TSV (see --example-params)")
    p.add_argument("-o", "--outdir", help="output directory")
    p.add_argument("-s", "--seed", type=int, default=42,
                   help="master seed used when a row has no seed (default: 42)")
    p.add_argument("--example-params", action="store_true",
                   help="print an example parameter TSV to stdout and exit")
    args = p.parse_args()

    if args.example_params:
        sys.stdout.write(EXAMPLE_PARAMS)
        return

    if not args.params or not args.outdir:
        p.error("--params and --outdir are required (or use --example-params)")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_path = outdir / "sequences.fasta"
    truth_path = outdir / "truth.tsv"
    monomers_path = outdir / "monomers.tsv"
    events_path = outdir / "events.tsv"

    with open(args.params) as fh, \
         open(fasta_path, "w") as fa, \
         open(truth_path, "w") as tr, \
         open(monomers_path, "w") as mo, \
         open(events_path, "w") as ev:

        tr.write("\t".join(TRUTH_FIELDS) + "\n")
        mo.write("case_id\tmonomer_idx\tblock_idx\tfounder_idx\tstart\tend\tlength\n")
        ev.write("case_id\tevent_order\tscope\tsource_idx\ttarget_idx\n")

        reader = csv.DictReader(fh, delimiter="\t")
        for raw in reader:
            case = parse_case(raw, args.seed)
            rng = random.Random(case["seed"])

            for name in ("indel_rate_intra", "indel_rate_inter"):
                if case[name] > 0.1:
                    sys.stderr.write(
                        f"[warn] {case['case_id']}: {name}={case[name]} is "
                        f"high; monomer lengths may degenerate\n"
                    )

            monomers, events = simulate_array(
                case["monomer_len"], case["hor_order"], case["n_blocks"],
                case["sub_rate_intra"], case["sub_rate_inter"],
                case["indel_rate_intra"], case["indel_rate_inter"],
                case["block_conversions"], case["monomer_conversions"],
                rng,
                submono_k=case["submono_k"],
            )
            metrics = diagnostic_metrics(monomers, case["hor_order"], rng)

            seq = "".join(m["seq"] for m in monomers)
            write_fasta_record(fa, case["case_id"], seq)

            truth_row = {**case, "array_length": len(seq),
                         "n_monomers": len(monomers), **metrics}
            tr.write("\t".join(fmt(truth_row[c]) for c in TRUTH_FIELDS) + "\n")

            offset = 0
            for idx, m in enumerate(monomers):
                length = len(m["seq"])
                mo.write(f"{case['case_id']}\t{idx}\t{m['block_idx']}\t"
                         f"{m['founder_idx']}\t{offset}\t{offset + length}\t{length}\n")
                offset += length

            for e in events:
                ev.write(f"{case['case_id']}\t{e['event_order']}\t{e['scope']}\t"
                         f"{e['source_idx']}\t{e['target_idx']}\n")

            sig = metrics["hor_signal"]
            sig_str = "NA" if sig != sig else f"{sig:.3f}"
            sys.stderr.write(
                f"[ok] {case['case_id']}: {len(monomers)} monomers, "
                f"{len(seq)} bp, hor_signal={sig_str}\n"
            )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate ~1600 v2-simulator YAML configs across 9 categories.

Output tree (relative to this script's parent):

    configs/01_simple_tr/         200 configs
    configs/02_hor_clean/         600
    configs/03_hor_wobble/        100
    configs/04_hor_shift/         200
    configs/05_hor_insertion/     100
    configs/06_hor_event/         200
    configs/07_mixed/             100
    configs/08_random/             50
    configs/09_gc_bias/            50
                                ----
                                 1600

Each config has a deterministic `seed` derived from its filename so the
corpus is byte-reproducible across regenerations. A `manifest.tsv` at
the corpus root lists every (case_id, category, key parameters).
"""

import itertools
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parent
CFG_ROOT = ROOT / "configs"

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def stable_seed(name: str) -> int:
    """Deterministic small-int seed from filename so the YAML is
    self-contained and reproducible."""
    h = hashlib.md5(name.encode()).digest()
    return int.from_bytes(h[:4], "little")  # in [0, 2^32)


def write_yaml(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seed = stable_seed(path.name)
    header = f"schema_version: 1\nseed: {seed}\n"
    path.write_text(header + body.lstrip("\n"))


MANIFEST_ROWS: list[dict] = []


def record(case_id: str, category: str, **kwargs) -> None:
    row = {"case_id": case_id, "category": category}
    row.update(kwargs)
    MANIFEST_ROWS.append(row)


# ----------------------------------------------------------------------
# 01: simple_TR  (200)
# ----------------------------------------------------------------------

def category_01_simple_tr() -> int:
    n = 0
    monomers = [100, 170, 250, 500, 1000]   # 5
    n_copies  = [200, 500, 1000, 2000]       # 4
    mut       = [0.0, 0.01, 0.02, 0.05, 0.10]  # 5
    indel     = [0.0, 0.01]                  # 2
    for L, nc, m, ir in itertools.product(monomers, n_copies, mut, indel):
        # 5 × 4 × 5 × 2 = 200
        name = f"st_L{L:04d}_n{nc:04d}_m{int(m*100):02d}_i{int(ir*100):02d}.yaml"
        body = f"""
global:
  mutation_rate: {m}
  indel_rate: {ir}
templates:
  m:
    type: monomer
    monomer_length_bp: {L}
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: {nc}
"""
        write_yaml(CFG_ROOT / "01_simple_tr" / name, body)
        record(name.replace(".yaml", ""), "simple_tr",
               monomer_length_bp=L, k=1, n_copies=nc, divergence=0.0,
               mutation_rate=m, indel_rate=ir)
        n += 1
    return n


# ----------------------------------------------------------------------
# 02: hor_clean  (600)
# ----------------------------------------------------------------------

def category_02_hor_clean() -> int:
    n = 0
    monomers   = [100, 171, 250]              # 3
    ks         = [3, 6, 8, 12, 16]            # 5
    n_copies   = [50, 100, 200, 500]          # 4
    divergence = [0.05, 0.10, 0.15, 0.25, 0.40]  # 5
    mut        = [0.01, 0.02]                 # 2
    for L, k, nc, d, m in itertools.product(monomers, ks, n_copies, divergence, mut):
        # 3 × 5 × 4 × 5 × 2 = 600
        name = (f"hc_L{L:04d}_k{k:02d}_n{nc:04d}"
                f"_d{int(d*100):02d}_m{int(m*100):02d}.yaml")
        body = f"""
global:
  mutation_rate: {m}
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: {d}
structure:
  - type: HOR
    template: t
    n_copies: {nc}
"""
        write_yaml(CFG_ROOT / "02_hor_clean" / name, body)
        record(name.replace(".yaml", ""), "hor_clean",
               monomer_length_bp=L, k=k, n_copies=nc, divergence=d,
               mutation_rate=m, indel_rate=0.0)
        n += 1
    return n


# ----------------------------------------------------------------------
# 03: hor_wobble  (100)
# ----------------------------------------------------------------------

def category_03_hor_wobble() -> int:
    n = 0
    monomers = [100, 171]               # 2
    ks       = [6, 12]                  # 2
    amps     = [0.5, 1.0, 2.0, 3.0, 5.0]  # 5
    periods  = [0, 200, 500, 1000, 2000]  # 5  (0 = aperiodic random_walk)
    for L, k, A, P in itertools.product(monomers, ks, amps, periods):
        model = "random_walk" if P == 0 else "sinusoidal"
        nc = 200
        name = (f"wb_L{L:04d}_k{k:02d}_A{int(A*10):02d}"
                f"_P{P:04d}.yaml")
        body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {nc}
modifiers:
  - wobble:
      amplitude_bp: {A}
      period_rows: {P}
      model: {model}
"""
        write_yaml(CFG_ROOT / "03_hor_wobble" / name, body)
        record(name.replace(".yaml", ""), "hor_wobble",
               monomer_length_bp=L, k=k, n_copies=nc, divergence=0.15,
               mutation_rate=0.02, indel_rate=0.0,
               wobble_amplitude=A, wobble_period_rows=P)
        n += 1
    return n


# ----------------------------------------------------------------------
# 04: hor_shift  (200)
# ----------------------------------------------------------------------

def category_04_hor_shift() -> int:
    n = 0
    monomers = [171, 250]                                        # 2
    ks       = [6, 12]                                           # 2
    blocks   = [25, 50, 100, 200, 500]                           # 5
    offsets  = [-50, -25, -10, 10, 25, 50, 100, 150, 200, 300]   # 10
    for L, k, nc, off in itertools.product(monomers, ks, blocks, offsets):
        # 2 × 2 × 5 × 10 = 200
        # Reject configs where negative offset violates the
        # monomer_length/2 rule (Q5).
        if off < 0 and abs(off) > L // 2:
            continue
        name = (f"sh_L{L:04d}_k{k:02d}_n{nc:04d}"
                f"_off{off:+05d}.yaml")
        body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {nc}
  - type: SHIFT
    offset_bp: {off}
  - type: HOR
    template: t
    n_copies: {nc}
"""
        write_yaml(CFG_ROOT / "04_hor_shift" / name, body)
        record(name.replace(".yaml", ""), "hor_shift",
               monomer_length_bp=L, k=k, n_copies=2 * nc, divergence=0.15,
               mutation_rate=0.02, indel_rate=0.0, shift_offset_bp=off)
        n += 1
    return n


# ----------------------------------------------------------------------
# 05: hor_insertion  (100)
# ----------------------------------------------------------------------

def category_05_hor_insertion() -> int:
    n = 0
    kinds    = ["random", "AT_rich", "GC_rich", "retro_like", "segdup_like"]  # 5
    lengths  = [500, 2000, 5000, 10000]                                       # 4
    pre_n    = [25, 50, 100, 200, 500]                                        # 5
    for kind, ln, nc in itertools.product(kinds, lengths, pre_n):
        # 5 × 4 × 5 = 100
        L = 171
        k = 12
        name = (f"in_{kind}_L{ln:05d}_n{nc:04d}.yaml")
        body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {nc}
  - type: INSERTION
    length_bp: {ln}
    kind: {kind}
  - type: HOR
    template: t
    n_copies: {nc}
"""
        write_yaml(CFG_ROOT / "05_hor_insertion" / name, body)
        record(name.replace(".yaml", ""), "hor_insertion",
               monomer_length_bp=L, k=k, n_copies=2 * nc, divergence=0.15,
               mutation_rate=0.02, indel_rate=0.0,
               insertion_kind=kind, insertion_length_bp=ln)
        n += 1
    return n


# ----------------------------------------------------------------------
# 06: hor_event  (200) — HYBRID / INVERSION / DUPLICATION / DELETION × 50 each
# ----------------------------------------------------------------------

def category_06_hor_event() -> int:
    n = 0
    L = 171
    k = 12
    base_n = 200

    # HYBRID: 50 = 5 at_copy × 5 slot × 2 source-pair
    at_copies   = [10, 50, 100, 150, 199]
    target_slots = [1, 3, 6, 8, 11]                            # 5 target slots
    slots_choice = [(1, 2), (3, 4), (6, 7), (7, 6), (1, 12)]   # 5 pairs (use first 2)
    for ac in at_copies:
        for s_target in target_slots:
            for pair in slots_choice[:2]:  # 2 source pairs
                if s_target > k:
                    continue
                # use pair as source_slots; s_target as the replaced slot
                ss1, ss2 = pair
                if ss1 > k or ss2 > k:
                    continue
                name = f"ev_HY_ac{ac:03d}_s{s_target:02d}_src{ss1}-{ss2}.yaml"
                body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: 0.20
structure:
  - type: HOR
    template: t
    n_copies: {base_n}
post_generation:
  - type: HYBRID
    block: 0
    at_copy: {ac}
    slot: {s_target}
    source_slots: [{ss1}, {ss2}]
"""
                write_yaml(CFG_ROOT / "06_hor_event" / name, body)
                record(name.replace(".yaml", ""), "hor_event_hybrid",
                       monomer_length_bp=L, k=k, n_copies=base_n,
                       divergence=0.20, mutation_rate=0.02, indel_rate=0.0,
                       event_type="HYBRID", at_copy=ac, slot=s_target,
                       source_slots=f"{ss1},{ss2}")
                n += 1

    # INVERSION: 50 = 5 start × 5 length × 2 (k=12 sizes)
    starts  = [10, 50, 100, 150, 180]
    lengths = [1, 2, 5, 10, 20]
    for st in starts:
        for ln in lengths:
            for variant in (0, 1):
                base = base_n if variant == 0 else 300
                if st + ln - 1 > base:
                    continue
                name = f"ev_IV_st{st:03d}_ln{ln:02d}_n{base:03d}.yaml"
                body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {base}
post_generation:
  - type: INVERSION
    block: 0
    start_copy: {st}
    length_copies: {ln}
"""
                write_yaml(CFG_ROOT / "06_hor_event" / name, body)
                record(name.replace(".yaml", ""), "hor_event_inversion",
                       monomer_length_bp=L, k=k, n_copies=base,
                       divergence=0.15, mutation_rate=0.02, indel_rate=0.0,
                       event_type="INVERSION", start_copy=st, length_copies=ln)
                n += 1

    # DUPLICATION: 50 = 5 starts × 5 lengths × 2 (k variants)
    for st in starts:
        for ln in lengths:
            for k_var in (12, 8):
                if st + ln - 1 > base_n:
                    continue
                name = f"ev_DP_st{st:03d}_ln{ln:02d}_k{k_var:02d}.yaml"
                body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k_var}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {base_n}
post_generation:
  - type: DUPLICATION
    block: 0
    start_copy: {st}
    length_copies: {ln}
"""
                write_yaml(CFG_ROOT / "06_hor_event" / name, body)
                record(name.replace(".yaml", ""), "hor_event_duplication",
                       monomer_length_bp=L, k=k_var, n_copies=base_n,
                       divergence=0.15, mutation_rate=0.02, indel_rate=0.0,
                       event_type="DUPLICATION", start_copy=st, length_copies=ln)
                n += 1

    # DELETION: 50 = 5 starts × 5 lengths × 2 (k variants)
    for st in starts:
        for ln in lengths:
            for k_var in (12, 8):
                if st + ln - 1 > base_n:
                    continue
                name = f"ev_DL_st{st:03d}_ln{ln:02d}_k{k_var:02d}.yaml"
                body = f"""
global:
  mutation_rate: 0.02
templates:
  t:
    type: HOR_slots
    monomer_length_bp: {L}
    k: {k_var}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: t
    n_copies: {base_n}
post_generation:
  - type: DELETION
    block: 0
    start_copy: {st}
    length_copies: {ln}
"""
                write_yaml(CFG_ROOT / "06_hor_event" / name, body)
                record(name.replace(".yaml", ""), "hor_event_deletion",
                       monomer_length_bp=L, k=k_var, n_copies=base_n,
                       divergence=0.15, mutation_rate=0.02, indel_rate=0.0,
                       event_type="DELETION", start_copy=st, length_copies=ln)
                n += 1
    return n


# ----------------------------------------------------------------------
# 07: mixed (coexisting periods) (100)
# ----------------------------------------------------------------------

def category_07_mixed() -> int:
    n = 0
    # 5 (L_a, k_a) × 5 (L_b, k_b) × 4 nc_pairs = 100
    combos_a = [(100, 4), (171, 6), (171, 12), (200, 8), (250, 6)]
    combos_b = [(120, 4), (200, 8), (250, 10), (300, 6), (400, 4)]
    nc_pairs = [(50, 50), (100, 50), (50, 100), (100, 100)]
    for (La, ka) in combos_a:
        for (Lb, kb) in combos_b:
            for (nca, ncb) in nc_pairs:
                name = (f"mx_a{La:03d}-{ka:02d}_b{Lb:03d}-{kb:02d}"
                        f"_n{nca:03d}-{ncb:03d}.yaml")
                body = f"""
global:
  mutation_rate: 0.02
templates:
  alpha:
    type: HOR_slots
    monomer_length_bp: {La}
    k: {ka}
    inter_slot_divergence: 0.15
  beta:
    type: HOR_slots
    monomer_length_bp: {Lb}
    k: {kb}
    inter_slot_divergence: 0.15
structure:
  - type: HOR
    template: alpha
    n_copies: {nca}
  - type: HOR
    template: beta
    n_copies: {ncb}
"""
                write_yaml(CFG_ROOT / "07_mixed" / name, body)
                record(name.replace(".yaml", ""), "mixed",
                       monomer_length_bp=La, k=ka, n_copies=nca + ncb,
                       divergence=0.15, mutation_rate=0.02, indel_rate=0.0,
                       second_monomer_bp=Lb, second_k=kb)
                n += 1
    return n


# ----------------------------------------------------------------------
# 08: random (negative control) (50)
# ----------------------------------------------------------------------

def category_08_random() -> int:
    n = 0
    lengths = [5_000, 10_000, 20_000, 50_000, 100_000]
    kinds   = ["random", "AT_rich", "GC_rich", "retro_like", "segdup_like"]
    for ln in lengths:
        for k in kinds:
            for rep in (0, 1):
                name = f"rn_L{ln:06d}_{k}_r{rep}.yaml"
                # segdup_like with no preceding bytes falls back to random.
                body = f"""
templates:
  dummy:
    type: monomer
    monomer_length_bp: 1
structure:
  - type: INSERTION
    length_bp: {ln}
    kind: {k}
"""
                write_yaml(CFG_ROOT / "08_random" / name, body)
                record(name.replace(".yaml", ""), "random",
                       monomer_length_bp=0, k=0, n_copies=0,
                       divergence=0.0, mutation_rate=0.0, indel_rate=0.0,
                       insertion_kind=k, insertion_length_bp=ln)
                n += 1
    return n


# ----------------------------------------------------------------------
# 09: gc_bias (AT-rich / GC-rich simple TRs) (50)
# ----------------------------------------------------------------------

def category_09_gc_bias() -> int:
    n = 0
    gcs       = [0.1, 0.2, 0.3, 0.7, 0.8]       # 5 GC values
    n_copies  = [500, 1000]                     # 2
    monomers  = [100, 170, 250, 500, 1000]      # 5 → 5×2×5 = 50
    combos = list(itertools.product(gcs, n_copies, monomers))
    for (gc, nc, L) in combos:
        name = f"gc_g{int(gc*10):02d}_L{L:04d}_n{nc:04d}.yaml"
        body = f"""
global:
  mutation_rate: 0.02
templates:
  m:
    type: monomer
    monomer_length_bp: {L}
    gc_content: {gc}
structure:
  - type: SIMPLE_TR
    template: m
    n_copies: {nc}
"""
        write_yaml(CFG_ROOT / "09_gc_bias" / name, body)
        record(name.replace(".yaml", ""), "gc_bias",
               monomer_length_bp=L, k=1, n_copies=nc, divergence=0.0,
               mutation_rate=0.02, indel_rate=0.0, gc_content=gc)
        n += 1
    return n


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

def main() -> None:
    counts = {
        "01_simple_tr": category_01_simple_tr(),
        "02_hor_clean": category_02_hor_clean(),
        "03_hor_wobble": category_03_hor_wobble(),
        "04_hor_shift": category_04_hor_shift(),
        "05_hor_insertion": category_05_hor_insertion(),
        "06_hor_event": category_06_hor_event(),
        "07_mixed": category_07_mixed(),
        "08_random": category_08_random(),
        "09_gc_bias": category_09_gc_bias(),
    }
    total = sum(counts.values())
    for cat, n in counts.items():
        print(f"{cat:20s} {n:5d}")
    print(f"{'TOTAL':20s} {total:5d}")

    # Manifest
    all_keys = set()
    for r in MANIFEST_ROWS:
        all_keys.update(r.keys())
    header_order = [
        "case_id", "category", "monomer_length_bp", "k", "n_copies",
        "divergence", "mutation_rate", "indel_rate",
    ]
    other = sorted(k for k in all_keys if k not in header_order)
    cols = header_order + other
    manifest_path = ROOT / "manifest.tsv"
    with manifest_path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in MANIFEST_ROWS:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"\nwrote {manifest_path}")


if __name__ == "__main__":
    main()

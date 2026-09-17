#!/usr/bin/env python3
"""Compare transported vs reconstructed persistent homology across L0-L3."""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from codec import EscDecoder, EscEncoder, strip_to_mode
from gray_scott import PRESETS, GrayScott
from homology import compare_homology
from manifold import attach_scaler, get_or_build_manifold

SEED, SIZE, STEPS = 42, 128, 900
CACHE, OUT = ROOT / "reports" / "manifold_s128.npz", ROOT / "reports"
CASES = [
    ("mitosis_clean", "mitosis", []),
    ("mitosis_three_defects", "mitosis", [(0.28, 0.30, 0.05, 0.9), (0.72, 0.68, 0.04, 0.85), (0.55, 0.22, 0.035, 0.7)]),
    ("mazes_clean", "mazes", []),
    ("worms_one_defect", "worms", [(0.33, 0.62, 0.045, 0.8)]),
    ("holes_one_defect", "holes", [(0.4, 0.4, 0.05, 0.75)]),
    ("coral_one_defect", "coral", [(0.7, 0.35, 0.04, 0.8)]),
]

def source(preset, defects):
    sim = GrayScott(SIZE, SIZE, PRESETS[preset])
    sim.reset((0.5, 0.5), noise=0.02, rng_seed=SEED)
    sim.step_n(max(0, STEPS - 50))
    for nx, ny, r, s in defects:
        sim.seed(nx, ny, r, s)
    sim.step_n(min(50, STEPS))
    return sim

def main():
    m = attach_scaler(get_or_build_manifold(CACHE, size=128, steps=700, snapshot_every=40, burn_in=280, seed=SEED))
    enc, dec = EscEncoder(m, size=SIZE, seed=SEED), EscDecoder(m)
    rows = []
    for name, preset, defects in CASES:
        src = source(preset, defects)
        encoded = enc.encode(src, time_step=STEPS, declared_class=preset)
        wire = encoded.packet["extropic_residual"]["homology"]
        print(f"=== {name} wire beta0={wire['beta0']} ===")
        for mode in ("L0", "L1", "L2", "L3"):
            rec = dec.decode(strip_to_mode(encoded.packet, mode))
            hr = compare_homology(src, rec)
            rows.append({"case": name, "mode": mode, "wire_beta0": wire["beta0"],
                         "src_beta0": hr.beta0_working[0], "rec_beta0": hr.beta0_working[1],
                         "beta0_match": int(hr.beta0_match), "bottleneck_h0": hr.bottleneck_h0})
            print(f"  {mode}: {hr.beta0_working} match={hr.beta0_match} dB={hr.bottleneck_h0:.4f}")
    OUT.mkdir(exist_ok=True)
    path = OUT / "homology_compare.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("wrote", path)

if __name__ == "__main__":
    main()

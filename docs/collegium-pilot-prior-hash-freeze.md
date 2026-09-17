# Collegium Pilot — prior_hash freeze stamp
**Season:** Extropy / ESC (C02 vs C06 Pilot)
**Frozen:** 2026-09-17 (America/New_York)
**Stamped by:** R3DB0T

## Packet
- File: Extropic Semantic Codecs (20 pp)
- sha256: `eb5372fc297b4b88d0a7e41a46bb5e00d860121abd2be07768ba3f093a14eb9f`
- Engine repo: `IAmM3ta/MORPH-OS` (hyphen) — **not** `MORPHOS`
- Engine HEAD: `bd40f4fdd8840681f6fa4baf27d500e81312979c`

## Shared prior (esc/SPEC.md)
- Engine id: `morphogen-gs-cpu-v2`
- Formula: `prior_hash = sha256(engine|name|F|k|Du|Dv)[:16]` with `:.6f` floats and **display name** (e.g. `Mitosis`)
- Kernel: periodic 5-point Laplacian, `dt=1`, clamp `[0,1]`, canonical square + 2% speckles

## Blob anchors
| Path | blob sha |
|------|----------|
| esc/gray_scott.py | c4ab2c2f4d4ad6564429a45a03904b710373a212 |
| esc/manifold.py | 325a8f2e7d171c6174bd2a097bc307eee2cfbc88 |
| esc/homology.py | e26c67c7b0bef5108300444a8b61bca73474d22f |

## prior_hash table (Pearson presets)

| class_id | prior_hash | payload |
|----------|------------|---------|
| mitosis | `40c8d4a53611fa4b` | morphogen-gs-cpu-v2\|Mitosis\|0.036700\|0.064900\|0.160000\|0.080000 |
| solitons | `918a0e61ba9e581d` | morphogen-gs-cpu-v2\|Solitons\|0.035300\|0.065300\|0.160000\|0.080000 |
| pulsing | `d2ef378ff90b5ead` | morphogen-gs-cpu-v2\|Pulsing\|0.025000\|0.060000\|0.140000\|0.070000 |
| holes | `96cbfe45ca1f1e50` | morphogen-gs-cpu-v2\|Holes\|0.039000\|0.058000\|0.160000\|0.080000 |
| mazes | `ecdca0f457a4b85c` | morphogen-gs-cpu-v2\|Mazes\|0.029000\|0.057000\|0.160000\|0.080000 |
| fingerprint | `6ae5164ac2bc2bad` | morphogen-gs-cpu-v2\|Fingerprint\|0.026000\|0.061000\|0.160000\|0.080000 |
| spirals | `33040ef2c0387c96` | morphogen-gs-cpu-v2\|Spirals\|0.018000\|0.051000\|0.160000\|0.080000 |
| worms | `24fabdd4fec7e363` | morphogen-gs-cpu-v2\|Worms\|0.046000\|0.063000\|0.160000\|0.080000 |
| coral | `10611f7fcc950bab` | morphogen-gs-cpu-v2\|Coral\|0.054500\|0.062000\|0.160000\|0.080000 |
| uskate | `87b066b9f4b0394a` | morphogen-gs-cpu-v2\|Skate\|0.062000\|0.060900\|0.160000\|0.080000 |

## Pilot cut anchors
- Mitosis prior_hash for L0/L1 homology cut: **`40c8d4a53611fa4b`**
- Cuts: LEI ≠ N; L1 vs L0 on mitosis+disks β₀

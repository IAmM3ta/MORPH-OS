# Camera-ready package — Collegium beta-03 (Extropy / ESC)

**Tag purpose:** workshop/arXiv priority artifact for MORPH L1 mutation β₀ recovery + homology ladder.  
**Board report:** Collegium `BOARD-REPORT-novelty-publishability.md` §4 / §8 (R3DB0T ops cut).

## Freeze locators (do not soft-extend)

| Item | Value |
|------|-------|
| Packet | *Extropic Semantic Codecs*, **20 pp** |
| Packet sha256 | `eb5372fc297b4b88d0a7e41a46bb5e00d860121abd2be07768ba3f093a14eb9f` |
| Engine | `morphogen-gs-cpu-v2` |
| Engine repo | `IAmM3ta/MORPH-OS` (hyphen) — **not** `MORPHOS` |
| Mitosis prior_hash | `40c8d4a53611fa4b` |
| Stamp doc | [`docs/collegium-pilot-prior-hash-freeze.md`](docs/collegium-pilot-prior-hash-freeze.md) |

**Excluded from this tag:** 32 pp editions (`b6b334a0…`, `73ead43c…`); V2 FPR/preservation % as evaluated science; IMAGE_8 / `MORPHOS` codec repo.

## Recompute (Table 4.4 / homology ladder)

```bash
cd esc
pip install -r requirements.txt
python run_homology.py
# writes reports/homology_compare.csv — includes mitosis_three_defects L0–L3
```

Supporting code: `esc/homology.py`, `esc/codec.py`, `esc/gray_scott.py`, `esc/manifold.py`.

## Claim surface for this tag

Publishable *as engineering contribution* (conditional on these scripts + reports shipping):

1. L0 miss vs L1 match on mitosis+3 disks working-β₀
2. Homology ladder tradeoffs (β₀ vs H₀ bottleneck)
3. Codec ≠ encryption hygiene (recipe / pointer + residual, not a lock)

Not claimed here: Extropy metaphysics; LEI≠N as frontier novelty; class accuracy 1.00 without holdout.

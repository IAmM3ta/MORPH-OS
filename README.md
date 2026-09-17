# MORPH-OS

**Extropic Semantic Codec for Morphogen Gray–Scott fields.**

Ship an ontology pointer and an extropic residual. Reconstruct the field locally. Preserve structural ontology and persistent homology — not pixel parity.

This repository is the implementable specification that sits between:

- [MORPHOS](https://github.com/IAmM3ta/MORPHOS) — generative image codec + Morphogen v2 instrument
- Thomas, *Extropy as a Metric of Cognitive Novelty* (2025)
- *Extropic Semantic Codecs* (ESC 1.0)

## Claim, narrowed

We do **not** transmit lossless bitstreams. We transmit:

1. `GrayScott.TuringField` + Pearson preset + `prior_hash`
2. Manifold anchor `(feed, kill, du, dv, seed, size, time_step)`
3. Extropy scalars N, LEI, C, L, Ex
4. Localized mutations and/or a coarse concentration-chart residual
5. A superlevel **H0 barcode** on the periodic torus

The decoder already has the engine. It instantiates the class, integrates, applies the residual, and checks the reconstructed barcode against the wire barcode.

## Residual ladder

| Mode | Wire | When |
|---|---|---|
| L0 | pointer + timestep | W1 below gate — on-manifold drift |
| L1 | + disks | default live packet — **beta0 exact** on the lab suite |
| L2 | + 8x8 dV-bar | lower H0 bottleneck, may split colonies |
| L3 | + 16x16 dV-bar | finer chart, still not pixels |

Lab, 128² / 900 steps (see `esc/reports/`):

- L0 ~ 412 B vs ~ 11 KB PNG (0.09x)
- L1 beta0 accuracy **1.00**, blob-count error **0**
- L2 mean H0 bottleneck **0.048** (best chart match, not best beta0)
- 256² mazes L0: **409 B** vs **60 KB** PNG, 121	o121 components

## Persistent homologous structures

`esc/homology.py` is a dependency-free cubical superlevel persistence engine:

- H0 bars via union-find on the periodic grid
- Betti profile (beta0, beta0 of complement) at five thresholds
- Bottleneck-style matching on the top bars
- Barcode rides on `extropic_residual.homology`

L1 is the mode that keeps working-threshold beta0. L2/L3 shrink barcode distance and can change the 0.25 isosurface. That tension is the point of the ladder.

## Layout

```
esc/                  Python reference codec
  gray_scott.py       Morphogen-aligned CPU kernel
  manifold.py         shared M, prior_hash, persisted catalog
  codec.py            ESC 1.0 encode / decode / Sync ingest
  homology.py         H0 persistence + Betti profile
  morphogen_esc.ts    L0 mapper for the Morphogen 20 Hz pump
  SPEC.md             closed equations
  reports/            rate-structure and homology CSVs
```

## Run

```bash
cd esc
pip install -r requirements.txt
python run_sweep.py
python run_homology.py
python encode_sync.py packets/incoming_sync_mitosis.json
```

## License

MIT

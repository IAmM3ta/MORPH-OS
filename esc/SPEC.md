# ESC 1.0 — Morphogen-native specification

Status: implementable reference. Preserves structural ontology and homology proxies. Does not claim bit parity.

## 1. Shared prior

Encoder and decoder share:

- Engine id `morphogen-gs-cpu-v2` (periodic 5-point Laplacian, `dt = 1`, clamp `[0,1]`)
- Preset table `(feed, kill, du, dv)`
- Canonical reset: `U=1, V=0`, central square `U=0.5, V=1` of half-width `0.08·min(h,w)`, plus 2% speckles from `rng_seed`
- `prior_hash = sha256(engine|name|F|k|Du|Dv)[:16]`

A hash mismatch is a version break, not a soft warning.

## 2. Manifold

Feature vector x in R^21:

`meanU, meanV, energy, cx, cy, edge, V_p10, V_p50, V_p90, hist_12`

Per class C, late-time snapshots form samples. Centroid mu_C, tangent T_p M = top-3 right singular vectors.

M_C = { x : ||(I-P)(x-mu_C)||_2 <= q_0.95 }

Classification readout (not the wire class) uses z-scored Euclidean to centroids.

## 3. Extropy maps

- N(x) = ||(I-P_{T_p M})(x-p)||_2
- LEI(x) = 1-cos theta(x, mu_C)
- C(x) = cos(x, mean samples_C)
- L(x) = constraint scalar on U,V in [0,1] and legal (F,k,Du,Dv)
- Ex(x) = N(x) C(x) L(x)
- Delta_u = 1/2 (W1(H^U_t, H^U_{t-1}) + W1(H^V_t, H^V_{t-1}))

## 4. Residual ladder

| Mode | Payload | Decode |
|---|---|---|
| L0 | class + anchor + time_step | exact prior replay |
| L1 | + localized disks | seed disks on the replay |
| L2 | + 8x8 Delta V-bar vs same-step reference | upsample and add the delta |
| L3 | + 16x16 Delta V-bar | same, finer chart |

L2/L3 are chart residuals, not alpha-blends. Clean L0 must remain numerically exact.

## 5. Packet

Required keys: esc_version, ontology_pointer, manifold_anchor, extropic_residual, temporal_delta. Optional residual_mode.

## 6. Metrics that count

- Packet-declared class accuracy (ontology transport)
- Working-threshold beta0 match and H0 bottleneck (homology)
- W1 on concentration histograms (temporal integrity)
- Payload bytes vs PNG / raw grid16

Pixel MSE and edge-energy relative error are diagnostics.

## 7. Persistent homologous structures

Superlevel H0 on the periodic torus (union-find, pixels appear high-V first):

- bar = (birth, death); death = 0 if the component survives to V = 0
- signature = top-8 bars by persistence
- Betti profile = (beta0, beta0 of complement) at t in {0.70, 0.50, 0.35, 0.25, 0.15}
- working beta0 = periodic components of {V >= 0.25}

These ride on extropic_residual.homology. Decode success for L1 is working-beta0 match. L2/L3 are scored on H0 bottleneck against the wire barcode.

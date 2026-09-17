# Chapter IV (excerpt): §§4.2–4.4

*Rewritten against the MORPH-OS 1.0 reference stack and laboratory reports.
A MORPH (Minimal Ontic Reference, Preserved Homology) is an ESC packet whose residual is sufficient to replay a Gray–Scott field up to working-threshold β0 on the periodic torus.*

See the companion Word excerpt in the local artifacts drop, or this file.

## 4.2 Mathematical Formalization of the RD-Specific Codec

To render the Extropic Semantic Codec computable inside a generative reaction–diffusion paradigm, the maps that were left schematic in the theoretical chapters are here given as closed operations on a Gray–Scott field. Encoder and decoder are assumed to share a prior: the engine identifier `morphogen-gs-cpu-v2` (periodic five-point Laplacian, unit Euler step, pointwise clamp of $U,V$ to $[0,1]$), a Pearson preset table $(F,k,D_U,D_V)$, and a canonical reset (central square of half-width $0.08\cdot\min(h,w)$ with $U=0.5$, $V=1$, plus two-percent deterministic speckles from `rng_seed`). The integrity of that agreement is a sixteen-hex prefix

$$\texttt{prior\_hash}=\mathrm{sha256}(\texttt{engine}\mid\texttt{name}\mid F\mid k\mid D_U\mid D_V)[:16].$$

A hash mismatch is a version break, not a soft warning.

### 4.2.1 The shared manifold $\mathcal{M}$

Let $\mathcal{C}$ be an agreed object ontology (e.g. `GrayScott.TuringField` with Pearson species `mitosis`, `mazes`, `worms`, `holes`, `coral`). From an ensemble of late-time snapshots of class $\mathcal{C}$ one extracts a twenty-one-dimensional feature

$$x=\bigl(\overline{U},\,\overline{V},\,E,\,c_x,\,c_y,\,e,\,V_{10},\,V_{50},\,V_{90},\,H_{12}\bigr)\in\mathbb{R}^{21},$$

where $E$ is field energy, $(c_x,c_y)$ is the concentration centroid, $e$ is Sobel edge energy, $V_{p}$ are percentiles of $V$, and $H_{12}$ is a twelve-bin histogram of $V$. Write $\mu_{\mathcal{C}}$ for the class centroid and $T_p\mathcal{M}$ for the three-dimensional tangent spanned by the leading right singular vectors of the centred ensemble. With $P_{T_p\mathcal{M}}$ the orthogonal projector onto that tangent, the probabilistic manifold of the class is the sublevel of residual radius at the empirical 95th percentile:

$$\mathcal{M}_{\mathcal{C}}=\bigl\{x:\ \bigl\|(I-P_{T_p\mathcal{M}})(x-\mu_{\mathcal{C}})\bigr\|_2\le q_{0.95}(\mathcal{C})\bigr\}.$$

Classification *readout* (nearest z-scored centroid) is a laboratory diagnostic. It is **not** the class that travels on the wire. The transported class is the ontology pointer. On the MORPH-OS sweep the packet-declared class is recovered at accuracy $1.00$; z-scored nearest-neighbour readout of the source feature is only $0.73$, because several Pearson species overlap in $\mathbb{R}^{21}$. Shipping the pointer is therefore load-bearing, not decorative.

### 4.2.2 Orthogonal novelty, LEI, coherence, integrity

$$N(x)=\bigl\|(I-P_{T_p\mathcal{M}})(x-p)\bigr\|_2,\qquad \mathrm{LEI}(x)=1-\cos\theta(x,\mu_{\mathcal{C}}).$$

$$C(x)=\cos\bigl(x,\ \mathrm{mean}(\{x_i\}_{\mathcal{C}})\bigr),\qquad L(x)=\mathbf{1}_{U,V\in[0,1]}\cdot\mathbf{1}_{(F,k,D_U,D_V)\ \mathrm{legal}}.$$

$$\mathrm{Ex}(x)=N(x)\cdot C(x)\cdot L(x).$$

### 4.2.3 Transformational delta $\Delta_u$

$$\Delta_u(t)=\tfrac12\Bigl(W_1(H^U_t,H^U_{t-1})+W_1(H^V_t,H^V_{t-1})\Bigr).$$

$W_1$ is computed on bin centres with histogram masses as weights.

### 4.2.4 Persistent homologous structures

Superlevel $H_0$ of $V$ on the periodic torus $\mathbb{T}^2$. Working threshold $\tau=0.25$:

$$\beta_0^{\mathrm{work}}(V)=\beta_0\bigl(\{V\ge 0.25\}\bigr)\quad\text{on }\mathbb{T}^2.$$

Engine: `esc/homology.py`.

### 4.2.5 Residual ladder

| Mode | On the wire | Decoder |
|---|---|---|
| L0 | ontology pointer + manifold anchor + `time_step` | exact prior replay |
| L1 | L0 + localized disks in normalized coordinates | seed disks on the replay, then integrate |
| L2 | + an $8\times 8$ chart of $\overline{\Delta V}$ | upsample and *add* the chart |
| L3 | as L2 with a $16\times 16$ chart | same, finer chart |

### 4.2.6 Packet

Required keys of ESC 1.0: `esc_version`, `ontology_pointer`, `manifold_anchor`, `extropic_residual`, `temporal_delta`.

## 4.3 Reference Implementation

MORPH-OS is the implementable specification. The short L1 listing in the Word excerpt is pedagogical only. Homology, manifold charts, L2/L3, and Sync ingest live in `esc/`.

## 4.4 Empirical Evaluation

See `esc/reports/sweep_summary.json` and `esc/reports/homology_compare.csv`.

L1 working-$\beta_0$ accuracy **1.00**. L2 mean H0 bottleneck **0.048**. L0 ~412 B vs ~11 KB PNG at $128^2$. Packet-declared class accuracy **1.00**; z-scored nearest neighbour **0.73**.

Full prose: repository companion Word file generated from this rewrite, and the conversation artefact `Chapter_IV_4.2-4.4_MORPH-OS.docx`.

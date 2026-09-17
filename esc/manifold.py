"""Shared probabilistic manifold for Morphogen Gray–Scott species."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from numpy.typing import NDArray

from gray_scott import PRESETS, GrayScott, RdParams

ENGINE_ID = "morphogen-gs-cpu-v2"
N_BINS = 12
FEATURE_NAMES = [
    "mean_u", "mean_v", "energy", "cx", "cy", "edge", "v_p10", "v_p50", "v_p90",
] + [f"hist_{i}" for i in range(N_BINS)]


def prior_hash(params: RdParams) -> str:
    payload = (
        f"{ENGINE_ID}|{params.name}|{params.feed:.6f}|{params.kill:.6f}|"
        f"{params.du:.6f}|{params.dv:.6f}"
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def feature_vector(sim: GrayScott) -> NDArray[np.float64]:
    st = sim.stats()
    V = sim.V
    hist, _ = np.histogram(V.ravel(), bins=N_BINS, range=(0.0, 1.0), density=True)
    hist = hist / (hist.sum() + 1e-12)
    return np.array(
        [st.mean_u, st.mean_v, st.energy, st.cx, st.cy, st.edge,
         float(np.percentile(V, 10)), float(np.percentile(V, 50)), float(np.percentile(V, 90)),
         *hist.tolist()],
        dtype=np.float64,
    )


def v_histogram(sim: GrayScott, bins: int = 32):
    hist, edges = np.histogram(sim.V.ravel(), bins=bins, range=(0.0, 1.0), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mass = hist * np.diff(edges)
    mass = mass / (mass.sum() + 1e-12)
    return centers, mass


@dataclass
class ClassManifold:
    class_id: str
    params: RdParams
    prior: str
    centroid: NDArray[np.float64]
    basis: NDArray[np.float64]
    samples: NDArray[np.float64]
    tau: float

    def orthogonal_residual(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        d = x - self.centroid
        if self.basis.size == 0:
            return d
        proj = self.basis @ (self.basis.T @ d)
        return d - proj

    def orthogonal_residual_norm(self, xs: NDArray[np.float64]) -> NDArray[np.float64]:
        d = xs - self.centroid
        if self.basis.size == 0:
            return np.linalg.norm(d, axis=1)
        proj = (d @ self.basis) @ self.basis.T
        return np.linalg.norm(d - proj, axis=1)

    def in_manifold(self, x: NDArray[np.float64]) -> bool:
        residuals = self.orthogonal_residual_norm(self.samples)
        thresh = float(np.quantile(residuals, self.tau)) if len(residuals) else 0.0
        return self.orthogonal_residual_norm(x[None, :])[0] <= max(thresh, 1e-6)


class SharedManifold:
    def __init__(self, classes: Dict[str, ClassManifold], tau: float = 0.95):
        self.classes = classes
        self.tau = tau

    def nearest(self, x: NDArray[np.float64]):
        best_id, best, best_dist = "", None, float("inf")
        for cid, cm in self.classes.items():
            dist = float(np.linalg.norm(x - cm.centroid))
            if dist < best_dist:
                best_id, best, best_dist = cid, cm, dist
        assert best is not None
        return best_id, best, best_dist


def build_manifold(size: int = 96, steps: int = 420, snapshot_every: int = 30,
                   burn_in: int = 180, seed: int = 42, tau: float = 0.95,
                   tangent_rank: int = 3) -> SharedManifold:
    rng = np.random.default_rng(seed)
    classes: Dict[str, ClassManifold] = {}
    for cid, params in PRESETS.items():
        sim = GrayScott(size, size, params)
        sim.reset((0.5, 0.5), radius=0.04)
        rows: List[NDArray[np.float64]] = []
        for t in range(1, steps + 1):
            sim.step()
            if t >= burn_in and t % snapshot_every == 0:
                rows.append(feature_vector(sim))
        for _ in range(4):
            sim.reset((0.5, 0.5), radius=0.04)
            sim.seed(float(rng.uniform(0.3, 0.7)), float(rng.uniform(0.3, 0.7)), 0.03, 0.8)
            sim.step_n(burn_in)
            rows.append(feature_vector(sim))
        samples = np.stack(rows, axis=0)
        centroid = samples.mean(axis=0)
        centered = samples - centroid
        try:
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            rank = min(tangent_rank, vt.shape[0])
            basis = vt[:rank].T
        except np.linalg.LinAlgError:
            basis = np.zeros((samples.shape[1], 0))
        classes[cid] = ClassManifold(cid, params, prior_hash(params), centroid, basis, samples, tau)
    return SharedManifold(classes, tau=tau)


def standardize_stats(manifold: "SharedManifold"):
    xs = np.concatenate([cm.samples for cm in manifold.classes.values()], axis=0)
    mu = xs.mean(axis=0)
    sd = xs.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return mu, sd


def attach_scaler(manifold: "SharedManifold") -> "SharedManifold":
    manifold.mu, manifold.sd = standardize_stats(manifold)
    return manifold


def scaled_nearest(manifold: "SharedManifold", x: np.ndarray):
    mu = getattr(manifold, "mu", None)
    sd = getattr(manifold, "sd", None)
    if mu is None:
        attach_scaler(manifold)
        mu, sd = manifold.mu, manifold.sd
    xz = (x - mu) / sd
    best_id, best, best_d = "", None, float("inf")
    for cid, cm in manifold.classes.items():
        d = float(np.linalg.norm(xz - (cm.centroid - mu) / sd))
        if d < best_d:
            best_id, best, best_d = cid, cm, d
    assert best is not None
    return best_id, best, best_d


def save_manifold(manifold: "SharedManifold", path) -> None:
    payload = {"tau": np.array([manifold.tau]), "class_ids": np.array(list(manifold.classes.keys()))}
    for cid, cm in manifold.classes.items():
        payload[f"{cid}__centroid"] = cm.centroid
        payload[f"{cid}__basis"] = cm.basis
        payload[f"{cid}__samples"] = cm.samples
        payload[f"{cid}__prior"] = np.array(cm.prior)
        payload[f"{cid}__params"] = np.array([cm.params.feed, cm.params.kill, cm.params.du, cm.params.dv])
    np.savez_compressed(path, **payload)


def load_manifold(path) -> "SharedManifold":
    z = np.load(path, allow_pickle=True)
    tau = float(z["tau"][0])
    classes = {}
    for cid in z["class_ids"]:
        cid = str(cid)
        classes[cid] = ClassManifold(
            cid, PRESETS[cid], str(z[f"{cid}__prior"]),
            z[f"{cid}__centroid"], z[f"{cid}__basis"], z[f"{cid}__samples"], tau,
        )
    return attach_scaler(SharedManifold(classes, tau=tau))


def get_or_build_manifold(path, **kwargs) -> "SharedManifold":
    from pathlib import Path
    path = Path(path)
    if path.exists():
        return load_manifold(path)
    m = attach_scaler(build_manifold(**kwargs))
    path.parent.mkdir(parents=True, exist_ok=True)
    save_manifold(m, path)
    return m

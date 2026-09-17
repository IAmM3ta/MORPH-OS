"""Extropic Semantic Codec — Morphogen-native reference encoder / decoder."""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray
from scipy import ndimage
from scipy.stats import wasserstein_distance

from gray_scott import DIFF_BOUNDS, FEED_BOUNDS, KILL_BOUNDS, PRESETS, GrayScott, RdParams
from homology import attach_homology_to_packet
from manifold import SharedManifold, feature_vector, prior_hash, v_histogram

ESC_VERSION = "1.0"
ONTOLOGY_CLASS = "GrayScott.TuringField"


def cosine(a, b) -> float:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def lei(x, centroid) -> float:
    return float(max(0.0, 1.0 - cosine(x, centroid)))


def orthogonal_novelty(x, residual) -> float:
    return float(np.linalg.norm(residual))


def contextual_coherence(x, neighborhood) -> float:
    return float(max(0.0, min(1.0, cosine(x, neighborhood.mean(axis=0)))))


def logical_integrity(sim: GrayScott) -> float:
    score = 1.0
    U, V = sim.U, sim.V
    if not np.isfinite(U).all() or not np.isfinite(V).all():
        return 0.0
    if (U < -1e-9).any() or (V < -1e-9).any() or (U > 1.0 + 1e-6).any() or (V > 1.0 + 1e-6).any():
        score *= 0.2
    p = sim.params
    if not (FEED_BOUNDS[0] <= p.feed <= FEED_BOUNDS[1]):
        score *= 0.5
    if not (KILL_BOUNDS[0] <= p.kill <= KILL_BOUNDS[1]):
        score *= 0.5
    if not (DIFF_BOUNDS[0] <= p.du <= DIFF_BOUNDS[1] and DIFF_BOUNDS[0] <= p.dv <= DIFF_BOUNDS[1]):
        score *= 0.5
    if p.dv > p.du:
        score *= 0.7
    return float(max(0.0, min(1.0, score)))


def v_histogram_of(arr, bins: int = 32):
    hist, edges = np.histogram(arr.ravel(), bins=bins, range=(0.0, 1.0), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mass = hist * np.diff(edges)
    mass = mass / (mass.sum() + 1e-12)
    return centers, mass


def wasserstein_uv(prev: GrayScott, curr: GrayScott) -> float:
    c0, m0 = v_histogram(prev)
    c1, m1 = v_histogram(curr)
    w_v = float(wasserstein_distance(c0, c1, m0, m1))
    cu0, mu0 = v_histogram_of(prev.U)
    cu1, mu1 = v_histogram_of(curr.U)
    return 0.5 * (w_v + float(wasserstein_distance(cu0, cu1, mu0, mu1)))


def extract_mutations(observed: GrayScott, reference: GrayScott, max_defects: int = 4,
                      threshold: float = 0.18, min_area: int = 12):
    delta = np.abs(observed.V - reference.V)
    mask = delta >= threshold
    labeled, nlab = ndimage.label(mask)
    mutations = []
    h, w = delta.shape
    for lab in range(1, nlab + 1):
        comp = labeled == lab
        area = int(comp.sum())
        if area < min_area:
            continue
        ys, xs = np.where(comp)
        mass = delta[comp]
        mag = float(mass.mean())
        mutations.append({
            "type": "localized_defect",
            "cx": round(float(np.average(xs, weights=mass) / w), 4),
            "cy": round(float(np.average(ys, weights=mass) / h), 4),
            "magnitude": round(mag, 4),
            "radius": round(max(float(np.sqrt(area / np.pi) / min(w, h)), 0.02), 4),
            "area_px": area,
        })
    mutations.sort(key=lambda m: m["magnitude"] * m["area_px"], reverse=True)
    return mutations[:max_defects]


def apply_mutations(sim: GrayScott, mutations) -> None:
    for m in mutations:
        if m.get("type") != "localized_defect":
            continue
        sim.seed(float(m["cx"]), float(m["cy"]), radius=float(m.get("radius", 0.04)),
                 strength=float(m.get("magnitude", 0.4)))


@dataclass
class EncodeResult:
    packet: dict[str, Any]
    class_id: str
    features: NDArray[np.float64]
    novelty: float
    lei_score: float
    c_score: float
    l_score: float


class EscEncoder:
    def __init__(self, manifold: SharedManifold, size: int, seed: int = 42):
        self.manifold = manifold
        self.size = size
        self.seed = seed

    def _clean_reference(self, class_id: str, steps: int) -> GrayScott:
        ref = GrayScott(self.size, self.size, PRESETS[class_id])
        ref.reset((0.5, 0.5), radius=0.04)
        ref.step_n(steps)
        return ref

    def encode(self, sim: GrayScott, time_step: int, prev: Optional[GrayScott] = None,
               declared_class: Optional[str] = None, include_moments: bool = True) -> EncodeResult:
        x = feature_vector(sim)
        if declared_class and declared_class in self.manifold.classes:
            class_id, cm = declared_class, self.manifold.classes[declared_class]
        else:
            class_id, cm, _ = self.manifold.nearest(x)
        residual = cm.orthogonal_residual(x)
        n_true = orthogonal_novelty(x, residual)
        n_lei = lei(x, cm.centroid)
        c_score = contextual_coherence(x, cm.samples)
        l_score = logical_integrity(sim)
        ref = self._clean_reference(class_id, time_step)
        mutations = extract_mutations(sim, ref)
        moments = sim.moments_grid(8)
        moments_delta = (moments - ref.moments_grid(8)).round(4)
        grid16_delta = (sim.grid16() - ref.grid16()).round(4)
        packet = {
            "esc_version": ESC_VERSION,
            "ontology_pointer": {"class": ONTOLOGY_CLASS, "preset": class_id, "prior_hash": cm.prior},
            "manifold_anchor": {
                "feed": cm.params.feed, "kill": cm.params.kill,
                "du": cm.params.du, "dv": cm.params.dv,
                "seed": self.seed, "size": self.size,
            },
            "extropic_residual": {
                "N_true": round(n_true, 6), "N_proxy_lei": round(n_lei, 6),
                "C_score": round(c_score, 6), "L_score": round(l_score, 6),
                "Ex": round(n_true * c_score * l_score, 6), "mutations": mutations,
            },
            "temporal_delta": {
                "time_step": time_step,
                "delta_w1": None if prev is None else round(wasserstein_uv(prev, sim), 6),
                "moments_grid": moments.round(4).tolist(),
                "moments_delta": moments_delta.tolist(),
                "grid16_delta": grid16_delta.tolist(),
            },
        }
        attach_homology_to_packet(packet, sim)
        return EncodeResult(packet, class_id, x, n_true, n_lei, c_score, l_score)


class EscDecoder:
    def __init__(self, manifold: SharedManifold):
        self.manifold = manifold

    def decode(self, packet: dict[str, Any], settle_steps: int = 8) -> GrayScott:
        if packet.get("esc_version") != ESC_VERSION:
            raise ValueError(f"unsupported esc_version {packet.get('esc_version')}")
        pointer = packet["ontology_pointer"]
        if pointer["class"] != ONTOLOGY_CLASS:
            raise ValueError(f"unknown ontology class {pointer['class']}")
        preset = pointer["preset"]
        if preset not in PRESETS:
            raise ValueError(f"unknown preset {preset}")
        anchor = packet["manifold_anchor"]
        size = int(anchor.get("size", 96))
        params = RdParams(
            feed=float(anchor["feed"]), kill=float(anchor["kill"]),
            du=float(anchor.get("du", PRESETS[preset].du)),
            dv=float(anchor.get("dv", PRESETS[preset].dv)),
            name=PRESETS[preset].name,
        )
        sim = GrayScott(size, size, params)
        sim.reset((0.5, 0.5), radius=0.04, noise=0.02, rng_seed=int(anchor.get("seed", 42)))
        sim.step_n(int(packet["temporal_delta"]["time_step"]))
        mode = packet.get("residual_mode", "auto")
        mutations = packet["extropic_residual"].get("mutations") or []
        td = packet["temporal_delta"]
        if mode == "L0":
            return sim
        if mode in ("L1", "auto") and mutations:
            apply_mutations(sim, mutations)
            if settle_steps and mode == "L1":
                sim.step_n(settle_steps)
        if mode in ("L2", "auto") and td.get("moments_delta") is not None:
            _apply_chart_delta(sim, np.array(td["moments_delta"], dtype=np.float64))
        if mode == "L3" and td.get("grid16_delta") is not None:
            _apply_chart_delta(sim, np.array(td["grid16_delta"], dtype=np.float64))
        return sim


def _apply_chart_delta(sim: GrayScott, delta: NDArray[np.float64], gain: float = 1.0) -> None:
    if delta.ndim == 1:
        side = int(np.sqrt(delta.size))
        delta = delta.reshape(side, side)
    h, w = sim.height, sim.width
    up = ndimage.zoom(delta, (h / delta.shape[0], w / delta.shape[1]), order=1)
    sim.V[:] = np.clip(sim.V + gain * up, 0.0, 1.0)
    sim.U[:] = np.clip(sim.U - 0.5 * gain * np.maximum(up, 0.0), 0.0, 1.0)


def strip_to_mode(packet: dict[str, Any], mode: str) -> dict[str, Any]:
    p = copy.deepcopy(packet)
    p["residual_mode"] = mode
    td = p["temporal_delta"]
    if mode == "L0":
        p["extropic_residual"]["mutations"] = []
        td.pop("moments_grid", None); td.pop("moments_delta", None); td.pop("grid16_delta", None)
    elif mode == "L1":
        td.pop("moments_grid", None); td.pop("moments_delta", None); td.pop("grid16_delta", None)
    elif mode == "L2":
        td.pop("grid16_delta", None); td.pop("moments_grid", None)
    elif mode == "L3":
        td.pop("moments_grid", None); td.pop("moments_delta", None)
    return p


def encode_from_sync_json(encoder: EscEncoder, sync: dict[str, Any], sim: GrayScott,
                          time_step: int, prev: Optional[GrayScott] = None) -> EncodeResult:
    params = sync.get("params", {})
    declared = None
    if "feed" in params and "kill" in params:
        best, best_d = None, 1e9
        for cid, p in PRESETS.items():
            d = abs(p.feed - params["feed"]) + abs(p.kill - params["kill"])
            if d < best_d:
                best, best_d = cid, d
        if best_d < 1e-4:
            declared = best
    result = encoder.encode(sim, time_step=time_step, prev=prev, declared_class=declared)
    result.packet["source_sync"] = {"t": sync.get("t"), "field": sync.get("field"), "loop": sync.get("loop")}
    return result


def packet_bytes(packet: dict[str, Any]) -> int:
    return len(json.dumps(packet, separators=(",", ":")).encode("utf-8"))

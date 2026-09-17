"""Persistent homologous structures on Gray–Scott fields.

Dependency-free cubical superlevel H0 on the periodic torus.
The barcode rides on the ESC wire so decode can be checked without pixels.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import NDArray
from scipy import ndimage

from gray_scott import GrayScott

THRESHOLDS = (0.70, 0.50, 0.35, 0.25, 0.15)


class UnionFind:
    def __init__(self, n: int):
        self.parent = np.arange(n, dtype=np.int32)
        self.rank = np.zeros(n, dtype=np.int8)
        self.birth = np.full(n, -np.inf, dtype=np.float64)

    def find(self, a: int) -> int:
        p = self.parent
        while p[a] != a:
            p[a] = p[p[a]]
            a = int(p[a])
        return a

    def union(self, a: int, b: int):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return None
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return ra, rb


def _periodic_neighbors(i: int, h: int, w: int):
    y, x = divmod(i, w)
    return [
        ((y - 1) % h) * w + x,
        ((y + 1) % h) * w + x,
        y * w + (x - 1) % w,
        y * w + (x + 1) % w,
    ]


def h0_persistence(field: NDArray[np.float64], min_persist: float = 0.04):
    v = np.ascontiguousarray(field, dtype=np.float64)
    h, w = v.shape
    n = h * w
    order = np.argsort(v.ravel())[::-1]
    uf = UnionFind(n)
    appeared = np.zeros(n, dtype=bool)
    bars = []
    flat = v.ravel()
    for idx in order:
        idx = int(idx)
        val = float(flat[idx])
        appeared[idx] = True
        uf.birth[idx] = val
        older = idx
        seen = []
        for nb in _periodic_neighbors(idx, h, w):
            if not appeared[nb]:
                continue
            root = uf.find(nb)
            if root not in seen:
                seen.append(root)
        for root in seen:
            merged = uf.union(older, root)
            if merged is None:
                continue
            keep, drop = merged
            if uf.birth[keep] < uf.birth[drop]:
                keep, drop = drop, keep
                uf.parent[drop] = keep
            persist = float(uf.birth[drop] - val)
            if persist >= min_persist:
                bars.append((float(uf.birth[drop]), val))
            older = keep
    alive = set(uf.find(i) for i in range(n) if appeared[i])
    for r in alive:
        bars.append((float(uf.birth[r]), 0.0))
    bars.sort(key=lambda b: (b[0] - b[1]), reverse=True)
    return bars


def periodic_components(mask: NDArray[np.bool_]) -> int:
    labeled, n = ndimage.label(mask)
    if n == 0:
        return 0
    h, w = mask.shape
    parent = {i: i for i in range(1, n + 1)}

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for x in range(w):
        a, b = labeled[0, x], labeled[h - 1, x]
        if a and b:
            union(int(a), int(b))
    for y in range(h):
        a, b = labeled[y, 0], labeled[y, w - 1]
        if a and b:
            union(int(a), int(b))
    return len({find(i) for i in range(1, n + 1)})


def betti_profile(field: NDArray[np.float64], thresholds: Iterable[float] = THRESHOLDS):
    out = []
    for t in thresholds:
        mask = field >= t
        out.append({"t": float(t), "beta0": int(periodic_components(mask)), "beta0_complement": int(periodic_components(~mask))})
    return out


def barcode_signature(bars, k: int = 8):
    return [{"birth": round(float(b), 4), "death": round(float(d), 4), "persist": round(float(b - d), 4)} for b, d in bars[:k]]


def bottleneck_h0(a, b, k: int = 12) -> float:
    def pts(bars):
        arr = np.array(bars[:k], dtype=np.float64) if bars else np.zeros((0, 2))
        if len(arr) < k:
            pad = np.zeros((k - len(arr), 2))
            arr = np.vstack([arr, pad]) if len(arr) else pad
        return arr
    return float(np.max(np.abs(pts(a) - pts(b))))


@dataclass
class HomologyReport:
    bars_src: list
    bars_rec: list
    signature_src: list
    signature_rec: list
    profile_src: list
    profile_rec: list
    bottleneck_h0: float
    beta0_working: tuple
    beta0_match: bool


def compare_homology(src: GrayScott, rec: GrayScott, working_t: float = 0.25) -> HomologyReport:
    bs = h0_persistence(src.V)
    br = h0_persistence(rec.V)
    return HomologyReport(
        bars_src=bs,
        bars_rec=br,
        signature_src=barcode_signature(bs),
        signature_rec=barcode_signature(br),
        profile_src=betti_profile(src.V),
        profile_rec=betti_profile(rec.V),
        bottleneck_h0=bottleneck_h0(bs, br),
        beta0_working=(periodic_components(src.V >= working_t), periodic_components(rec.V >= working_t)),
        beta0_match=periodic_components(src.V >= working_t) == periodic_components(rec.V >= working_t),
    )


def attach_homology_to_packet(packet: dict, sim: GrayScott) -> dict:
    bars = h0_persistence(sim.V)
    packet.setdefault("extropic_residual", {})
    packet["extropic_residual"]["homology"] = {
        "kind": "superlevel_H0_periodic",
        "working_t": 0.25,
        "beta0": periodic_components(sim.V >= 0.25),
        "barcode": barcode_signature(bars, k=8),
        "betti_profile": betti_profile(sim.V),
    }
    return packet

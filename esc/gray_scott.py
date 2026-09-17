"""CPU Gray–Scott kernel aligned with Morphogen v2 (IAmM3ta/MORPHOS).

Periodic 5-point Laplacian, dt=1, clamp to [0, 1], same preset table.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class RdParams:
    feed: float
    kill: float
    du: float
    dv: float
    name: str = ""


PRESETS: Dict[str, RdParams] = {
    "mitosis": RdParams(0.0367, 0.0649, 0.16, 0.08, "Mitosis"),
    "solitons": RdParams(0.0353, 0.0653, 0.16, 0.08, "Solitons"),
    "pulsing": RdParams(0.025, 0.06, 0.14, 0.07, "Pulsing"),
    "holes": RdParams(0.039, 0.058, 0.16, 0.08, "Holes"),
    "mazes": RdParams(0.029, 0.057, 0.16, 0.08, "Mazes"),
    "fingerprint": RdParams(0.026, 0.061, 0.16, 0.08, "Fingerprint"),
    "spirals": RdParams(0.018, 0.051, 0.16, 0.08, "Spirals"),
    "worms": RdParams(0.046, 0.063, 0.16, 0.08, "Worms"),
    "coral": RdParams(0.0545, 0.062, 0.16, 0.08, "Coral"),
    "uskate": RdParams(0.062, 0.0609, 0.16, 0.08, "Skate"),
}

FEED_BOUNDS = (0.0, 0.1)
KILL_BOUNDS = (0.0, 0.1)
DIFF_BOUNDS = (0.0, 0.5)


@dataclass
class FieldStats:
    mean_u: float
    mean_v: float
    energy: float
    cx: float
    cy: float
    edge: float


class GrayScott:
    def __init__(self, width: int, height: int, params: RdParams):
        self.width = width
        self.height = height
        self.params = params
        n = width * height
        self.u = np.ones(n, dtype=np.float64)
        self.v = np.zeros(n, dtype=np.float64)

    @property
    def U(self) -> NDArray[np.float64]:
        return self.u.reshape(self.height, self.width)

    @property
    def V(self) -> NDArray[np.float64]:
        return self.v.reshape(self.height, self.width)

    def reset(
        self,
        seed_xy: tuple[float, float] = (0.5, 0.5),
        radius: float = 0.04,
        noise: float = 0.02,
        rng_seed: int = 42,
    ) -> None:
        """Canonical shared prior: central square + deterministic speckles."""
        self.u.fill(1.0)
        self.v.fill(0.0)
        U, V = self.U, self.V
        h, w = V.shape
        half = max(6, int(0.08 * min(h, w)))
        cy = int(seed_xy[1] * h)
        cx = int(seed_xy[0] * w)
        U[max(0, cy - half): cy + half, max(0, cx - half): cx + half] = 0.50
        V[max(0, cy - half): cy + half, max(0, cx - half): cx + half] = 1.00
        if noise > 0:
            rng = np.random.default_rng(rng_seed)
            speck = rng.random(V.shape) < noise
            V[speck] = np.maximum(V[speck], 0.80)
            U[speck] = np.minimum(U[speck], 0.50)

    def seed(self, nx: float, ny: float, radius: float = 0.05, strength: float = 1.0) -> None:
        w, h = self.width, self.height
        cx, cy = nx * w, ny * h
        r = radius * min(w, h)
        r2 = r * r
        ys, xs = np.indices((h, w))
        d2 = (xs - cx) ** 2 + (ys - cy) ** 2
        mask = d2 <= r2
        t = 1.0 - d2[mask] / r2
        s = strength * t * t
        V = self.V
        U = self.U
        V[mask] = np.minimum(1.0, V[mask] + s)
        U[mask] = np.maximum(0.0, U[mask] - s * 0.5)

    def step(self, dt: float = 1.0) -> None:
        f, k, du, dv = self.params.feed, self.params.kill, self.params.du, self.params.dv
        U = self.U
        V = self.V
        lap_u = (
            np.roll(U, 1, 0) + np.roll(U, -1, 0) + np.roll(U, 1, 1) + np.roll(U, -1, 1) - 4.0 * U
        )
        lap_v = (
            np.roll(V, 1, 0) + np.roll(V, -1, 0) + np.roll(V, 1, 1) + np.roll(V, -1, 1) - 4.0 * V
        )
        uvv = U * V * V
        U_n = U + (du * lap_u - uvv + f * (1.0 - U)) * dt
        V_n = V + (dv * lap_v + uvv - (k + f) * V) * dt
        np.clip(U_n, 0.0, 1.0, out=self.U)
        np.clip(V_n, 0.0, 1.0, out=self.V)

    def step_n(self, n: int, dt: float = 1.0) -> None:
        for _ in range(n):
            self.step(dt)

    def stats(self) -> FieldStats:
        U, V = self.U, self.V
        h, w = V.shape
        n = h * w
        energy = float(np.mean(V * V))
        mass = float(np.sum(V))
        ys, xs = np.indices((h, w))
        if mass > 1e-6:
            cx = float(np.sum(xs * V) / mass / w)
            cy = float(np.sum(ys * V) / mass / h)
        else:
            cx = cy = 0.5
        gx = np.zeros_like(V)
        gy = np.zeros_like(V)
        gx[:, 1:-1] = V[:, 2:] - V[:, :-2]
        gy[1:-1, :] = V[2:, :] - V[:-2, :]
        edge = float(np.mean(np.hypot(gx[1:-1, 1:-1], gy[1:-1, 1:-1])))
        return FieldStats(
            mean_u=float(U.mean()),
            mean_v=float(V.mean()),
            energy=energy,
            cx=cx,
            cy=cy,
            edge=edge,
        )

    def grid16(self) -> NDArray[np.float64]:
        V = self.V
        h, w = V.shape
        out = np.zeros((16, 16), dtype=np.float64)
        for gy in range(16):
            y0, y1 = int(gy * h / 16), int((gy + 1) * h / 16)
            for gx in range(16):
                x0, x1 = int(gx * w / 16), int((gx + 1) * w / 16)
                block = V[y0:y1, x0:x1]
                out[gy, gx] = float(block.mean()) if block.size else 0.0
        return out

    def moments_grid(self, side: int = 8) -> NDArray[np.float64]:
        V = self.V
        h, w = V.shape
        out = np.zeros((side, side), dtype=np.float64)
        for gy in range(side):
            y0, y1 = int(gy * h / side), int((gy + 1) * h / side)
            for gx in range(side):
                x0, x1 = int(gx * w / side), int((gx + 1) * w / side)
                block = V[y0:y1, x0:x1]
                out[gy, gx] = float(block.mean()) if block.size else 0.0
        return out

    def copy_from(self, other: "GrayScott") -> None:
        self.u[:] = other.u
        self.v[:] = other.v
        self.params = other.params

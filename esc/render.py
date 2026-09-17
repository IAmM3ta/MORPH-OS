"""Field to PNG helpers (imperative baseline and visual diffs)."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image


def field_to_rgb(v: np.ndarray, u: np.ndarray | None = None) -> Image.Image:
    v = np.clip(v, 0.0, 1.0)
    h, w = v.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[..., 1] = (40 + 200 * v).astype(np.uint8)
    rgb[..., 2] = (30 + 220 * v).astype(np.uint8)
    if u is not None:
        rgb[..., 0] = (20 + 40 * np.clip(1.0 - u, 0, 1)).astype(np.uint8)
    else:
        rgb[..., 0] = (12 + 20 * v).astype(np.uint8)
    return Image.fromarray(rgb, mode="RGB")


def save_field(path: Path, v: np.ndarray, u: np.ndarray | None = None) -> int:
    img = field_to_rgb(v, u)
    img.save(path, format="PNG", optimize=True)
    return path.stat().st_size


def png_bytes_of(v: np.ndarray, u: np.ndarray | None = None) -> int:
    img = field_to_rgb(v, u)
    buf = BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return len(buf.getvalue())


def side_by_side(src_v: np.ndarray, rec_v: np.ndarray, path: Path) -> None:
    a = field_to_rgb(src_v)
    b = field_to_rgb(rec_v)
    delta = np.clip(np.abs(src_v - rec_v) * 3.0, 0, 1)
    c = field_to_rgb(delta)
    w, h = a.size
    canvas = Image.new("RGB", (w * 3 + 8, h))
    canvas.paste(a, (0, 0))
    canvas.paste(b, (w + 4, 0))
    canvas.paste(c, (2 * w + 8, 0))
    canvas.save(path, format="PNG", optimize=True)

"""Structural metrics: ontology / homology proxies, not pixel PSNR."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from gray_scott import GrayScott
from manifold import SharedManifold, feature_vector, scaled_nearest


def blob_count(sim: GrayScott, thresh: float = 0.25) -> int:
    mask = sim.V >= thresh
    labeled, n = ndimage.label(mask)
    if n == 0:
        return 0
    areas = ndimage.sum(mask, labeled, index=range(1, n + 1))
    if np.isscalar(areas):
        areas = [areas]
    return int(sum(1 for a in areas if a >= 8))


@dataclass
class StructuralReport:
    predicted_class: str
    true_class: str
    class_match: bool
    edge_src: float
    edge_rec: float
    edge_rel_error: float
    mean_v_src: float
    mean_v_rec: float
    mean_v_rel_error: float
    centroid_error: float
    blobs_src: int
    blobs_rec: int
    blob_abs_error: int
    mse_v: float
    packet_bytes: int
    png_bytes: int
    grid16_bytes: int
    ratio_vs_png: float


def compare(src, rec, manifold, true_class, packet_bytes, png_bytes, grid16_bytes):
    x_rec = feature_vector(rec)
    pred, _, _ = scaled_nearest(manifold, x_rec)
    st_s, st_r = src.stats(), rec.stats()
    edge_rel = abs(st_r.edge - st_s.edge) / max(st_s.edge, 1e-8)
    mv_rel = abs(st_r.mean_v - st_s.mean_v) / max(st_s.mean_v, 1e-8)
    c_err = float(np.hypot(st_r.cx - st_s.cx, st_r.cy - st_s.cy))
    b_s, b_r = blob_count(src), blob_count(rec)
    mse = float(np.mean((src.V - rec.V) ** 2))
    return StructuralReport(
        predicted_class=pred,
        true_class=true_class,
        class_match=pred == true_class,
        edge_src=st_s.edge,
        edge_rec=st_r.edge,
        edge_rel_error=edge_rel,
        mean_v_src=st_s.mean_v,
        mean_v_rec=st_r.mean_v,
        mean_v_rel_error=mv_rel,
        centroid_error=c_err,
        blobs_src=b_s,
        blobs_rec=b_r,
        blob_abs_error=abs(b_r - b_s),
        mse_v=mse,
        packet_bytes=packet_bytes,
        png_bytes=png_bytes,
        grid16_bytes=grid16_bytes,
        ratio_vs_png=packet_bytes / max(png_bytes, 1),
    )

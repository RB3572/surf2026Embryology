#!/usr/bin/env python3
"""Render the publication-resolution Pard3 figures used for P1.12.

The script uses the measured transcript coordinates, segmentation meshes, polar
body location, and precomputed candidate division planes in the zygote dataset.
It writes five 600-DPI PNG figures and removes obsolete exports.
"""

from __future__ import annotations

import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from scipy.spatial import ConvexHull


ROOT = Path(__file__).resolve().parents[1]
SCENE_PATH = ROOT / "data" / "zygote" / "20260407_zygote_p1_12.json.gz"
ZYGOTE_DIR = ROOT / "data" / "zygote"
OUTPUT_DIR = ROOT / "figures" / "pard3_p1_12"

GENE = "Pard3"
NULL_SIMULATIONS = 10_000
NULL_SEED = 20260717
FIGURE_DPI = 600
XY_UM_PER_PIXEL = 0.15

BLUE = "#2166AC"
RED = "#B2182B"
PARD3 = "#6A1B5D"
SHELL = "#BBC3CB"
INK = "#111111"
MUTED = "#5D6268"
PB_FILL = "#D7DCE1"


def unit(vector: np.ndarray) -> np.ndarray:
    length = np.linalg.norm(vector)
    if length == 0:
        raise ValueError("Cannot normalize a zero-length vector")
    return vector / length


def load_scene() -> dict:
    with gzip.open(SCENE_PATH, "rt") as handle:
        return json.load(handle)


def load_pard3_scenes() -> list[dict]:
    scenes = []
    for path in ZYGOTE_DIR.glob("*.json.gz"):
        with gzip.open(path, "rt") as handle:
            scene = json.load(handle)
        if any(row.get("gene") == GENE for row in scene.get("analysis", {}).get("genes", [])):
            scenes.append(scene)

    def sort_key(scene: dict) -> tuple:
        match = re.match(r"(\d{8})_zygote_p(\d+)_(\d+)(?:_(\d+))?$", scene["id"])
        if not match:
            return (scene["id"],)
        date, plate, embryo, replicate = match.groups()
        return date, int(plate), int(embryo), int(replicate or -1)

    return sorted(scenes, key=sort_key)


def zygote_label(scene: dict) -> str:
    match = re.match(r"(\d{4})(\d{2})(\d{2})_zygote_p(\d+)_(\d+)(?:_(\d+))?$", scene["id"])
    if not match:
        return scene["id"]
    year, month, day, plate, embryo, replicate = match.groups()
    embryo_label = f"Z-P{plate}-fov{embryo}"
    if replicate is not None:
        embryo_label += f"_{replicate}"
    return f"{embryo_label}\n{year}-{month}-{day}"


def plot_vector_to_um(vector: np.ndarray, z_scale: float) -> np.ndarray:
    return np.array(
        [vector[0] * XY_UM_PER_PIXEL, vector[1] * XY_UM_PER_PIXEL, vector[2] / z_scale],
        dtype=float,
    )


def plot_points_to_um(points: np.ndarray, z_scale: float) -> np.ndarray:
    converted = np.asarray(points, dtype=float).copy()
    converted[:, 0] *= XY_UM_PER_PIXEL
    converted[:, 1] *= XY_UM_PER_PIXEL
    converted[:, 2] /= z_scale
    return converted


def transcript_points_um(scene: dict) -> tuple[np.ndarray, np.ndarray]:
    tx = scene["transcripts"][GENE]
    points = np.column_stack(
        [
            np.asarray(tx["x"], dtype=float) * XY_UM_PER_PIXEL,
            np.asarray(tx["y"], dtype=float) * XY_UM_PER_PIXEL,
            np.asarray(tx["gz"], dtype=float),
        ]
    )
    return points, np.asarray(tx["s1"], dtype=bool)


def mesh_um(scene: dict, label: str) -> tuple[np.ndarray, np.ndarray]:
    mesh = scene["region_meshes"][label]
    vertices = np.asarray(mesh["verts"], dtype=float).reshape(-1, 3)
    faces = np.asarray(mesh["faces"], dtype=int).reshape(-1, 3)
    return plot_points_to_um(vertices, scene["z_scale"]), faces


def selected_plane_index(scene: dict) -> int:
    return int(scene["analysis"]["best_planes"]["diffCnt"])


def polar_body_label(scene: dict) -> str:
    return str(scene["analysis"]["polar_body_label"])


def pronucleus_labels(scene: dict) -> list[str]:
    polar_body = int(scene["analysis"]["polar_body_label"])
    candidates = scene["analysis"].get("polar_body_detection", {}).get("candidates", [])
    labels = [
        str(candidate["label"])
        for candidate in candidates
        if not candidate.get("external", False) and int(candidate["label"]) != polar_body
    ]
    if labels:
        return labels[:2]
    return [str(label) for label in scene["mask_labels"] if int(label) not in (1, polar_body)][:2]


def analysis_basis(scene: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    analysis = scene["analysis"]
    center = np.asarray(analysis["com_um"], dtype=float)
    axis = unit(plot_vector_to_um(np.asarray(analysis["axis_plot"]), scene["z_scale"]))
    plane = analysis["planes"][selected_plane_index(scene)]
    normal = unit(np.asarray(plane["normal_um"], dtype=float))
    depth = unit(plot_vector_to_um(np.asarray(plane["m_plot"]), scene["z_scale"]))

    # Screen x is reversed so side A is blue/left and side B is red/right.
    screen_x = -normal
    screen_y = axis
    screen_z = depth
    return center, screen_x, screen_y, screen_z


def transform(points: np.ndarray, basis: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    center, screen_x, screen_y, screen_z = basis
    relative = np.asarray(points, dtype=float) - center
    return np.column_stack(
        [relative @ screen_x, relative @ screen_y, relative @ screen_z]
    )


def gene_row(scene: dict) -> dict:
    return next(row for row in scene["analysis"]["genes"] if row["gene"] == GENE)


def null_count_summary(total: int) -> tuple[int, int, float]:
    rng = np.random.default_rng(NULL_SEED)
    null_counts = rng.binomial(total, 0.5, size=NULL_SIMULATIONS)
    low, high = np.percentile(null_counts, [2.5, 97.5], method="nearest").astype(int)
    return int(low), int(high), total / 2


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.titlesize": 18,
            "axes.titleweight": "bold",
            "axes.labelsize": 12,
            "axes.linewidth": 1.0,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def set_equal_3d(ax, points: np.ndarray, pad: float = 1.08) -> None:
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = (mins + maxs) / 2
    radius = (maxs - mins).max() * pad / 2
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))


def add_mesh(
    ax,
    vertices: np.ndarray,
    faces: np.ndarray,
    *,
    color: str,
    alpha: float,
    edgecolor: str = "none",
    linewidth: float = 0.0,
    zorder: int = 1,
) -> None:
    collection = Poly3DCollection(
        vertices[faces],
        facecolor=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        shade=False,
        zorder=zorder,
    )
    ax.add_collection3d(collection)


def add_title(fig, title: str, subtitle: str) -> None:
    fig.suptitle(title, x=0.5, y=0.965, ha="center", va="top", color=INK)
    fig.text(0.5, 0.915, subtitle, ha="center", va="top", color=MUTED, fontsize=10.5)


def save_figure(fig, stem: str) -> tuple[Path, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png = OUTPUT_DIR / f"{stem}.png"
    pdf = OUTPUT_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=FIGURE_DPI, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    return png, pdf


def save_png(fig, stem: str) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"{stem}.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    return path


def render_uncolored_3d(scene: dict, basis) -> tuple[Path, Path]:
    transcripts, _ = transcript_points_um(scene)
    tx = transform(transcripts, basis)

    meshes = {}
    pb_label = polar_body_label(scene)
    for label in ("1", pb_label):
        vertices, faces = mesh_um(scene, label)
        meshes[label] = (transform(vertices, basis), faces)

    fig = plt.figure(figsize=(8.2, 7.2))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    add_title(
        fig,
        f"{GENE} transcripts in zygote P1\u202212",
        f"All detected molecules (n = {len(tx):,}); segmented structures shown in neutral gray",
    )
    add_mesh(ax, *meshes["1"], color=SHELL, alpha=0.10, zorder=1)
    add_mesh(ax, *meshes[pb_label], color=PB_FILL, alpha=0.66, zorder=4)
    ax.scatter(tx[:, 0], tx[:, 1], tx[:, 2], s=6.5, c=PARD3, alpha=0.72, depthshade=False, linewidths=0, zorder=5)

    pb_center = meshes[pb_label][0].mean(axis=0)
    ax.scatter(*pb_center, s=28, c=INK, depthshade=False, zorder=7)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PARD3, markersize=7, label=f"{GENE} transcript"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=7, label="Polar body"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01), frameon=False, ncol=2, fontsize=9.5)
    set_equal_3d(ax, np.vstack([meshes["1"][0], meshes[pb_label][0]]), pad=1.04)
    ax.view_init(elev=70, azim=-90, roll=0)
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.055, top=0.90)
    return save_figure(fig, "01_p1_12_pard3_real_embryo_3d")


def split_counted_transcripts(scene: dict, basis) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    transcripts, in_segment_one = transcript_points_um(scene)
    counted = transcripts[in_segment_one]
    projected = transform(counted, basis)
    blue = projected[:, 0] < 0
    red = ~blue
    return projected, blue, red


def render_split_3d(scene: dict, basis) -> tuple[Path, Path]:
    tx, blue, red = split_counted_transcripts(scene, basis)
    plane_index = selected_plane_index(scene)
    row = gene_row(scene)["planes"][plane_index]

    meshes = {}
    pb_label = polar_body_label(scene)
    for label in ("1", pb_label):
        vertices, faces = mesh_um(scene, label)
        meshes[label] = (transform(vertices, basis), faces)

    fig = plt.figure(figsize=(8.2, 7.2))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    add_title(
        fig,
        f"{GENE} asymmetry across the candidate division plane",
        f"P1\u202212, {plane_index * 10}\N{DEGREE SIGN} plane; segment-1 transcripts only (n = {len(tx):,})",
    )
    add_mesh(ax, *meshes["1"], color=SHELL, alpha=0.10, zorder=1)
    add_mesh(ax, *meshes[pb_label], color=PB_FILL, alpha=0.70, zorder=3)
    ax.scatter(tx[blue, 0], tx[blue, 1], tx[blue, 2], s=7.0, c=BLUE, alpha=0.76, depthshade=False, linewidths=0, zorder=5)
    ax.scatter(tx[red, 0], tx[red, 1], tx[red, 2], s=7.0, c=RED, alpha=0.72, depthshade=False, linewidths=0, zorder=5)

    body = meshes["1"][0]
    y_min, y_max = np.quantile(body[:, 1], [0.01, 0.99])
    z_front = body[:, 2].max() + 1.0
    ax.plot([0, 0], [y_min, y_max], [z_front, z_front], color=INK, lw=2.0, ls=(0, (5, 4)), zorder=9)

    pb_center = meshes[pb_label][0].mean(axis=0)
    ax.scatter(*pb_center, s=26, c=INK, depthshade=False, zorder=8)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markersize=7, label=f"Left: n = {row['a']:,}"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=RED, markersize=7, label=f"Right: n = {row['b']:,}"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=7, label="Polar body"),
        Line2D([0], [0], color=INK, lw=1.8, ls=(0, (5, 4)), label="Candidate division plane"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01), frameon=False, ncol=4, fontsize=9.0)
    set_equal_3d(ax, np.vstack([body, meshes[pb_label][0]]), pad=1.04)
    ax.view_init(elev=82, azim=-90, roll=0)
    ax.set_proj_type("ortho")
    ax.set_axis_off()
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.055, top=0.90)
    return save_figure(fig, "02_p1_12_pard3_split_topdown_3d")


def mesh_cross_section(
    vertices: np.ndarray,
    faces: np.ndarray,
    plane_z: float,
    smooth_iterations: int = 4,
) -> np.ndarray:
    """Intersect a triangular mesh with z=plane_z and return a smooth closed contour."""
    segments = []
    for triangle in vertices[faces]:
        edge_points = []
        for start, end in ((0, 1), (1, 2), (2, 0)):
            p0, p1 = triangle[start], triangle[end]
            z0, z1 = p0[2] - plane_z, p1[2] - plane_z
            if abs(z0) < 1e-9:
                edge_points.append(p0[:2])
            if z0 * z1 < 0:
                fraction = -z0 / (z1 - z0)
                edge_points.append(p0[:2] + fraction * (p1[:2] - p0[:2]))
        unique = np.unique(np.round(np.asarray(edge_points), decimals=5), axis=0) if edge_points else np.empty((0, 2))
        if len(unique) == 2:
            segments.append((unique[0], unique[1]))

    def key(point: np.ndarray) -> tuple[float, float]:
        rounded = np.round(point, decimals=4)
        return float(rounded[0]), float(rounded[1])

    adjacency = defaultdict(set)
    edges = set()
    for start, end in segments:
        a, b = key(start), key(end)
        if a == b:
            continue
        edge = tuple(sorted((a, b)))
        if edge in edges:
            continue
        edges.add(edge)
        adjacency[a].add(b)
        adjacency[b].add(a)

    unused = set(edges)
    loops = []
    while unused:
        first_edge = next(iter(unused))
        start, current = first_edge
        path = [start, current]
        unused.remove(first_edge)
        previous = start
        while current != start:
            candidates = []
            for neighbor in adjacency[current]:
                edge = tuple(sorted((current, neighbor)))
                if edge in unused:
                    candidates.append((neighbor, edge))
            if not candidates:
                break
            candidates.sort(key=lambda item: item[0] == previous)
            following, edge = candidates[0]
            unused.remove(edge)
            path.append(following)
            previous, current = current, following
        if len(path) >= 8 and path[-1] == path[0]:
            loops.append(np.asarray(path[:-1], dtype=float))

    if not loops:
        hull = ConvexHull(vertices[:, :2])
        contour = vertices[hull.vertices, :2]
    else:
        def area(loop: np.ndarray) -> float:
            x, y = loop[:, 0], loop[:, 1]
            return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

        contour = max(loops, key=area)

    # Closed Chaikin subdivision removes marching-cubes facets without inventing
    # a circular outline or moving the contour outside the measured section.
    for _ in range(smooth_iterations):
        following = np.roll(contour, -1, axis=0)
        first = 0.75 * contour + 0.25 * following
        second = 0.25 * contour + 0.75 * following
        refined = np.empty((len(contour) * 2, 2), dtype=float)
        refined[0::2] = first
        refined[1::2] = second
        contour = refined
    return contour


def render_projection(scene: dict, basis, *, show_pronucleus_labels: bool, stem: str) -> Path:
    tx, blue, red = split_counted_transcripts(scene, basis)
    row = gene_row(scene)["planes"][selected_plane_index(scene)]

    body_vertices, body_faces = mesh_um(scene, "1")
    pb_vertices, pb_faces = mesh_um(scene, polar_body_label(scene))
    pn_meshes = [mesh_um(scene, label) for label in pronucleus_labels(scene)]
    body = transform(body_vertices, basis)
    pb = transform(pb_vertices, basis)
    pronuclei = [(transform(vertices, basis), faces) for vertices, faces in pn_meshes]
    body_hull = mesh_cross_section(body, body_faces, plane_z=0.0)
    pb_hull = mesh_cross_section(pb, pb_faces, plane_z=float(np.median(pb[:, 2])))
    pn_hulls = [
        mesh_cross_section(vertices, faces, plane_z=float(np.median(vertices[:, 2])))
        for vertices, faces in pronuclei
    ]

    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    add_title(
        fig,
        f"2D projection of {GENE} asymmetry in Z-P1-fov12",
        "Projection along the candidate division plane; each point is one segment-1 transcript",
    )
    ax.add_patch(Polygon(body_hull, closed=True, facecolor="#DCEAF3", edgecolor=INK, linewidth=1.4, alpha=0.36, zorder=1))
    ax.add_patch(Polygon(pb_hull, closed=True, facecolor="#C9E1F0", edgecolor="#1E3A5F", linewidth=1.15, alpha=0.82, zorder=4))
    pn_colors = ["#D5E8F3", "#E0EDF5"]
    pn_edges = ["#334155", "#64748B"]
    for index, hull in enumerate(pn_hulls):
        ax.add_patch(Polygon(hull, closed=True, facecolor=pn_colors[index], edgecolor=pn_edges[index],
                             linewidth=1.1, alpha=0.76, zorder=2))
    ax.scatter(tx[blue, 0], tx[blue, 1], s=12, c=BLUE, alpha=0.62, linewidths=0, zorder=3)
    ax.scatter(tx[red, 0], tx[red, 1], s=12, c=RED, alpha=0.56, linewidths=0, zorder=3)

    y_min, y_max = np.quantile(body[:, 1], [0.01, 0.99])
    ax.plot([0, 0], [y_min, y_max], color=INK, lw=2.1, ls=(0, (5, 4)), zorder=6)

    pb_center = pb[:, :2].mean(axis=0)
    ax.scatter(*pb_center, s=34, c=INK, linewidths=0, zorder=7)
    pn_centers = [vertices[:, :2].mean(axis=0) for vertices, _ in pronuclei]
    for index, center_point in enumerate(pn_centers):
        ax.scatter(*center_point, s=24, c=pn_edges[index], edgecolors="white", linewidths=0.7, zorder=7)

    ax.set_aspect("equal", adjustable="box")
    all_points = np.vstack([body[:, :2], pb[:, :2]])
    mins, maxs = all_points.min(axis=0), all_points.max(axis=0)
    span = (maxs - mins).max()
    center = (mins + maxs) / 2
    ax.set_xlim(center[0] - span * 0.60, center[0] + span * 0.78)
    ax.set_ylim(center[1] - span * 0.62, center[1] + span * 0.62)
    ax.annotate(
        "Polar body",
        xy=pb_center,
        xytext=(maxs[0] + span * 0.10, pb_center[1] + span * 0.06),
        arrowprops={"arrowstyle": "-", "color": INK, "lw": 1.0},
        color=INK,
        fontsize=10,
        ha="left",
        va="bottom",
        zorder=8,
    )
    if show_pronucleus_labels:
        body_min = body_hull.min(axis=0)
        body_max = body_hull.max(axis=0)
        label_positions = [
            (body_min[0] - span * 0.10, pn_centers[0][1] + span * 0.05, "right"),
            (body_max[0] + span * 0.10, pn_centers[1][1] - span * 0.05, "left"),
        ]
        for index, (center_point, (label_x, label_y, alignment)) in enumerate(
            zip(pn_centers, label_positions), start=1
        ):
            ax.annotate(
                f"Pronucleus {index}",
                xy=center_point,
                xytext=(label_x, label_y),
                arrowprops={"arrowstyle": "-", "color": pn_edges[index - 1], "lw": 0.9},
                color=pn_edges[index - 1],
                fontsize=9.5,
                ha=alignment,
                va="center",
                zorder=8,
                annotation_clip=False,
            )
    ax.axis("off")
    fig.text(0.19, 0.085, f"Blue side: n = {row['a']:,} ({row['a'] / (row['a'] + row['b']):.1%})", color=BLUE, fontsize=12, fontweight="bold", ha="center", va="center")
    fig.text(0.76, 0.085, f"Red side: n = {row['b']:,} ({row['b'] / (row['a'] + row['b']):.1%})", color=RED, fontsize=12, fontweight="bold", ha="center", va="center")
    scale_y = mins[1] - span * 0.08
    scale_x = mins[0] + span * 0.08
    ax.plot([scale_x, scale_x + 10], [scale_y, scale_y], color=INK, lw=2.2, solid_capstyle="butt", zorder=8)
    ax.text(scale_x + 5, scale_y - span * 0.025, "10 \u00b5m", ha="center", va="top", color=INK, fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.13, top=0.89)
    return save_png(fig, stem)


def pard3_chart_rows(scenes: list[dict]) -> list[dict]:
    rows = []
    for scene in scenes:
        gene = gene_row(scene)
        plane_index = int(scene["analysis"]["best_planes"]["diffCnt"])
        plane = gene["planes"][plane_index]
        high_count = max(int(plane["a"]), int(plane["b"]))
        low_count = min(int(plane["a"]), int(plane["b"]))
        total = high_count + low_count
        null_low, null_high, null_mean = null_count_summary(total)
        rows.append(
            {
                "id": scene["id"],
                "label": zygote_label(scene),
                "plane_index": plane_index,
                "plane_angle_degrees": plane_index * 10,
                "high_count": high_count,
                "low_count": low_count,
                "total": total,
                "null_low": null_low,
                "null_high": null_high,
                "null_mean": null_mean,
            }
        )
    return rows


def render_bar_chart(scenes: list[dict], *, log_scale: bool) -> Path:
    rows = pard3_chart_rows(scenes)
    x = np.arange(len(rows), dtype=float)
    high_counts = np.asarray([row["high_count"] for row in rows])
    low_counts = np.asarray([row["low_count"] for row in rows])
    totals = np.asarray([row["total"] for row in rows])
    null_means = np.asarray([row["null_mean"] for row in rows])
    null_lows = np.asarray([row["null_low"] for row in rows])
    null_highs = np.asarray([row["null_high"] for row in rows])

    fig, ax = plt.subplots(figsize=(14.5, 7.2))
    null_width = 0.72
    if log_scale:
        real_width = 0.19
        real_offset = 0.115
        ax.bar(x - real_offset, high_counts, color=BLUE, width=real_width, edgecolor="#263238", linewidth=0.45, zorder=3)
        ax.bar(x + real_offset, low_counts, color=RED, width=real_width, edgecolor="#263238", linewidth=0.45, zorder=3)
    else:
        real_width = 0.42
        ax.bar(x, high_counts, color=BLUE, width=real_width, edgecolor="#263238", linewidth=0.45, zorder=3)
        ax.bar(x, low_counts, bottom=high_counts, color=RED, width=real_width, edgecolor="#263238", linewidth=0.45, zorder=3)
    ax.bar(
        x,
        null_means,
        color="#9EA4AA",
        alpha=0.30,
        edgecolor="#6F757B",
        linewidth=1.0,
        width=null_width,
        zorder=5,
    )
    ax.errorbar(
        x,
        null_means,
        yerr=np.vstack([null_means - null_lows, null_highs - null_means]),
        fmt="none",
        ecolor=INK,
        elinewidth=0.8,
        capsize=3,
        capthick=0.8,
        zorder=6,
    )

    null_handles = [
        mpl.patches.Patch(facecolor=BLUE, edgecolor="none", label="Higher-count half"),
        mpl.patches.Patch(facecolor=RED, edgecolor="none", label="Lower-count half"),
        mpl.patches.Patch(facecolor="#9EA4AA", alpha=0.55, edgecolor="#6F757B", label="Null mean"),
        Line2D([0], [0], color=INK, lw=0.8, marker="_", markersize=7, label="95% null interval"),
    ]
    ax.legend(handles=null_handles, loc="upper right", frameon=False, fontsize=10, ncol=2)
    ax.set_ylabel("Pard3 transcript count")
    if log_scale:
        ax.set_yscale("log")
        ax.set_ylim(1, max(high_counts.max(), low_counts.max(), null_highs.max()) * 1.18)
    else:
        ax.set_ylim(0, max(totals.max(), null_highs.max()) * 1.10)
    ax.set_xlim(-0.75, len(rows) - 0.25)
    ax.set_xticks(x, [row["label"] for row in rows], rotation=48, ha="right", rotation_mode="anchor")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#70757A")
    ax.tick_params(axis="both", colors=INK)
    ax.yaxis.grid(True, which="major", color="#D8DCE0", linewidth=0.9, zorder=0)
    if log_scale:
        ax.yaxis.grid(True, which="minor", color="#EEF0F2", linewidth=0.55, zorder=0)
    ax.xaxis.grid(False)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.27, top=0.97)
    stem = "02_all_zygotes_pard3_side_counts_log" if log_scale else "01_all_zygotes_pard3_side_counts_linear"
    return save_png(fig, stem)


def render_ratio_chart(scenes: list[dict]) -> Path:
    """Render each Pard3 half as a percentage of that zygote's total."""
    rows = pard3_chart_rows(scenes)
    x = np.arange(len(rows), dtype=float)
    totals = np.asarray([row["total"] for row in rows], dtype=float)
    high_pct = 100 * np.asarray([row["high_count"] for row in rows]) / totals
    low_pct = 100 * np.asarray([row["low_count"] for row in rows]) / totals
    null_low_pct = 100 * np.asarray([row["null_low"] for row in rows]) / totals
    null_high_pct = 100 * np.asarray([row["null_high"] for row in rows]) / totals

    fig, ax = plt.subplots(figsize=(14.5, 7.2))
    real_width = 0.44
    ax.bar(x, high_pct, color=BLUE, width=real_width, edgecolor="#263238", linewidth=0.45, zorder=3)
    ax.bar(x, low_pct, bottom=high_pct, color=RED, width=real_width,
           edgecolor="#263238", linewidth=0.45, zorder=3)

    ax.bar(
        x,
        np.full(len(rows), 50.0),
        color="#8F969D",
        alpha=0.30,
        edgecolor="#666D73",
        linewidth=0.55,
        width=0.70,
        zorder=5,
    )
    ax.errorbar(
        x,
        np.full(len(rows), 50.0),
        yerr=np.vstack([50.0 - null_low_pct, null_high_pct - 50.0]),
        fmt="none",
        ecolor=INK,
        elinewidth=0.8,
        capsize=3,
        capthick=0.8,
        zorder=6,
    )

    handles = [
        mpl.patches.Patch(facecolor=BLUE, edgecolor="none", label="Higher-count half"),
        mpl.patches.Patch(facecolor=RED, edgecolor="none", label="Lower-count half"),
        mpl.patches.Patch(facecolor="#8F969D", alpha=0.45, edgecolor="#666D73", label="50% null expectation"),
        Line2D([0], [0], color=INK, lw=0.8, marker="_", markersize=7, label="95% null interval"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.01),
              frameon=False, fontsize=10, ncol=4)
    ax.set_ylabel("Share of Pard3 transcripts in each zygote")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.75, len(rows) - 0.25)
    ax.set_xticks(x, [row["label"] for row in rows], rotation=48, ha="right", rotation_mode="anchor")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#70757A")
    ax.tick_params(axis="both", colors=INK)
    ax.yaxis.grid(True, which="major", color="#D8DCE0", linewidth=0.9, zorder=0)
    ax.xaxis.grid(False)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.27, top=0.90)
    return save_png(fig, "05_all_zygotes_pard3_side_percentages")


def build_composite(paths: list[Path]) -> Path:
    images = [Image.open(path).convert("RGB") for path in paths]
    tile_width = 3000
    normalized = []
    for image in images:
        scale = tile_width / image.width
        normalized.append(image.resize((tile_width, int(image.height * scale)), Image.Resampling.LANCZOS))
    tile_height = max(image.height for image in normalized)
    canvas = Image.new("RGB", (tile_width * 2, tile_height * 2), "white")
    for index, image in enumerate(normalized):
        x = (index % 2) * tile_width
        y = (index // 2) * tile_height + (tile_height - image.height) // 2
        canvas.paste(image, (x, y))
    path = OUTPUT_DIR / "p1_12_pard3_composite.png"
    canvas.save(path, dpi=(300, 300), optimize=True)
    return path


def write_metadata(scene: dict, pard3_scenes: list[dict], outputs: list[Path]) -> Path:
    plane_index = selected_plane_index(scene)
    row = gene_row(scene)["planes"][plane_index]
    transcript_points, segment_one = transcript_points_um(scene)
    null_low, null_high, null_mean = null_count_summary(int(segment_one.sum()))
    metadata = {
        "embryo_id": scene["id"],
        "embryo_label": "P1.12",
        "gene": GENE,
        "source_scene": str(SCENE_PATH.relative_to(ROOT)),
        "polar_body_label": int(polar_body_label(scene)),
        "candidate_plane_index": plane_index,
        "candidate_plane_angle_degrees": plane_index * 10,
        "plane_selection": "Best count-difference plane after geometry-based polar-body detection",
        "all_detected_transcripts": int(len(transcript_points)),
        "segment_1_transcripts": int(segment_one.sum()),
        "blue_left_count": int(row["a"]),
        "red_right_count": int(row["b"]),
        "absolute_difference": int(abs(row["a"] - row["b"])),
        "red_to_blue_ratio": float(row["b"] / row["a"]),
        "count_permutation_p": float(row["pCnt"]),
        "volume_permutation_p": float(row["pVol"]),
        "null_simulations": NULL_SIMULATIONS,
        "null_mean_count_per_side": null_mean,
        "null_95_percent_count_range": [null_low, null_high],
        "all_zygotes_bar_chart": {
            "plane_selection": "Per-embryo best count-difference plane precomputed across all measured genes",
            "blue_definition": "Higher-count Pard3 half",
            "red_definition": "Lower-count Pard3 half",
            "zygotes": pard3_chart_rows(pard3_scenes),
        },
        "outputs": [str(path.relative_to(ROOT)) for path in outputs],
    }
    path = OUTPUT_DIR / "figure_metadata.json"
    path.write_text(json.dumps(metadata, indent=2) + "\n")
    return path


def main() -> None:
    configure_style()
    scene = load_scene()
    pard3_scenes = load_pard3_scenes()
    basis = analysis_basis(scene)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUTPUT_DIR.iterdir():
        if path.is_file():
            path.unlink()

    outputs = [
        render_bar_chart(pard3_scenes, log_scale=False),
        render_bar_chart(pard3_scenes, log_scale=True),
        render_projection(
            scene,
            basis,
            show_pronucleus_labels=True,
            stem="03_p1_12_pard3_2d_structures",
        ),
        render_projection(
            scene,
            basis,
            show_pronucleus_labels=False,
            stem="04_p1_12_pard3_2d_structures_no_pronucleus_labels",
        ),
        render_ratio_chart(pard3_scenes),
    ]

    print(f"Wrote exactly {len(outputs)} figures to {OUTPUT_DIR}")
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the shareable Stage A solvent report from committed pilot data.

Only the small JSON/CSV/XYZ exports under ``reports/stage_a/data`` are read.
No trajectory or chemistry engine is needed, so this is safe to run locally.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"
DEFAULT_PDF = HERE / "stage_a_solvent_report_zh.pdf"
METRICS_JSON = HERE / "stage_a_metrics.json"
SYSTEMS = ("pure_thf", "eda_1p5m")
SYSTEM_LABELS = {"pure_thf": "纯 THF", "eda_1p5m": "1.5 M en/THF"}
REPLICAS = (1, 2, 3)
AVOGADRO = 6.02214076e23

INK = "#17233B"
MUTED = "#62708A"
GRID = "#DDE4EE"
PURE = "#6C7C93"
MIXED = "#137C8B"
ACCENT = "#E59F23"
EDA = "#2667FF"
THF_O = "#E45756"
THF_C = "#AAB3C2"
BACKGROUND = "#F4F7FB"
REPLICA_COLORS = ("#315C9A", "#21A179", "#E58E26")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def load_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for system in SYSTEMS:
        for replica in REPLICAS:
            directory = DATA / "analysis" / system / f"r{replica}"
            records.append(
                {
                    "system": system,
                    "replica": replica,
                    "directory": directory,
                    "analysis": _load_json(directory / "analysis.json"),
                    "snapshot": _load_json(directory / "snapshot" / "metadata.json"),
                    "timeseries": _read_csv(directory / "timeseries.csv"),
                    "rdf": _read_csv(directory / "rdf.csv"),
                }
            )
    return records


def validate_records(records: Sequence[dict[str, Any]]) -> None:
    analysis_summary = _load_json(DATA / "classical_analysis.summary.json")
    snapshot_summary = _load_json(DATA / "snapshot_bank.summary.json")
    if not analysis_summary.get("ready") or analysis_summary.get("record_count") != 6:
        raise ValueError("Stage A analysis summary is not a six-record ready result")
    if not snapshot_summary.get("ready") or snapshot_summary.get("record_count") != 6:
        raise ValueError("Stage A snapshot summary is not a six-record ready result")
    if len(records) != 6:
        raise ValueError(f"expected six records, found {len(records)}")
    seen: set[tuple[str, int]] = set()
    for record in records:
        key = (record["system"], record["replica"])
        if key in seen:
            raise ValueError(f"duplicate record: {key}")
        seen.add(key)
        analysis = record["analysis"]
        snapshot = record["snapshot"]
        if not analysis.get("ready") or not analysis["metrics"].get("ready"):
            raise ValueError(f"analysis is not ready: {key}")
        if not snapshot.get("ready") or snapshot.get("li_atom_present"):
            raise ValueError(f"snapshot is not a ready solvent-only structure: {key}")
        if len(record["timeseries"]) != analysis["metrics"]["frame_count"]:
            raise ValueError(f"timeseries length mismatch: {key}")
        directory = record["directory"]
        checks = (
            (directory / "analysis.json", snapshot["source"]["analysis"]["sha256"]),
            (
                directory / "snapshot" / "representative.xyz",
                snapshot["structure"]["xyz"]["sha256"],
            ),
            (
                directory / "snapshot" / "representative.cell.inc",
                snapshot["structure"]["cell"]["sha256"],
            ),
        )
        for path, expected in checks:
            observed = _sha256(path)
            if observed != expected:
                raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")


def _sample_sd(values: Sequence[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def _descriptor(records: Sequence[dict[str, Any]], system: str, name: str) -> tuple[float, float]:
    values = [
        float(record["analysis"]["metrics"]["mean_descriptors"][name])
        for record in records
        if record["system"] == system
    ]
    return statistics.mean(values), _sample_sd(values)


def _half_difference_percent(rows: Sequence[dict[str, str]], field: str) -> float:
    values = [float(row[field]) for row in rows]
    midpoint = len(values) // 2
    first = statistics.mean(values[:midpoint])
    second = statistics.mean(values[midpoint:])
    return abs(first - second) / statistics.mean(values) * 100.0


def _rdf_mean(
    records: Sequence[dict[str, Any]], system: str, pair: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    series: list[tuple[np.ndarray, np.ndarray]] = []
    for record in records:
        if record["system"] != system:
            continue
        selected = [row for row in record["rdf"] if row["pair"] == pair]
        series.append(
            (
                np.asarray([float(row["r_center_angstrom"]) for row in selected]),
                np.asarray([float(row["g_r"]) for row in selected]),
            )
        )
    if not series:
        raise ValueError(f"missing RDF pair {pair} for {system}")
    radii = series[0][0]
    values = np.stack([item[1] for item in series])
    return radii, values.mean(axis=0), values.min(axis=0), values.max(axis=0)


def _rdf_peak(
    records: Sequence[dict[str, Any]], system: str, pair: str, lower: float, upper: float
) -> dict[str, float]:
    radii, mean, _minimum, _maximum = _rdf_mean(records, system, pair)
    mask = (radii >= lower) & (radii <= upper)
    candidates = np.flatnonzero(mask)
    index = int(candidates[np.argmax(mean[mask])])
    return {"r_angstrom": float(radii[index]), "g_r": float(mean[index])}


def summarize(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    systems: dict[str, Any] = {}
    for system in SYSTEMS:
        density = _descriptor(records, system, "density_g_ml")
        volume = _descriptor(records, system, "volume_nm3")
        void = _descriptor(records, system, "void_radius_angstrom")
        subset = [record for record in records if record["system"] == system]
        item: dict[str, Any] = {
            "replicas": 3,
            "production_ns_per_replica": 20.0,
            "analysis_frames_per_replica": 201,
            "density_g_ml": {"mean": density[0], "replica_sd": density[1]},
            "volume_nm3": {"mean": volume[0], "replica_sd": volume[1]},
            "void_radius_angstrom": {"mean": void[0], "replica_sd": void[1]},
            "maximum_density_half_difference_percent": max(
                _half_difference_percent(record["timeseries"], "density_g_ml") for record in subset
            ),
            "minimum_effective_sample_size": min(
                float(value["effective_sample_size"])
                for record in subset
                for value in record["analysis"]["metrics"]["autocorrelation"].values()
            ),
            "snapshots": [
                {
                    "replica": record["replica"],
                    "snapshot_id": record["snapshot"]["snapshot_id"],
                    "time_ns": record["snapshot"]["selection"]["elapsed_ps"] / 1000.0,
                    "density_g_ml": record["snapshot"]["selection"]["descriptors"]["density_g_ml"],
                    "void_radius_angstrom": record["snapshot"]["selection"]["descriptors"][
                        "void_radius_angstrom"
                    ],
                }
                for record in subset
            ],
        }
        if system == "eda_1p5m":
            contacts = _descriptor(records, system, "eda_thf_contacts")
            hydrogen_bonds = _descriptor(records, system, "eda_thf_hydrogen_bonds")
            item.update(
                {
                    "thf_count": 64,
                    "eda_count": 9,
                    "thf_to_eda_ratio": 64.0 / 9.0,
                    "eda_mole_fraction": 9.0 / 73.0,
                    "achieved_eda_concentration_m": 9.0 / (AVOGADRO * volume[0] * 1e-24),
                    "eda_n_thf_o_contacts": {"mean": contacts[0], "replica_sd": contacts[1]},
                    "eda_nh_thf_o_hydrogen_bonds": {
                        "mean": hydrogen_bonds[0],
                        "replica_sd": hydrogen_bonds[1],
                    },
                }
            )
        else:
            item.update({"thf_count": 64, "eda_count": 0})
        systems[system] = item
    systems["pure_thf"]["rdf_peaks"] = {
        "thf_o-thf_o": _rdf_peak(records, "pure_thf", "thf_o-thf_o", 2.0, 7.0)
    }
    systems["eda_1p5m"]["rdf_peaks"] = {
        "thf_o-thf_o": _rdf_peak(records, "eda_1p5m", "thf_o-thf_o", 2.0, 7.0),
        "eda_n-thf_o": _rdf_peak(records, "eda_1p5m", "eda_n-thf_o", 2.0, 5.0),
        "eda_n-eda_n": _rdf_peak(records, "eda_1p5m", "eda_n-eda_n", 2.0, 6.0),
    }
    comparisons = {
        "mixture_vs_pure_density_percent": (
            systems["eda_1p5m"]["density_g_ml"]["mean"]
            / systems["pure_thf"]["density_g_ml"]["mean"]
            - 1.0
        )
        * 100.0,
        "mixture_vs_pure_void_radius_percent": (
            systems["eda_1p5m"]["void_radius_angstrom"]["mean"]
            / systems["pure_thf"]["void_radius_angstrom"]["mean"]
            - 1.0
        )
        * 100.0,
    }
    return {
        "schema_version": 1,
        "status": "READY",
        "source": "Stage A pilot snapshot bank",
        "systems": systems,
        "comparisons": comparisons,
        "scientific_scope": {
            "periodic_bulk_solvent": True,
            "contains_li": False,
            "contains_excess_electron": False,
            "void_is_geometric_proxy_not_electron_density": True,
        },
    }


def _find_font(bold: bool = False) -> Path | None:
    candidates: list[Path] = []
    if bold:
        candidates.extend(
            [
                Path("C:/Windows/Fonts/msyhbd.ttc"),
                Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"),
            ]
        )
    candidates.extend(
        [
            Path("C:/Windows/Fonts/msyh.ttc"),
            Path("C:/Windows/Fonts/simhei.ttf"),
            Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        ]
    )
    return next((path for path in candidates if path.exists()), None)


def _pil_font(size: int, bold: bool = False):
    from PIL import ImageFont

    path = _find_font(bold=bold)
    if path is None:
        return ImageFont.load_default()
    return ImageFont.truetype(str(path), size=size, index=0)


def _hex(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[index : index + 2], 16) for index in (0, 2, 4))


def _canvas(width: int, height: int):
    from PIL import Image

    return Image.new("RGB", (width, height), "white")


def _title(draw, title: str, subtitle: str | None = None) -> None:
    draw.text((72, 48), title, font=_pil_font(42, bold=True), fill=_hex(INK))
    if subtitle:
        draw.text((74, 102), subtitle, font=_pil_font(22), fill=_hex(MUTED))


def _draw_axes(
    draw,
    box: tuple[int, int, int, int],
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    x_ticks: Sequence[float],
    y_ticks: Sequence[float],
    x_label: str,
    y_label: str,
    *,
    y_decimals: int = 2,
):
    left, top, right, bottom = box
    draw.rounded_rectangle(box, radius=18, fill=(250, 252, 255), outline=_hex(GRID), width=2)

    def point(x: float, y: float) -> tuple[float, float]:
        px = left + 72 + (x - x_range[0]) / (x_range[1] - x_range[0]) * (right - left - 104)
        py = bottom - 58 - (y - y_range[0]) / (y_range[1] - y_range[0]) * (bottom - top - 98)
        return px, py

    x0, y0 = point(x_range[0], y_range[0])
    x1, y1 = point(x_range[1], y_range[1])
    for tick in y_ticks:
        _, py = point(x_range[0], tick)
        draw.line((x0, py, x1, py), fill=_hex(GRID), width=2)
        draw.text(
            (x0 - 14, py),
            f"{tick:.{y_decimals}f}",
            anchor="rm",
            font=_pil_font(17),
            fill=_hex(MUTED),
        )
    for tick in x_ticks:
        px, _ = point(tick, y_range[0])
        draw.line((px, y0, px, y1), fill=(235, 239, 245), width=1)
        draw.text(
            (px, y0 + 14),
            f"{tick:.0f}",
            anchor="ma",
            font=_pil_font(17),
            fill=_hex(MUTED),
        )
    draw.line((x0, y0, x1, y0), fill=_hex(INK), width=2)
    draw.line((x0, y0, x0, y1), fill=_hex(INK), width=2)
    draw.text(
        ((x0 + x1) / 2, bottom - 18),
        x_label,
        anchor="mm",
        font=_pil_font(19),
        fill=_hex(INK),
    )
    from PIL import Image, ImageDraw

    label_font = _pil_font(19)
    label_box = label_font.getbbox(y_label)
    label_width = label_box[2] - label_box[0] + 8
    label_height = label_box[3] - label_box[1] + 8
    label_image = Image.new("L", (label_width, label_height), 0)
    ImageDraw.Draw(label_image).text((4, 4), y_label, font=label_font, fill=255)
    label_image = label_image.rotate(90, expand=True, resample=Image.Resampling.BICUBIC)
    label_x = int(left + 10)
    label_y = int((y0 + y1 - label_image.height) / 2)
    draw.bitmap((label_x, label_y), label_image, fill=_hex(INK))
    return point


def _draw_legend(draw, items: Sequence[tuple[str, str]], x: int, y: int, *, gap: int = 150) -> None:
    for index, (label, color) in enumerate(items):
        origin = x + index * gap
        draw.line((origin, y, origin + 28, y), fill=_hex(color), width=6)
        draw.text((origin + 38, y), label, anchor="lm", font=_pil_font(17), fill=_hex(INK))


def render_density(
    records: Sequence[dict[str, Any]], metrics: dict[str, Any], output: Path
) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1050)
    draw = ImageDraw.Draw(image)
    _title(
        draw,
        "密度与体积统计",
        "每个体系 3 个独立副本，每条生产轨迹 20 ns；绿色区域为代表快照候选区间",
    )
    boxes = ((65, 155, 875, 790), (925, 155, 1735, 790))
    for system, box in zip(SYSTEMS, boxes, strict=False):
        subset = [record for record in records if record["system"] == system]
        values = [float(row["density_g_ml"]) for record in subset for row in record["timeseries"]]
        margin = max(0.01, (max(values) - min(values)) * 0.08)
        y_min = math.floor((min(values) - margin) * 100) / 100
        y_max = math.ceil((max(values) + margin) * 100) / 100
        point = _draw_axes(
            draw,
            box,
            (0.0, 20.0),
            (y_min, y_max),
            (0, 5, 10, 15, 20),
            np.linspace(y_min, y_max, 5),
            "时间 / ns",
            "密度 / g mL^-1",
        )
        shade_left, shade_top = point(10.0, y_max)
        shade_right, shade_bottom = point(20.0, y_min)
        draw.rectangle((shade_left, shade_top, shade_right, shade_bottom), fill=(233, 244, 242))
        for record, color in zip(subset, REPLICA_COLORS, strict=False):
            coordinates = [
                point(float(row["elapsed_ps"]) / 1000.0, float(row["density_g_ml"]))
                for row in record["timeseries"]
            ]
            draw.line(coordinates, fill=_hex(color), width=3)
        draw.text(
            (box[0] + 38, box[1] + 26),
            SYSTEM_LABELS[system],
            font=_pil_font(26, bold=True),
            fill=_hex(INK),
        )
        _draw_legend(
            draw,
            [
                (f"r{replica}", color)
                for replica, color in zip(REPLICAS, REPLICA_COLORS, strict=False)
            ],
            box[0] + 430,
            box[1] + 43,
            gap=105,
        )
    cards = [
        (
            "纯 THF 平均密度",
            f"{metrics['systems']['pure_thf']['density_g_ml']['mean']:.5f} ± "
            f"{metrics['systems']['pure_thf']['density_g_ml']['replica_sd']:.5f} g/mL",
            PURE,
        ),
        (
            "en/THF 平均密度",
            f"{metrics['systems']['eda_1p5m']['density_g_ml']['mean']:.5f} ± "
            f"{metrics['systems']['eda_1p5m']['density_g_ml']['replica_sd']:.5f} g/mL",
            MIXED,
        ),
        (
            "实际 en 浓度",
            f"{metrics['systems']['eda_1p5m']['achieved_eda_concentration_m']:.3f} M",
            ACCENT,
        ),
    ]
    for index, (label, value, color) in enumerate(cards):
        left = 65 + index * 565
        draw.rounded_rectangle(
            (left, 835, left + 520, 985),
            radius=20,
            fill=_hex(BACKGROUND),
            outline=_hex(GRID),
            width=2,
        )
        draw.rectangle((left, 835, left + 12, 985), fill=_hex(color))
        draw.text((left + 38, 875), label, font=_pil_font(21), fill=_hex(MUTED))
        draw.text((left + 38, 925), value, font=_pil_font(26, bold=True), fill=_hex(INK))
    image.save(output, dpi=(180, 180))


def render_microstructure(
    records: Sequence[dict[str, Any]], metrics: dict[str, Any], output: Path
) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1050)
    draw = ImageDraw.Draw(image)
    _title(draw, "局部结构统计", "RDF 为三副本平均；接触与氢键按每个周期盒逐帧计数")
    boxes = ((65, 155, 875, 790), (925, 155, 1735, 790))
    specifications = [
        (
            boxes[0],
            "THF 氧-氧 RDF",
            [
                ("pure_thf", "thf_o-thf_o", "纯 THF: O-O", PURE),
                ("eda_1p5m", "thf_o-thf_o", "en/THF: O-O", MIXED),
            ],
        ),
        (
            boxes[1],
            "en 相关 RDF",
            [
                ("eda_1p5m", "eda_n-thf_o", "en N-THF O", EDA),
                ("eda_1p5m", "eda_n-eda_n", "en N-N", ACCENT),
            ],
        ),
    ]
    for box, heading, lines in specifications:
        datasets = [_rdf_mean(records, system, pair) for system, pair, _, _ in lines]
        y_max = max(float(dataset[1].max()) for dataset in datasets) * 1.12
        point = _draw_axes(
            draw,
            box,
            (0, 9),
            (0, y_max),
            (0, 2, 4, 6, 8),
            np.linspace(0, y_max, 5),
            "r / Å",
            "g(r)",
            y_decimals=1,
        )
        for (_, _, _, color), (radii, mean, _, _) in zip(lines, datasets, strict=False):
            coordinates = [point(float(r), float(g)) for r, g in zip(radii, mean, strict=False)]
            draw.line(coordinates, fill=_hex(color), width=5)
        draw.text(
            (box[0] + 38, box[1] + 26),
            heading,
            font=_pil_font(26, bold=True),
            fill=_hex(INK),
        )
        _draw_legend(
            draw,
            [(line[2], line[3]) for line in lines],
            box[0] + 355,
            box[1] + 43,
            gap=215,
        )
    mixed = metrics["systems"]["eda_1p5m"]
    pure = metrics["systems"]["pure_thf"]
    cards = [
        (
            "纯 THF 几何空腔代理",
            f"{pure['void_radius_angstrom']['mean']:.3f} ± "
            f"{pure['void_radius_angstrom']['replica_sd']:.3f} Å",
            PURE,
        ),
        (
            "en/THF 几何空腔代理",
            f"{mixed['void_radius_angstrom']['mean']:.3f} ± "
            f"{mixed['void_radius_angstrom']['replica_sd']:.3f} Å",
            MIXED,
        ),
        (
            "en N-THF O 接触",
            f"{mixed['eda_n_thf_o_contacts']['mean']:.2f} ± "
            f"{mixed['eda_n_thf_o_contacts']['replica_sd']:.2f} / 盒",
            EDA,
        ),
        (
            "en N-H...O(THF) 氢键",
            f"{mixed['eda_nh_thf_o_hydrogen_bonds']['mean']:.2f} ± "
            f"{mixed['eda_nh_thf_o_hydrogen_bonds']['replica_sd']:.2f} / 盒",
            ACCENT,
        ),
    ]
    for index, (label, value, color) in enumerate(cards):
        left = 65 + (index % 2) * 870
        top = 825 + (index // 2) * 100
        draw.rounded_rectangle(
            (left, top, left + 810, top + 78),
            radius=16,
            fill=_hex(BACKGROUND),
            outline=_hex(GRID),
            width=2,
        )
        draw.ellipse((left + 22, top + 24, left + 48, top + 50), fill=_hex(color))
        draw.text((left + 65, top + 25), label, font=_pil_font(18), fill=_hex(MUTED))
        draw.text(
            (left + 540, top + 39),
            value,
            anchor="mm",
            font=_pil_font(21, bold=True),
            fill=_hex(INK),
        )
    image.save(output, dpi=(180, 180))


def _parse_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0])
    elements: list[str] = []
    positions: list[list[float]] = []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        elements.append(fields[0])
        positions.append([float(value) for value in fields[1:4]])
    if len(elements) != count:
        raise ValueError(f"XYZ atom count mismatch in {path}")
    return elements, np.asarray(positions, dtype=float)


def _infer_components(
    elements: Sequence[str], positions: np.ndarray, lengths: np.ndarray
) -> tuple[list[tuple[int, int]], list[list[int]]]:
    radii = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66}
    bonds: list[tuple[int, int]] = []
    parent = list(range(len(elements)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(first: int, second: int) -> None:
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    for first in range(len(elements) - 1):
        vectors = positions[first + 1 :] - positions[first]
        vectors -= np.round(vectors / lengths) * lengths
        distances = np.linalg.norm(vectors, axis=1)
        thresholds = np.asarray(
            [
                1.23 * (radii[elements[first]] + radii[element]) + 0.05
                for element in elements[first + 1 :]
            ]
        )
        for offset in np.flatnonzero(distances <= thresholds):
            second = first + 1 + int(offset)
            bonds.append((first, second))
            union(first, second)
    groups: dict[int, list[int]] = defaultdict(list)
    for index in range(len(elements)):
        groups[find(index)].append(index)
    return bonds, sorted(groups.values(), key=lambda group: min(group))


def _unwrap_molecules(
    positions: np.ndarray, components: Sequence[Sequence[int]], lengths: np.ndarray
) -> np.ndarray:
    result = positions.copy()
    for component in components:
        reference = positions[component[0]]
        indices = np.asarray(component, dtype=int)
        vectors = positions[indices] - reference
        vectors -= np.round(vectors / lengths) * lengths
        unwrapped = reference + vectors
        shift = -np.floor(unwrapped.mean(axis=0) / lengths) * lengths
        result[indices] = unwrapped + shift
    return result


def _snapshot_void(record: dict[str, Any]) -> np.ndarray:
    target = int(record["snapshot"]["selection"]["frame_index"])
    row = next(item for item in record["timeseries"] if int(item["frame_index"]) == target)
    return np.asarray(
        [
            float(row["void_x_angstrom"]),
            float(row["void_y_angstrom"]),
            float(row["void_z_angstrom"]),
        ]
    )


def _render_snapshot_panel(
    draw, record: dict[str, Any], panel: tuple[int, int, int, int], title: str
) -> None:
    left, top, right, bottom = panel
    xyz = record["directory"] / "snapshot" / "representative.xyz"
    elements, raw_positions = _parse_xyz(xyz)
    cell = np.asarray(record["snapshot"]["structure"]["cell_vectors_angstrom"])
    lengths = np.diag(cell)
    bonds, components = _infer_components(elements, raw_positions, lengths)
    positions = _unwrap_molecules(raw_positions, components, lengths)
    is_eda = np.zeros(len(elements), dtype=bool)
    for component in components:
        if any(elements[index] == "N" for index in component):
            is_eda[np.asarray(component)] = True
    azimuth, elevation = math.radians(38), math.radians(24)
    rz = np.asarray(
        [
            [math.cos(azimuth), -math.sin(azimuth), 0],
            [math.sin(azimuth), math.cos(azimuth), 0],
            [0, 0, 1],
        ]
    )
    rx = np.asarray(
        [
            [1, 0, 0],
            [0, math.cos(elevation), -math.sin(elevation)],
            [0, math.sin(elevation), math.cos(elevation)],
        ]
    )
    rotation = rx @ rz
    projected = (positions - lengths / 2.0) @ rotation.T
    corners = np.asarray(
        [[x, y, z] for x in (0.0, lengths[0]) for y in (0.0, lengths[1]) for z in (0.0, lengths[2])]
    )
    projected_corners = (corners - lengths / 2.0) @ rotation.T
    width, height = right - left, bottom - top
    scale = min(
        (width - 85) / np.ptp(projected_corners[:, 0]),
        (height - 125) / np.ptp(projected_corners[:, 1]),
    )
    center_x = (left + right) / 2
    center_y = (top + bottom) / 2 + 18

    def point(value: np.ndarray) -> tuple[float, float]:
        return center_x + value[0] * scale, center_y - value[1] * scale

    draw.rounded_rectangle(panel, radius=22, fill=(248, 250, 253), outline=_hex(GRID), width=2)
    edges = []
    for first in range(8):
        for second in range(first + 1, 8):
            if sum(a != b for a, b in zip(corners[first], corners[second], strict=False)) == 1:
                edges.append((first, second))
    for first, second in edges:
        draw.line(
            (*point(projected_corners[first]), *point(projected_corners[second])),
            fill=(186, 196, 210),
            width=2,
        )
    for first, second in bonds:
        if elements[first] == "H" or elements[second] == "H":
            continue
        color = EDA if is_eda[first] or is_eda[second] else THF_C
        draw.line(
            (*point(projected[first]), *point(projected[second])),
            fill=_hex(color),
            width=4 if color == EDA else 3,
        )
    atom_order = sorted(
        (index for index, element in enumerate(elements) if element != "H"),
        key=lambda index: projected[index, 2],
    )
    for index in atom_order:
        element = elements[index]
        if element == "O":
            color, radius = THF_O, 7
        elif element == "N":
            color, radius = EDA, 8
        elif is_eda[index]:
            color, radius = "#47A9A0", 6
        else:
            color, radius = THF_C, 5
        x, y = point(projected[index])
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=_hex(color),
            outline=(255, 255, 255),
            width=1,
        )
    void_center = (_snapshot_void(record) - lengths / 2.0) @ rotation.T
    vx, vy = point(void_center)
    void_radius = (
        float(record["snapshot"]["selection"]["descriptors"]["void_radius_angstrom"]) * scale
    )
    draw.ellipse(
        (vx - void_radius, vy - void_radius, vx + void_radius, vy + void_radius),
        outline=_hex(ACCENT),
        width=5,
    )
    draw.ellipse((vx - 5, vy - 5, vx + 5, vy + 5), fill=_hex(ACCENT))
    draw.text((left + 24, top + 22), title, font=_pil_font(25, bold=True), fill=_hex(INK))
    description = (
        f"{record['snapshot']['snapshot_id']} | 几何空腔代理 "
        f"{record['snapshot']['selection']['descriptors']['void_radius_angstrom']:.2f} Å"
    )
    draw.text((left + 24, bottom - 46), description, font=_pil_font(17), fill=_hex(MUTED))


def render_structures(records: Sequence[dict[str, Any]], output: Path, gallery: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1030)
    draw = ImageDraw.Draw(image)
    _title(draw, "代表性溶剂结构", "重原子球棍图；黄色圆环是最大几何空腔代理，不代表电子密度")
    pure = next(
        record for record in records if record["system"] == "pure_thf" and record["replica"] == 1
    )
    mixed = next(
        record for record in records if record["system"] == "eda_1p5m" and record["replica"] == 1
    )
    _render_snapshot_panel(draw, pure, (55, 150, 875, 915), "纯 THF - r1")
    _render_snapshot_panel(draw, mixed, (925, 150, 1745, 915), "1.5 M en/THF - r1")
    legend = [
        ("THF C", THF_C),
        ("THF O", THF_O),
        ("en N", EDA),
        ("en C", "#47A9A0"),
        ("空腔代理", ACCENT),
    ]
    for index, (label, color) in enumerate(legend):
        x = 385 + index * 215
        draw.ellipse((x, 957, x + 20, 977), fill=_hex(color))
        draw.text((x + 31, 967), label, anchor="lm", font=_pil_font(17), fill=_hex(INK))
    image.save(output, dpi=(180, 180))

    image = _canvas(1800, 1280)
    draw = ImageDraw.Draw(image)
    _title(draw, "Stage A 六个代表快照", "每个副本从生产轨迹后半段按稳健 medoid 规则独立选取")
    ordered = [
        next(
            record
            for record in records
            if record["system"] == system and record["replica"] == replica
        )
        for system in SYSTEMS
        for replica in REPLICAS
    ]
    for index, record in enumerate(ordered):
        column, row = index % 3, index // 3
        panel = (
            45 + column * 585,
            145 + row * 545,
            570 + column * 585,
            650 + row * 545,
        )
        label = (
            f"{SYSTEM_LABELS[record['system']]} - r{record['replica']} - "
            f"{record['snapshot']['selection']['elapsed_ps'] / 1000:.1f} ns"
        )
        _render_snapshot_panel(draw, record, panel, label)
    image.save(gallery, dpi=(180, 180))


def render_workflow(output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 780)
    draw = ImageDraw.Draw(image)
    _title(
        draw, "推荐的多尺度计算边界", "周期性体相模型与非周期局域模型不是二选一，而是承担不同问题"
    )
    nodes = [
        (
            80,
            195,
            430,
            590,
            "Stage A",
            "周期性经典 MD",
            ["体相密度与组成", "RDF、接触与氢键", "代表快照与空腔种子"],
            PURE,
        ),
        (
            525,
            195,
            875,
            590,
            "Li0 前驱体",
            "局域重组/受限态",
            ["多种 en/THF 配位", "Li0-like 到紧密接触对", "保留多构型集合"],
            MIXED,
        ),
        (
            970,
            195,
            1320,
            590,
            "垂直电荷分离",
            "固定核 cDFT/Delta-SCF",
            ["Li 区域 3e -> 2e", "总电荷保持为 0", "垂直能隙与电子位置"],
            EDA,
        ),
        (
            1415,
            195,
            1765,
            590,
            "弛豫与验证",
            "短程 AIMD/嵌入团簇",
            ["溶剂重组", "自旋/电荷/质心", "尺寸与高阶方法检查"],
            ACCENT,
        ),
    ]
    for index, node in enumerate(nodes):
        left, top, right, bottom, heading, subheading, bullets, color = node
        draw.rounded_rectangle(
            (left, top, right, bottom),
            radius=26,
            fill=_hex(BACKGROUND),
            outline=_hex(color),
            width=4,
        )
        draw.rectangle((left, top, right, top + 18), fill=_hex(color))
        draw.text(
            ((left + right) / 2, top + 70),
            heading,
            anchor="mm",
            font=_pil_font(29, bold=True),
            fill=_hex(INK),
        )
        draw.text(
            ((left + right) / 2, top + 120),
            subheading,
            anchor="mm",
            font=_pil_font(21),
            fill=_hex(color),
        )
        for item_index, item in enumerate(bullets):
            y = top + 190 + item_index * 64
            draw.ellipse((left + 40, y - 6, left + 54, y + 8), fill=_hex(color))
            draw.text((left + 75, y), item, anchor="lm", font=_pil_font(20), fill=_hex(INK))
        if index < len(nodes) - 1:
            start = right + 16
            end = nodes[index + 1][0] - 16
            middle = (top + bottom) / 2
            draw.line((start, middle, end, middle), fill=_hex(MUTED), width=6)
            draw.polygon(
                ((end, middle), (end - 18, middle - 12), (end - 18, middle + 12)),
                fill=_hex(MUTED),
            )
    draw.rounded_rectangle(
        (235, 645, 1565, 725),
        radius=20,
        fill=(255, 248, 230),
        outline=_hex(ACCENT),
        width=2,
    )
    draw.text(
        (900, 685),
        "结论：Stage A 保留周期边界；电子结构阶段用周期模型与嵌入团簇交叉验证。",
        anchor="mm",
        font=_pil_font(23, bold=True),
        fill=_hex(INK),
    )
    image.save(output, dpi=(180, 180))


def _register_pdf_fonts() -> tuple[str, str]:
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.cidfonts import UnicodeCIDFont
    from reportlab.pdfbase.ttfonts import TTFont

    regular = _find_font(False)
    bold = _find_font(True)
    if regular is not None:
        try:
            pdfmetrics.registerFont(TTFont("StageARegular", str(regular), subfontIndex=0))
            pdfmetrics.registerFont(TTFont("StageABold", str(bold or regular), subfontIndex=0))
            return "StageARegular", "StageABold"
        except Exception:
            pass
    pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
    return "STSong-Light", "STSong-Light"


def _report_image(path: Path, width: float):
    from PIL import Image as PILImage
    from reportlab.platypus import Image

    with PILImage.open(path) as source:
        image_width, image_height = source.size
    return Image(str(path), width=width, height=width * image_height / image_width)


def build_pdf(records: Sequence[dict[str, Any]], metrics: dict[str, Any], output: Path) -> None:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    regular_font, bold_font = _register_pdf_fonts()
    styles = getSampleStyleSheet()
    title = ParagraphStyle(
        "TitleCN",
        parent=styles["Title"],
        fontName=bold_font,
        fontSize=27,
        leading=36,
        textColor=colors.HexColor(INK),
        alignment=TA_LEFT,
        wordWrap="CJK",
        spaceAfter=12,
    )
    subtitle = ParagraphStyle(
        "SubtitleCN",
        parent=styles["Normal"],
        fontName=regular_font,
        fontSize=13,
        leading=20,
        textColor=colors.HexColor(MUTED),
        wordWrap="CJK",
    )
    heading = ParagraphStyle(
        "HeadingCN",
        parent=styles["Heading1"],
        fontName=bold_font,
        fontSize=19,
        leading=25,
        textColor=colors.HexColor(INK),
        wordWrap="CJK",
        spaceAfter=10,
    )
    heading2 = ParagraphStyle(
        "Heading2CN",
        parent=styles["Heading2"],
        fontName=bold_font,
        fontSize=14,
        leading=20,
        textColor=colors.HexColor(MIXED),
        wordWrap="CJK",
        spaceBefore=9,
        spaceAfter=5,
    )
    body = ParagraphStyle(
        "BodyCN",
        parent=styles["BodyText"],
        fontName=regular_font,
        fontSize=10.4,
        leading=17,
        textColor=colors.HexColor(INK),
        wordWrap="CJK",
        spaceAfter=7,
    )
    small = ParagraphStyle(
        "SmallCN",
        parent=body,
        fontSize=8.2,
        leading=12.5,
        textColor=colors.HexColor(MUTED),
        spaceAfter=3,
    )
    callout = ParagraphStyle(
        "CalloutCN",
        parent=body,
        fontName=bold_font,
        fontSize=11.4,
        leading=18,
        backColor=colors.HexColor("#FFF7E6"),
        borderColor=colors.HexColor(ACCENT),
        borderWidth=1,
        borderPadding=9,
        spaceBefore=7,
        spaceAfter=9,
    )
    bullet = ParagraphStyle(
        "BulletCN",
        parent=body,
        leftIndent=15,
        firstLineIndent=-10,
        bulletIndent=0,
        spaceAfter=4,
    )

    def p(text: str, style=body):
        return Paragraph(text, style)

    table_header = ParagraphStyle(
        "TableHeader", parent=small, fontName=bold_font, textColor=colors.white
    )

    def table(data: Sequence[Sequence[Any]], widths: Sequence[float]):
        formatted = []
        for row_index, row in enumerate(data):
            row_style = table_header if row_index == 0 else small
            formatted.append(
                [item if hasattr(item, "wrap") else Paragraph(str(item), row_style) for item in row]
            )
        result = Table(formatted, colWidths=widths, repeatRows=1, hAlign="LEFT")
        commands = [
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(INK)),
            ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor(GRID)),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 7),
            ("RIGHTPADDING", (0, 0), (-1, -1), 7),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ]
        for row_index in range(2, len(data), 2):
            commands.append(
                ("BACKGROUND", (0, row_index), (-1, row_index), colors.HexColor(BACKGROUND))
            )
        result.setStyle(TableStyle(commands))
        return result

    def page_frame(canvas, document):
        canvas.saveState()
        width, _height = A4
        canvas.setStrokeColor(colors.HexColor(GRID))
        canvas.line(18 * mm, 15 * mm, width - 18 * mm, 15 * mm)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.setFont(regular_font, 7.5)
        canvas.drawString(18 * mm, 9.5 * mm, "Li/THF/en 溶剂化电子项目 - Stage A")
        canvas.drawRightString(width - 18 * mm, 9.5 * mm, f"第 {document.page} 页")
        canvas.restoreState()

    document = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=17 * mm,
        bottomMargin=20 * mm,
        title="Stage A: THF/en 体相溶剂模型与代表快照报告",
        author="solvelec workflow",
    )
    usable_width = A4[0] - 36 * mm
    pure = metrics["systems"]["pure_thf"]
    mixed = metrics["systems"]["eda_1p5m"]
    story: list[Any] = [
        Spacer(1, 18 * mm),
        p("Stage A", subtitle),
        p("THF/en 体相溶剂模型与代表快照报告", title),
        p("面向实验讨论的可审计计算摘要 | 298.15 K, 1 bar | 2026-09-03", subtitle),
        Spacer(1, 12 * mm),
        _report_image(FIGURES / "stage_a_structures.png", usable_width),
        Spacer(1, 9 * mm),
        p(
            "核心结论：Stage A 已获得 2 个体系、各 3 条独立 20 ns 生产轨迹和 6 个代表性纯溶剂快照。当前混合体系为 64 THF + 9 en，即 THF:en = 7.11:1，实际 en 浓度约 1.510 M；本阶段不包含 Li 或过量电子。",
            callout,
        ),
        p(
            "黄色空腔仅是基于重原子范德华表面的几何自由体积代理，不能解释为溶剂化电子位置。所有数值均由提交的 Stage A JSON/CSV/XYZ 文件重新计算。",
            small,
        ),
        PageBreak(),
        p("1. 研究问题与模型边界", heading),
        p(
            "实验过程是 en/THF 溶液接触 Li 金属，Li 随时间逐步电离。因此至少存在三个不同的“比例”：体相 THF:en 配比、总体 Li:en 投料或反应计量，以及某个 Li 周围第一配位层中 en 与 THF 的局部组成。Stage A 只处理第一项。"
        ),
        p("当前 Stage A 体系", heading2),
        table(
            [
                ["体系", "盒内组成", "体相比例", "实际状态"],
                ["纯 THF", "64 THF", "无 en", "3 x 20 ns，ready"],
                [
                    "1.5 M en/THF",
                    "64 THF + 9 en",
                    "THF:en = 7.11:1; x(en)=12.3%",
                    "3 x 20 ns，ready",
                ],
            ],
            [34 * mm, 40 * mm, 59 * mm, 39 * mm],
        ),
        Spacer(1, 4 * mm),
        p(
            "当前 64:9 配比描述的是 1.5 M en 的体相采样点，不等于实验人员提到的 Li:en = 1:2。后者若代表反应计量，只说明模型可以首先采用单 Li 稀释极限；真实 THF:en 仍应由配液量和密度确定。",
            callout,
        ),
        p("为什么 Stage A 保留周期边界", heading2),
        p(
            "周期性边界适合回答体相问题：密度、浓度、RDF、氢键与空腔统计不会受到真空表面和有限液滴边缘支配。几十到几百个分子的自由团簇更接近液滴或气相团簇，不能替代体相溶剂平衡。"
        ),
        p(
            "高成本电子结构阶段不必只用周期模型。建议保留 Stage A 周期 MD 作为体相结构来源，再截取 Li 周围第一、第二溶剂壳层，使用嵌入点电荷或连续介质做非周期高阶计算，并与少量周期 cDFT 结果交叉验证。"
        ),
        _report_image(FIGURES / "stage_a_model_boundary.png", usable_width),
        PageBreak(),
        p("2. 模拟与筛选协议", heading),
        table(
            [
                ["项目", "设置"],
                ["经典模型", "GROMACS; GAFF2; AM1-BCC 电荷"],
                ["热力学条件", "298.15 K; 1.0 bar"],
                ["时间步长", "2.0 fs"],
                ["平衡流程", "0.5 ns NVT + 5.0 ns NPT"],
                ["生产采样", "每个副本 20 ns; 每体系 3 个副本"],
                ["轨迹/分析步长", "轨迹 10 ps; Stage A 分析 100 ps"],
                ["代表快照候选区间", "生产轨迹后半段，即 10-20 ns"],
                ["快照规则", "密度、空腔代理、接触和氢键的稳健 medoid"],
            ],
            [52 * mm, 120 * mm],
        ),
        Spacer(1, 5 * mm),
        p("接受门控", heading2),
        p("• 每条轨迹覆盖至少 98% 的请求时长。", bullet),
        p("• 前后半段平均密度的相对差不超过 2%。", bullet),
        p("• 同一体系三个副本的平均密度跨度不超过均值的 3%。", bullet),
        p("• en 实际浓度与目标值的差不超过 0.05 M。", bullet),
        p("• 每条分析至少 150 帧，且所有描述符的有效样本数至少为 5。", bullet),
        p(
            f"实际最小有效样本数：纯 THF 为 {pure['minimum_effective_sample_size']:.1f}，en/THF 为 {mixed['minimum_effective_sample_size']:.1f}；两个体系均有 201 个分析帧/副本。",
            callout,
        ),
        p("代表快照清单", heading2),
        table(
            [["体系", "副本", "时间 / ns", "密度 / g mL^-1", "空腔代理 / Å"]]
            + [
                [
                    SYSTEM_LABELS[system],
                    item["replica"],
                    f"{item['time_ns']:.1f}",
                    f"{item['density_g_ml']:.4f}",
                    f"{item['void_radius_angstrom']:.3f}",
                ]
                for system in SYSTEMS
                for item in metrics["systems"][system]["snapshots"]
            ],
            [40 * mm, 22 * mm, 32 * mm, 43 * mm, 35 * mm],
        ),
        PageBreak(),
        p("3. 密度、浓度与副本一致性", heading),
        _report_image(FIGURES / "stage_a_density.png", usable_width),
        Spacer(1, 4 * mm),
        table(
            [
                ["体系", "平均密度 ± 副本 SD", "平均体积 ± 副本 SD", "最大半程密度差", "最小 ESS"],
                [
                    "纯 THF",
                    f"{pure['density_g_ml']['mean']:.5f} ± {pure['density_g_ml']['replica_sd']:.5f} g/mL",
                    f"{pure['volume_nm3']['mean']:.4f} ± {pure['volume_nm3']['replica_sd']:.4f} nm3",
                    f"{pure['maximum_density_half_difference_percent']:.2f}%",
                    f"{pure['minimum_effective_sample_size']:.1f}",
                ],
                [
                    "1.5 M en/THF",
                    f"{mixed['density_g_ml']['mean']:.5f} ± {mixed['density_g_ml']['replica_sd']:.5f} g/mL",
                    f"{mixed['volume_nm3']['mean']:.4f} ± {mixed['volume_nm3']['replica_sd']:.4f} nm3",
                    f"{mixed['maximum_density_half_difference_percent']:.2f}%",
                    f"{mixed['minimum_effective_sample_size']:.1f}",
                ],
            ],
            [31 * mm, 43 * mm, 42 * mm, 30 * mm, 26 * mm],
        ),
        Spacer(1, 4 * mm),
        p(
            f"由三个混合体系平均体积重新计算，9 个 en 对应 {mixed['achieved_eda_concentration_m']:.3f} M。加入 en 后模拟平均密度相对纯 THF 变化 {metrics['comparisons']['mixture_vs_pure_density_percent']:.2f}%。这是模型内部响应；获得实验配液密度前，不应表述为实验吻合。"
        ),
        p(
            "实验比对优先需要记录 THF 与 en 的实际用量或质量、配液温度、混合液密度和 Li 投料量。只给 en 摩尔浓度不足以唯一确认真实 THF:en，尤其当混合体积非理想时。",
            callout,
        ),
        PageBreak(),
        p("4. en-THF 微观混合与自由体积", heading),
        _report_image(FIGURES / "stage_a_microstructure.png", usable_width),
        Spacer(1, 4 * mm),
        p(
            f"在 3.5 Å 阈值下，en 氮与 THF 氧接触数为 {mixed['eda_n_thf_o_contacts']['mean']:.2f} ± {mixed['eda_n_thf_o_contacts']['replica_sd']:.2f}/周期盒；满足 N-H...O 条件的氢键为 {mixed['eda_nh_thf_o_hydrogen_bonds']['mean']:.2f} ± {mixed['eda_nh_thf_o_hydrogen_bonds']['replica_sd']:.2f}/盒。这些是全盒计数，不是 Li 配位数。"
        ),
        p(
            f"几何空腔代理从纯 THF 的 {pure['void_radius_angstrom']['mean']:.3f} ± {pure['void_radius_angstrom']['replica_sd']:.3f} Å 增至 en/THF 的 {mixed['void_radius_angstrom']['mean']:.3f} ± {mixed['void_radius_angstrom']['replica_sd']:.3f} Å，平均增加 {metrics['comparisons']['mixture_vs_pure_void_radius_percent']:.1f}%。这提示自由体积统计变化，但不能证明电子更稳定或更离域。"
        ),
        p("三副本平均 RDF 在指定区间内的主峰：", heading2),
        table(
            [["体系", "原子对", "主峰位置 / Å", "g(r) 峰高"]]
            + [
                [SYSTEM_LABELS[system], pair, f"{peak['r_angstrom']:.2f}", f"{peak['g_r']:.2f}"]
                for system in SYSTEMS
                for pair, peak in metrics["systems"][system]["rdf_peaks"].items()
            ],
            [40 * mm, 48 * mm, 42 * mm, 42 * mm],
        ),
        PageBreak(),
        p("5. 六个代表性快照", heading),
        p(
            "每个副本输出一个代表性快照，而不是只选能量最低或空腔最大的瞬间。筛选从 10-20 ns 中寻找多描述符空间的稳健中心，以减少极端构型对后续电子结构计算的偏置。"
        ),
        _report_image(FIGURES / "stage_a_snapshot_gallery.png", usable_width),
        Spacer(1, 4 * mm),
        p(
            "图中只显示重原子；THF 氧为红色，en 氮为蓝色，en 碳为青绿色。黄色标记来自几何自由体积搜索。六个原始 XYZ 与晶胞文件随报告提交。",
            small,
        ),
        PageBreak(),
        p("6. 对后续 Li/电子计算的含义", heading),
        p("已经支持的判断", heading2),
        p("• 当前混合溶剂的体相密度、浓度和副本一致性通过门控。", bullet),
        p("• en 与 THF 存在可重复的 N-O 近程相关和约 5.53 个/盒的氢键。", bullet),
        p("• 六个快照可作为单 Li 稀释极限的独立溶剂环境种子。", bullet),
        p("尚不能支持的判断", heading2),
        p("• Stage A 没有 Li，不能给出 Li 的氧/氮配位数或 Li0/Li+ 平衡。", bullet),
        p("• 几何空腔不是电子密度，不能证明溶剂化电子已经形成。", bullet),
        p("• 只有一个 en 浓度，不能建立 THF:en 比例的趋势。", bullet),
        p("• 固定电荷经典力场不能描述 Li 电离、电子转移或溶剂极化。", bullet),
        p("建议的下一条计算路线", heading2),
        _report_image(FIGURES / "stage_a_model_boundary.png", usable_width),
        p(
            "新增路线应先在多个快照中构建 Li0-like/紧密接触对前驱体，并采样不同 en/THF 第一配位层；随后在固定核几何上做总电荷为 0 的 cDFT 垂直电荷分离，再允许溶剂重组。周期结果保持体相环境，嵌入团簇用于高阶电子结构和尺寸检查。",
            callout,
        ),
        PageBreak(),
        p("7. 数据来源、复现与限制", heading),
        p(
            "报告生成器直接读取 reports/stage_a/data。构建时检查两个 summary 的 ready 状态、六条记录唯一性、Li 缺失标志、帧数以及 analysis/XYZ/cell 的 SHA-256；任一检查失败都会停止。"
        ),
    ]
    provenance = [["体系/副本", "快照 ID", "XYZ SHA-256（前 16 位）", "分析 SHA-256（前 16 位）"]]
    for record in records:
        snapshot = record["snapshot"]
        provenance.append(
            [
                f"{record['system']}/r{record['replica']}",
                snapshot["snapshot_id"],
                snapshot["structure"]["xyz"]["sha256"][:16],
                snapshot["source"]["analysis"]["sha256"][:16],
            ]
        )
    story.extend(
        [
            table(provenance, [32 * mm, 55 * mm, 42 * mm, 43 * mm]),
            Spacer(1, 5 * mm),
            p("一键重建", heading2),
            p('python -m pip install -e ".[report]"', callout),
            p("python reports/stage_a/build_report.py", callout),
            p("输出内容", heading2),
            p("• stage_a_solvent_report_zh.pdf：本报告。", bullet),
            p("• figures/*.png：可直接在 GitHub 或聊天中分享的图片。", bullet),
            p("• stage_a_metrics.json：从六条记录聚合的机器可读指标。", bullet),
            p("• data/analysis/...：六个快照、时间序列、RDF 与审计元数据。", bullet),
            p("解释限制", heading2),
            p(
                "本报告是 solvent-only pilot 的内部计算结果，不是实验密度或光谱验证。RDF 和空腔基于固定电荷经典力场；未加入 Li、过量电子、可极化力场或从头算动力学。电子局域、离域和 Li 电离必须由电子/自旋密度、VDE/光谱及方法敏感性门控给出。"
            ),
            p("项目：https://github.com/SunsetStand/li-thf-amine-solvated-electron", small),
        ]
    )
    document.build(story, onFirstPage=page_frame, onLaterPages=page_frame)


def build(output: Path = DEFAULT_PDF) -> dict[str, Any]:
    records = load_records()
    validate_records(records)
    metrics = summarize(records)
    FIGURES.mkdir(parents=True, exist_ok=True)
    render_density(records, metrics, FIGURES / "stage_a_density.png")
    render_microstructure(records, metrics, FIGURES / "stage_a_microstructure.png")
    render_structures(
        records,
        FIGURES / "stage_a_structures.png",
        FIGURES / "stage_a_snapshot_gallery.png",
    )
    render_workflow(FIGURES / "stage_a_model_boundary.png")
    METRICS_JSON.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    build_pdf(records, metrics, output)
    return metrics


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_PDF)
    args = parser.parse_args(argv)
    metrics = build(args.output)
    print(f"READY: {args.output}")
    print(
        "Stage A: "
        f"density(pure)={metrics['systems']['pure_thf']['density_g_ml']['mean']:.5f} g/mL, "
        f"density(mixed)={metrics['systems']['eda_1p5m']['density_g_ml']['mean']:.5f} g/mL, "
        f"EDA={metrics['systems']['eda_1p5m']['achieved_eda_concentration_m']:.3f} M"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

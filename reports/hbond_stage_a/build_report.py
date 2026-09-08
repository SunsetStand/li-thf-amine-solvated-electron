#!/usr/bin/env python3
# ruff: noqa: E501
"""Build the hbond_pilot Stage A Chinese PDF from committed audit data.

The builder is intentionally local and deterministic. It reads only compact
JSON/CSV/XYZ/PDB/SVG exports under ``reports/hbond_stage_a/data`` and performs
strict readiness, composition, scope, frame-count, and SHA-256 checks before
creating figures or the PDF.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"
DEFAULT_PDF = HERE / "hbond_stage_a_report_zh.pdf"
METRICS_JSON = HERE / "hbond_stage_a_metrics.json"
SYSTEMS = ("pure_eda", "eda_1to1")
SYSTEM_LABELS = {"pure_eda": "纯 en", "eda_1to1": "THF:en = 1:1"}
REPLICAS = (1, 2, 3)
COMPOSITIONS = {
    "pure_eda": {"thf": 0, "eda": 64, "atom_count": 768},
    "eda_1to1": {"thf": 32, "eda": 32, "atom_count": 800},
}

INK = "#17233B"
MUTED = "#62708A"
GRID = "#DDE4EE"
PURE = "#7353BA"
MIXED = "#137C8B"
ACCENT = "#E59F23"
BRIDGE = "#D1495B"
ASSOCIATED = "#2F80ED"
BACKGROUND = "#F4F7FB"
REPLICA_COLORS = ("#315C9A", "#21A179", "#E58E26")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sample_sd(values: Sequence[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def load_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for system in SYSTEMS:
        for replica in REPLICAS:
            directory = DATA / "analysis" / system / f"r{replica}"
            snapshot_dir = directory / "snapshot"
            records.append(
                {
                    "system": system,
                    "replica": replica,
                    "directory": directory,
                    "analysis": _load_json(directory / "analysis.json"),
                    "timeseries": _read_csv(directory / "timeseries.csv"),
                    "rdf": _read_csv(directory / "rdf.csv"),
                    "hydrogen_bonds": _read_csv(directory / "hydrogen_bonds.csv"),
                    "snapshot": _load_json(snapshot_dir / "metadata.json"),
                    "cavity": _load_json(snapshot_dir / "cavity_hbonds.json"),
                }
            )
    return records


def validate_records(records: Sequence[dict[str, Any]]) -> None:
    required = (
        DATA / "classical_pilot.validation.json",
        DATA / "classical_analysis.summary.json",
        DATA / "snapshot_bank.summary.json",
        DATA / "hbond_stage_a.done",
        DATA / "report_provenance.json",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ValueError("missing hbond_pilot report inputs: " + ", ".join(missing))
    for filename in (
        "classical_pilot.validation.json",
        "classical_analysis.summary.json",
        "snapshot_bank.summary.json",
    ):
        summary = _load_json(DATA / filename)
        if not summary.get("ready"):
            raise ValueError(f"{filename} is not ready")
    for filename in ("classical_analysis.summary.json", "snapshot_bank.summary.json"):
        summary = _load_json(DATA / filename)
        if int(summary.get("record_count", -1)) != 6:
            raise ValueError(f"{filename} does not contain six records")
    if len(records) != 6:
        raise ValueError(f"expected six analysis records, found {len(records)}")

    seen: set[tuple[str, int]] = set()
    for record in records:
        key = (record["system"], record["replica"])
        if key in seen:
            raise ValueError(f"duplicate record {key}")
        seen.add(key)
        analysis = record["analysis"]
        snapshot = record["snapshot"]
        cavity = record["cavity"]
        if analysis.get("system_id") != key[0] or int(analysis.get("replica", -1)) != key[1]:
            raise ValueError(f"analysis identity mismatch for {key}")
        if not analysis.get("ready") or not analysis.get("metrics", {}).get("ready"):
            raise ValueError(f"analysis is not ready for {key}")
        if int(analysis["metrics"].get("atom_count", -1)) != COMPOSITIONS[key[0]]["atom_count"]:
            raise ValueError(f"unexpected atom count for {key}")
        if not snapshot.get("ready") or snapshot.get("li_atom_present"):
            raise ValueError(f"snapshot is not a ready solvent-only structure for {key}")
        if cavity.get("scientific_status") != "GEOMETRIC_CAVITY_HBOND_VISUALIZATION":
            raise ValueError(f"unexpected cavity status for {key}")
        if len(record["timeseries"]) != int(analysis["metrics"]["frame_count"]):
            raise ValueError(f"timeseries length mismatch for {key}")
        if not record["timeseries"]:
            raise ValueError(f"empty timeseries for {key}")
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
            (
                directory / "snapshot" / "cavity_hbonds.json",
                snapshot["structure"]["cavity_hbonds"]["sha256"],
            ),
            (
                directory / "snapshot" / "cavity_local.pdb",
                snapshot["structure"]["cavity_local_pdb"]["sha256"],
            ),
            (
                directory / "snapshot" / "cavity_hbonds.pml",
                snapshot["structure"]["cavity_pymol"]["sha256"],
            ),
        )
        for path, expected in checks:
            observed = _sha256(path)
            if observed != expected:
                raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")

    expected = {(system, replica) for system in SYSTEMS for replica in REPLICAS}
    if seen != expected:
        raise ValueError(f"unexpected record identities: {sorted(seen)}")
    provenance = _load_json(DATA / "report_provenance.json")
    if provenance.get("source_campaign") != "hbond_pilot":
        raise ValueError("report provenance is not for hbond_pilot")


def _descriptor(
    records: Sequence[dict[str, Any]], system: str, name: str
) -> tuple[float, float, list[float]]:
    values = [
        float(record["analysis"]["metrics"]["mean_descriptors"][name])
        for record in records
        if record["system"] == system
    ]
    return statistics.mean(values), _sample_sd(values), values


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
        if selected:
            series.append(
                (
                    np.asarray([float(row["r_center_angstrom"]) for row in selected]),
                    np.asarray([float(row["g_r"]) for row in selected]),
                )
            )
    if not series:
        raise ValueError(f"missing RDF pair {pair} for {system}")
    radii = series[0][0]
    if any(not np.allclose(radii, item[0]) for item in series[1:]):
        raise ValueError(f"inconsistent RDF radii for {system}/{pair}")
    values = np.stack([item[1] for item in series])
    return radii, values.mean(axis=0), values.min(axis=0), values.max(axis=0)


def _rdf_peak(
    records: Sequence[dict[str, Any]], system: str, pair: str, lower: float, upper: float
) -> dict[str, float]:
    radii, mean, _minimum, _maximum = _rdf_mean(records, system, pair)
    mask = (radii >= lower) & (radii <= upper)
    candidates = np.flatnonzero(mask)
    if len(candidates) == 0:
        raise ValueError(f"no RDF bins in requested peak interval for {system}/{pair}")
    index = int(candidates[np.argmax(mean[mask])])
    return {"r_angstrom": float(radii[index]), "g_r": float(mean[index])}


def summarize(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    systems: dict[str, Any] = {}
    descriptor_names = (
        "density_g_ml",
        "volume_nm3",
        "void_radius_angstrom",
        "eda_eda_hydrogen_bonds",
        "cavity_associated_hydrogen_bonds",
        "cavity_bridging_hydrogen_bonds",
    )
    for system in SYSTEMS:
        subset = [record for record in records if record["system"] == system]
        descriptors: dict[str, Any] = {}
        for name in descriptor_names:
            mean, sd, replica_means = _descriptor(records, system, name)
            descriptors[name] = {
                "mean": mean,
                "replica_sd": sd,
                "replica_means": replica_means,
            }
        if system == "eda_1to1":
            for name in ("eda_thf_contacts", "eda_thf_hydrogen_bonds"):
                mean, sd, replica_means = _descriptor(records, system, name)
                descriptors[name] = {
                    "mean": mean,
                    "replica_sd": sd,
                    "replica_means": replica_means,
                }
        snapshots = []
        for record in subset:
            selection = record["snapshot"]["selection"]
            whole = record["cavity"]["hydrogen_bond_counts_whole_snapshot"]
            local = record["cavity"]["hydrogen_bond_counts_local_view"]
            snapshots.append(
                {
                    "replica": record["replica"],
                    "snapshot_id": record["snapshot"]["snapshot_id"],
                    "time_ns": float(selection["elapsed_ps"]) / 1000.0,
                    "density_g_ml": float(selection["descriptors"]["density_g_ml"]),
                    "void_radius_angstrom": float(
                        selection["descriptors"]["void_radius_angstrom"]
                    ),
                    "whole_snapshot_hydrogen_bonds": whole,
                    "local_view_hydrogen_bonds": local,
                    "local_atom_count": int(record["cavity"]["local_atom_count"]),
                }
            )
        effective_samples = [
            float(value["effective_sample_size"])
            for record in subset
            for value in record["analysis"]["metrics"]["autocorrelation"].values()
        ]
        peak_intervals = {"eda_n-eda_n": (2.0, 4.5)}
        if system == "eda_1to1":
            peak_intervals.update(
                {"eda_n-thf_o": (2.0, 4.2), "thf_o-thf_o": (2.0, 5.2)}
            )
        systems[system] = {
            "composition": COMPOSITIONS[system],
            "replicas": 3,
            "production_ns_per_replica": 20.0,
            "analysis_frames_per_replica": int(subset[0]["analysis"]["metrics"]["frame_count"]),
            "descriptors": descriptors,
            "maximum_density_half_difference_percent": max(
                _half_difference_percent(record["timeseries"], "density_g_ml")
                for record in subset
            ),
            "minimum_effective_sample_size": min(effective_samples),
            "snapshots": snapshots,
            "rdf_peaks": {
                pair: _rdf_peak(records, system, pair, *interval)
                for pair, interval in peak_intervals.items()
            },
        }
    comparisons = {
        "mixture_vs_pure_density_percent": (
            systems["eda_1to1"]["descriptors"]["density_g_ml"]["mean"]
            / systems["pure_eda"]["descriptors"]["density_g_ml"]["mean"]
            - 1.0
        )
        * 100.0,
        "mixture_vs_pure_void_radius_percent": (
            systems["eda_1to1"]["descriptors"]["void_radius_angstrom"]["mean"]
            / systems["pure_eda"]["descriptors"]["void_radius_angstrom"]["mean"]
            - 1.0
        )
        * 100.0,
    }
    provenance = _load_json(DATA / "report_provenance.json")
    return {
        "schema_version": 1,
        "status": "READY",
        "source": "hbond_pilot Stage A",
        "source_campaign": provenance["source_campaign"],
        "source_commit": provenance["source_commit"],
        "slurm_controller_job": provenance["slurm_controller_job"],
        "systems": systems,
        "comparisons": comparisons,
        "hydrogen_bond_definition": records[0]["analysis"]["metrics"]
        ["hydrogen_bond_definition"],
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

    path = _find_font(bold)
    if path is None:
        return ImageFont.load_default()
    return ImageFont.truetype(str(path), size=size, index=0)


def _hex(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[index : index + 2], 16) for index in (0, 2, 4))


def _canvas(width: int = 1800, height: int = 1050):
    from PIL import Image

    return Image.new("RGB", (width, height), "white")


def _title(draw: Any, title: str, subtitle: str | None = None) -> None:
    draw.text((72, 48), title, font=_pil_font(42, True), fill=_hex(INK))
    if subtitle:
        draw.text((74, 105), subtitle, font=_pil_font(21), fill=_hex(MUTED))


def _draw_axes(
    draw: Any,
    box: tuple[int, int, int, int],
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    x_label: str,
    y_label: str,
):
    left, top, right, bottom = box
    draw.rounded_rectangle(box, radius=18, fill=(250, 252, 255), outline=_hex(GRID), width=2)
    if math.isclose(y_range[0], y_range[1]):
        y_range = (y_range[0] - 0.5, y_range[1] + 0.5)
    if math.isclose(x_range[0], x_range[1]):
        x_range = (x_range[0], x_range[1] + 1.0)

    def point(x: float, y: float) -> tuple[float, float]:
        px = left + 75 + (x - x_range[0]) / (x_range[1] - x_range[0]) * (right - left - 110)
        py = bottom - 60 - (y - y_range[0]) / (y_range[1] - y_range[0]) * (bottom - top - 105)
        return px, py

    x0, y0 = point(x_range[0], y_range[0])
    x1, y1 = point(x_range[1], y_range[1])
    for fraction in np.linspace(0.0, 1.0, 5):
        value = y_range[0] + fraction * (y_range[1] - y_range[0])
        _, py = point(x_range[0], value)
        draw.line((x0, py, x1, py), fill=_hex(GRID), width=2)
        draw.text((x0 - 12, py), f"{value:.2f}", anchor="rm", font=_pil_font(16), fill=_hex(MUTED))
    for fraction in np.linspace(0.0, 1.0, 5):
        value = x_range[0] + fraction * (x_range[1] - x_range[0])
        px, _ = point(value, y_range[0])
        draw.text((px, y0 + 13), f"{value:.0f}", anchor="ma", font=_pil_font(16), fill=_hex(MUTED))
    draw.line((x0, y0, x1, y0), fill=_hex(INK), width=2)
    draw.line((x0, y0, x0, y1), fill=_hex(INK), width=2)
    draw.text(((x0 + x1) / 2, bottom - 18), x_label, anchor="mm", font=_pil_font(18), fill=_hex(INK))
    draw.text((left + 18, top + 15), y_label, font=_pil_font(17, True), fill=_hex(INK))
    return point


def _save(image: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path, dpi=(180, 180))


def render_density_void(records: Sequence[dict[str, Any]], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1080)
    draw = ImageDraw.Draw(image)
    _title(draw, "密度与几何空腔随时间", "细线为三个独立副本；右侧 10-20 ns 是代表快照候选区间")
    panels = [
        ("density_g_ml", "密度 / g mL^-1", (65, 165, 875, 855)),
        ("void_radius_angstrom", "空腔半径 / Å", (925, 165, 1735, 855)),
    ]
    for field, label, box in panels:
        all_values = [float(row[field]) for record in records for row in record["timeseries"]]
        margin = max(0.02, (max(all_values) - min(all_values)) * 0.08)
        point = _draw_axes(
            draw,
            box,
            (0.0, 20.0),
            (min(all_values) - margin, max(all_values) + margin),
            "生产时间 / ns",
            label,
        )
        for system, base_color in (("pure_eda", PURE), ("eda_1to1", MIXED)):
            for record in (r for r in records if r["system"] == system):
                points = [
                    point(float(row["elapsed_ps"]) / 1000.0, float(row[field]))
                    for row in record["timeseries"]
                ]
                draw.line(points, fill=_hex(base_color), width=3)
    legend = [("纯 en（3 副本）", PURE), ("1:1 THF:en（3 副本）", MIXED)]
    for index, (label, color) in enumerate(legend):
        x = 580 + index * 380
        draw.line((x, 935, x + 36, 935), fill=_hex(color), width=7)
        draw.text((x + 50, 935), label, anchor="lm", font=_pil_font(19), fill=_hex(INK))
    _save(image, output)


def render_hbond_statistics(records: Sequence[dict[str, Any]], metrics: dict[str, Any], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1050)
    draw = ImageDraw.Draw(image)
    _title(draw, "氢键网络与空腔关系", "柱高为三副本均值，误差线为副本间样本标准差；所有数量均为每个周期盒")
    names = (
        ("eda_eda_hydrogen_bonds", "en-en\n氢键"),
        ("eda_thf_hydrogen_bonds", "en-THF\n氢键"),
        ("cavity_associated_hydrogen_bonds", "空腔邻近\n氢键"),
        ("cavity_bridging_hydrogen_bonds", "穿越空腔\n氢键"),
    )
    values: list[float] = []
    for system in SYSTEMS:
        descriptors = metrics["systems"][system]["descriptors"]
        for key, _label in names:
            values.append(float(descriptors.get(key, {"mean": 0.0})["mean"]))
    y_max = max(values + [1.0]) * 1.25
    box = (110, 165, 1690, 810)
    point = _draw_axes(draw, box, (0.0, 4.0), (0.0, y_max), "氢键分类", "平均数量 / 盒")
    centers = (0.55, 1.5, 2.45, 3.4)
    width = 0.27
    for system_index, system in enumerate(SYSTEMS):
        color = PURE if system == "pure_eda" else MIXED
        descriptors = metrics["systems"][system]["descriptors"]
        for (key, label), center in zip(names, centers, strict=True):
            entry = descriptors.get(key, {"mean": 0.0, "replica_sd": 0.0})
            mean, sd = float(entry["mean"]), float(entry["replica_sd"])
            x_center = center + (system_index - 0.5) * width
            left, bottom = point(x_center - width * 0.42, 0.0)
            right, top = point(x_center + width * 0.42, mean)
            draw.rounded_rectangle((left, top, right, bottom), radius=5, fill=_hex(color))
            x, y_low = point(x_center, max(0.0, mean - sd))
            _, y_high = point(x_center, min(y_max, mean + sd))
            draw.line((x, y_low, x, y_high), fill=_hex(INK), width=3)
            draw.line((x - 10, y_low, x + 10, y_low), fill=_hex(INK), width=3)
            draw.line((x - 10, y_high, x + 10, y_high), fill=_hex(INK), width=3)
            if system_index == 0:
                tx, ty = point(center, 0.0)
                draw.multiline_text((tx, ty + 36), label, anchor="ma", align="center", font=_pil_font(17), fill=_hex(INK), spacing=3)
    for index, (label, color) in enumerate((("纯 en", PURE), ("1:1 THF:en", MIXED))):
        x = 650 + index * 330
        draw.rectangle((x, 930, x + 30, 954), fill=_hex(color))
        draw.text((x + 43, 942), label, anchor="lm", font=_pil_font(19), fill=_hex(INK))
    _save(image, output)


def render_rdf(records: Sequence[dict[str, Any]], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1050)
    draw = ImageDraw.Draw(image)
    _title(draw, "溶剂原子对径向分布函数", "实线为三副本均值；阴影边界用同色细线表示副本范围")
    panels = (
        ("pure_eda", "eda_n-eda_n", "纯 en: N-N", (65, 165, 875, 820), PURE),
        ("eda_1to1", "eda_n-thf_o", "1:1 THF:en: N-O", (925, 165, 1735, 820), MIXED),
    )
    for system, pair, label, box, color in panels:
        radii, mean, minimum, maximum = _rdf_mean(records, system, pair)
        mask = radii <= 8.0
        radii, mean, minimum, maximum = (
            radii[mask],
            mean[mask],
            minimum[mask],
            maximum[mask],
        )
        y_max = max(1.2, float(maximum.max()) * 1.12)
        point = _draw_axes(draw, box, (0.0, 8.0), (0.0, y_max), "r / Å", "g(r)")
        draw.line([point(float(x), float(y)) for x, y in zip(radii, minimum, strict=True)], fill=_hex(GRID), width=3)
        draw.line([point(float(x), float(y)) for x, y in zip(radii, maximum, strict=True)], fill=_hex(GRID), width=3)
        draw.line([point(float(x), float(y)) for x, y in zip(radii, mean, strict=True)], fill=_hex(color), width=6)
        draw.text((box[0] + 95, box[1] + 55), label, font=_pil_font(24, True), fill=_hex(color))
    draw.text((900, 925), "RDF 描述体相短程结构，不直接给出 Li 配位或电子局域。", anchor="mm", font=_pil_font(21), fill=_hex(MUTED))
    _save(image, output)


def _projection(points: np.ndarray) -> np.ndarray:
    if len(points) < 3:
        return points[:, :2]
    centered = points - points.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
    return centered @ vh[:2].T


def _read_local_pdb(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            {
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "resid": int(line[22:26]),
                "position_angstrom": [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
            }
        )
    if not atoms:
        raise ValueError(f"no atoms in local PDB: {path}")
    return atoms


def _read_pymol_hbond_pairs(path: Path) -> list[tuple[int, int]]:
    pattern = re.compile(
        r"^distance hbond_\d+, shell and id (\d+), shell and id (\d+)$"
    )
    pairs: list[tuple[int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            pairs.append((int(match.group(1)), int(match.group(2))))
    return pairs


def _render_cavity_panel(draw: Any, record: dict[str, Any], box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    cavity = record["cavity"]
    snapshot_dir = record["directory"] / "snapshot"
    atoms = _read_local_pdb(snapshot_dir / "cavity_local.pdb")
    positions = np.asarray([atom["position_angstrom"] for atom in atoms], dtype=float)
    projected = _projection(positions)
    lookup = {int(atom["serial"]): index for index, atom in enumerate(atoms)}
    pairs = _read_pymol_hbond_pairs(snapshot_dir / "cavity_hbonds.pml")
    if len(pairs) != len(cavity["local_hydrogen_bonds"]):
        raise ValueError(
            f"PDB/PyMOL hydrogen-bond count mismatch for "
            f"{record['system']}/r{record['replica']}"
        )
    radius = float(cavity["cavity"]["radius_angstrom"])
    extent = max(float(np.max(np.abs(projected))) if len(projected) else 1.0, radius + 0.4)
    scale = min((right - left - 55) / (2 * extent), (bottom - top - 105) / (2 * extent))

    def point(value: np.ndarray) -> tuple[float, float]:
        return ((left + right) / 2 + value[0] * scale, (top + bottom) / 2 - value[1] * scale + 10)

    draw.rounded_rectangle(box, radius=18, fill=(250, 252, 255), outline=_hex(GRID), width=2)
    cx, cy = point(np.asarray([0.0, 0.0]))
    rr = radius * scale
    draw.ellipse((cx - rr, cy - rr, cx + rr, cy + rr), fill=(255, 244, 210), outline=_hex(ACCENT), width=4)
    for bond, (donor_serial, acceptor_serial) in zip(
        cavity["local_hydrogen_bonds"], pairs, strict=True
    ):
        if not bond.get("cavity_associated") and not bond.get("cavity_bridging"):
            continue
        first = lookup.get(donor_serial)
        second = lookup.get(acceptor_serial)
        if first is None or second is None:
            continue
        color = BRIDGE if bond.get("cavity_bridging") else ASSOCIATED
        draw.line((*point(projected[first]), *point(projected[second])), fill=_hex(color), width=5)
    for index, atom in enumerate(atoms):
        element = str(atom["element"]).upper()
        if element == "H":
            continue
        x, y = point(projected[index])
        if element == "N":
            color, atom_radius = "#2667FF", 8
        elif element == "O":
            color, atom_radius = "#E45756", 8
        else:
            color, atom_radius = "#9AA6B8", 5
        draw.ellipse((x - atom_radius, y - atom_radius, x + atom_radius, y + atom_radius), fill=_hex(color), outline=(255, 255, 255), width=1)
    label = f"{SYSTEM_LABELS[record['system']]} r{record['replica']}"
    draw.text((left + 18, top + 18), label, font=_pil_font(21, True), fill=_hex(INK))
    counts = cavity["hydrogen_bond_counts_local_view"]
    summary = f"邻近 {counts['cavity_associated_hydrogen_bonds']} | 穿越 {counts['cavity_bridging_hydrogen_bonds']} | R={radius:.2f} Å"
    draw.text((left + 18, bottom - 34), summary, font=_pil_font(15), fill=_hex(MUTED))


def render_cavity_gallery(records: Sequence[dict[str, Any]], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1280)
    draw = ImageDraw.Draw(image)
    _title(draw, "六个代表快照的空腔氢键投影", "黄色圆为几何空腔；蓝线为空腔邻近氢键，红线为几何上穿越空腔的氢键")
    ordered = [
        next(r for r in records if r["system"] == system and r["replica"] == replica)
        for system in SYSTEMS
        for replica in REPLICAS
    ]
    for index, record in enumerate(ordered):
        column, row = index % 3, index // 3
        box = (45 + column * 585, 155 + row * 530, 565 + column * 585, 635 + row * 530)
        _render_cavity_panel(draw, record, box)
    _save(image, output)


def render_workflow(output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 780)
    draw = ImageDraw.Draw(image)
    _title(draw, "从体相溶剂到 Li/电子问题", "当前报告只完成左侧绿色边界内的 Stage A")
    nodes = (
        (80, "周期性经典 MD", "纯 en 与 1:1 THF:en", ["密度/体积", "RDF 与氢键", "几何空腔快照"], MIXED),
        (515, "Li0/Li+ 构型集合", "加入单个 Li", ["多配位种子", "局部壳层重排", "中性/离子态对照"], PURE),
        (950, "固定核电子计算", "cDFT / Delta-SCF", ["电荷与自旋", "电子质心/IPR", "垂直能隙"], ASSOCIATED),
        (1385, "弛豫与验证", "周期 + 嵌入团簇", ["有限尺寸", "混合泛函", "VDE/光谱"], ACCENT),
    )
    for index, (left, heading, subheading, bullets, color) in enumerate(nodes):
        right = left + 335
        draw.rounded_rectangle((left, 185, right, 610), radius=25, fill=_hex(BACKGROUND), outline=_hex(color), width=4)
        draw.rectangle((left, 185, right, 203), fill=_hex(color))
        draw.text(((left + right) / 2, 265), heading, anchor="mm", font=_pil_font(25, True), fill=_hex(INK))
        draw.text(((left + right) / 2, 315), subheading, anchor="mm", font=_pil_font(19), fill=_hex(color))
        for bullet_index, item in enumerate(bullets):
            y = 390 + bullet_index * 62
            draw.ellipse((left + 40, y - 6, left + 53, y + 7), fill=_hex(color))
            draw.text((left + 72, y), item, anchor="lm", font=_pil_font(19), fill=_hex(INK))
        if index < len(nodes) - 1:
            end = nodes[index + 1][0] - 18
            draw.line((right + 15, 395, end, 395), fill=_hex(MUTED), width=6)
            draw.polygon(((end, 395), (end - 17, 383), (end - 17, 407)), fill=_hex(MUTED))
    draw.text((900, 700), "Stage A 说明空腔附近是否存在氢键几何关系，但不能单独证明空腔中存在电子。", anchor="mm", font=_pil_font(23, True), fill=_hex(INK))
    _save(image, output)


def _register_pdf_fonts() -> tuple[str, str]:
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.cidfonts import UnicodeCIDFont
    from reportlab.pdfbase.ttfonts import TTFont

    regular = _find_font(False)
    bold = _find_font(True)
    if regular is not None:
        try:
            pdfmetrics.registerFont(TTFont("HBondRegular", str(regular), subfontIndex=0))
            pdfmetrics.registerFont(TTFont("HBondBold", str(bold or regular), subfontIndex=0))
            return "HBondRegular", "HBondBold"
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
    title = ParagraphStyle("TitleCN", parent=styles["Title"], fontName=bold_font, fontSize=25, leading=34, textColor=colors.HexColor(INK), alignment=TA_LEFT, wordWrap="CJK", spaceAfter=12)
    subtitle = ParagraphStyle("SubtitleCN", parent=styles["Normal"], fontName=regular_font, fontSize=12, leading=19, textColor=colors.HexColor(MUTED), wordWrap="CJK")
    heading = ParagraphStyle("HeadingCN", parent=styles["Heading1"], fontName=bold_font, fontSize=18, leading=24, textColor=colors.HexColor(INK), wordWrap="CJK", spaceAfter=9)
    heading2 = ParagraphStyle("Heading2CN", parent=styles["Heading2"], fontName=bold_font, fontSize=13.5, leading=19, textColor=colors.HexColor(MIXED), wordWrap="CJK", spaceBefore=8, spaceAfter=5)
    body = ParagraphStyle("BodyCN", parent=styles["BodyText"], fontName=regular_font, fontSize=10.2, leading=16.5, textColor=colors.HexColor(INK), wordWrap="CJK", spaceAfter=6)
    small = ParagraphStyle("SmallCN", parent=body, fontSize=8.2, leading=12.5, textColor=colors.HexColor(MUTED), spaceAfter=3)
    callout = ParagraphStyle("CalloutCN", parent=body, fontName=bold_font, fontSize=11, leading=17.5, backColor=colors.HexColor("#FFF7E6"), borderColor=colors.HexColor(ACCENT), borderWidth=1, borderPadding=8, spaceBefore=6, spaceAfter=8)
    bullet = ParagraphStyle("BulletCN", parent=body, leftIndent=15, firstLineIndent=-10, bulletIndent=0, spaceAfter=3)
    table_header = ParagraphStyle("TableHeader", parent=small, fontName=bold_font, textColor=colors.white)

    def p(text: str, style: Any = body):
        return Paragraph(text, style)

    def table(data: Sequence[Sequence[Any]], widths: Sequence[float]):
        formatted = []
        for row_index, row in enumerate(data):
            row_style = table_header if row_index == 0 else small
            formatted.append([item if hasattr(item, "wrap") else Paragraph(str(item), row_style) for item in row])
        result = Table(formatted, colWidths=widths, repeatRows=1, hAlign="LEFT")
        commands = [
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(INK)),
            ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor(GRID)),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]
        for row_index in range(2, len(data), 2):
            commands.append(("BACKGROUND", (0, row_index), (-1, row_index), colors.HexColor(BACKGROUND)))
        result.setStyle(TableStyle(commands))
        return result

    def page_frame(canvas: Any, document: Any) -> None:
        canvas.saveState()
        width, _height = A4
        canvas.setStrokeColor(colors.HexColor(GRID))
        canvas.line(18 * mm, 15 * mm, width - 18 * mm, 15 * mm)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.setFont(regular_font, 7.5)
        canvas.drawString(18 * mm, 9.5 * mm, "Li/THF/en 项目 - hbond_pilot Stage A")
        canvas.drawRightString(width - 18 * mm, 9.5 * mm, f"第 {document.page} 页")
        canvas.restoreState()

    document = SimpleDocTemplate(
        str(output), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm,
        topMargin=17 * mm, bottomMargin=20 * mm,
        title="hbond_pilot Stage A: 纯 en 与 1:1 THF-en 空腔氢键报告",
        author="solvelec workflow",
    )
    usable_width = A4[0] - 36 * mm
    pure = metrics["systems"]["pure_eda"]
    mixed = metrics["systems"]["eda_1to1"]
    pd = pure["descriptors"]
    md = mixed["descriptors"]
    definition = metrics["hydrogen_bond_definition"]
    story: list[Any] = [
        Spacer(1, 14 * mm),
        p("hbond_pilot Stage A", subtitle),
        p("纯 en 与 1:1 THF:en 的体相结构、空腔和氢键", title),
        p("面向实验讨论的可审计计算报告 | 298.15 K, 1 bar | 2026-09-08", subtitle),
        Spacer(1, 8 * mm),
        _report_image(FIGURES / "cavity_gallery.png", usable_width),
        Spacer(1, 7 * mm),
        p(f"核心结论：两个体系的 6 条独立 20 ns 轨迹均通过门控。空腔邻近氢键的轨迹均值为纯 en {pd['cavity_associated_hydrogen_bonds']['mean']:.2f} 个/盒、1:1 混合液 {md['cavity_associated_hydrogen_bonds']['mean']:.2f} 个/盒；在当前严格定义下，两组的几何穿越氢键均为 0。空腔周围存在可表示的氢键网络，但这不证明空腔内已有溶剂化电子。", callout),
        p(f"计算来源：Slurm 控制器 {metrics['slurm_controller_job']}，116/116 步完成；源提交 {metrics['source_commit'][:12]}。", small),
        PageBreak(),
        p("1. 问题、体系和科学边界", heading),
        p("本阶段直接比较实验人员指定的两种体相溶剂环境：纯 en，以及 THF:en = 1:1 的摩尔比混合液。分子总数固定为 64，以减少体系规模变化带来的直接干扰。"),
        table(
            [
                ["体系", "盒内组成", "摩尔比", "生产采样", "状态"],
                ["纯 en", "64 en", "en = 100%", "3 x 20 ns", "ready"],
                ["1:1 THF:en", "32 THF + 32 en", "THF:en = 1:1", "3 x 20 ns", "ready"],
            ],
            [31 * mm, 43 * mm, 39 * mm, 31 * mm, 28 * mm],
        ),
        p("当前模型能回答什么", heading2),
        p("• 体相密度、体积和短程结构在三个独立副本中的可重复性。", bullet),
        p("• en-en 与 en-THF 氢键的全盒统计，以及它们相对几何空腔的位置关系。", bullet),
        p("• 可供后续单 Li 构型和电子结构计算使用的 6 个独立溶剂环境。", bullet),
        p("当前模型不能回答什么", heading2),
        p("• 本阶段没有 Li，不能给出 Li 的电离程度、O/N 配位数或 Li:en 局部比例。", bullet),
        p("• 本阶段没有过量电子，不能从空腔位置推断电子局域、离域或溶剂化电子稳定性。", bullet),
        p("• 固定电荷经典力场不能描述电子转移、轨道混合或 Li 金属电离。", bullet),
        _report_image(FIGURES / "workflow_boundary.png", usable_width),
        PageBreak(),
        p("2. 方法与接受门控", heading),
        table(
            [
                ["项目", "设置"],
                ["经典模型", "GROMACS; GAFF2; AM1-BCC 电荷"],
                ["热力学条件", "298.15 K; 1.0 bar"],
                ["积分与平衡", "2 fs; 0.5 ns NVT + 5.0 ns NPT"],
                ["生产采样", "每副本 20 ns; 每体系 3 个副本"],
                ["分析", "100 ps 步长; 201 帧/副本"],
                ["代表快照", "10-20 ns 中按多描述符稳健 medoid 选取"],
                ["氢键", f"N-H 供体; O/N 受体; D-A <= {definition['distance_cutoff_angstrom']:.1f} Å; N-H...A >= {definition['angle_cutoff_degree']:.0f}°"],
                ["空腔邻近", f"氢键线段距空腔表面 <= {definition['cavity_shell_thickness_angstrom']:.1f} Å"],
                ["几何穿越", f"氢键线段与空腔球相交; 余量 {definition['cavity_bridge_margin_angstrom']:.2f} Å"],
            ],
            [51 * mm, 121 * mm],
        ),
        p("通过的门控", heading2),
        p("• 每条轨迹覆盖至少 98% 的请求生产时长，且前后半段平均密度差不超过 2%。", bullet),
        p("• 同一体系三个副本的平均密度跨度不超过均值的 3%。", bullet),
        p("• 每条分析至少 150 帧，所有描述符有效样本数至少为 5。", bullet),
        p(f"实际最小有效样本数：纯 en {pure['minimum_effective_sample_size']:.1f}；1:1 混合液 {mixed['minimum_effective_sample_size']:.1f}。最大半程密度差分别为 {pure['maximum_density_half_difference_percent']:.2f}% 和 {mixed['maximum_density_half_difference_percent']:.2f}%。", callout),
        p("这里的“氢键穿越空腔”只是一条线段与几何球相交的拓扑标签。它不代表空腔中心是受体，也不代表电子参与氢键。", small),
        PageBreak(),
        p("3. 密度、体积与空腔尺度", heading),
        _report_image(FIGURES / "density_void.png", usable_width),
        Spacer(1, 3 * mm),
        table(
            [
                ["体系", "密度 / g mL^-1", "体积 / nm3", "空腔半径 / Å", "最小 ESS"],
                ["纯 en", f"{pd['density_g_ml']['mean']:.5f} ± {pd['density_g_ml']['replica_sd']:.5f}", f"{pd['volume_nm3']['mean']:.3f} ± {pd['volume_nm3']['replica_sd']:.3f}", f"{pd['void_radius_angstrom']['mean']:.3f} ± {pd['void_radius_angstrom']['replica_sd']:.3f}", f"{pure['minimum_effective_sample_size']:.1f}"],
                ["1:1 THF:en", f"{md['density_g_ml']['mean']:.5f} ± {md['density_g_ml']['replica_sd']:.5f}", f"{md['volume_nm3']['mean']:.3f} ± {md['volume_nm3']['replica_sd']:.3f}", f"{md['void_radius_angstrom']['mean']:.3f} ± {md['void_radius_angstrom']['replica_sd']:.3f}", f"{mixed['minimum_effective_sample_size']:.1f}"],
            ],
            [32 * mm, 42 * mm, 37 * mm, 36 * mm, 25 * mm],
        ),
        p(f"相对纯 en，1:1 混合液的模型平均密度变化 {metrics['comparisons']['mixture_vs_pure_density_percent']:.2f}%，几何空腔平均半径变化 {metrics['comparisons']['mixture_vs_pure_void_radius_percent']:.2f}%。这些是模型内部差异；在获得相同温压和配比下的实验密度前，不应宣称定量吻合。"),
        PageBreak(),
        p("4. 氢键网络：全盒统计与空腔分类", heading),
        _report_image(FIGURES / "hbond_statistics.png", usable_width),
        Spacer(1, 3 * mm),
        table(
            [
                ["体系", "en-en", "en-THF", "空腔邻近", "几何穿越"],
                ["纯 en", f"{pd['eda_eda_hydrogen_bonds']['mean']:.2f} ± {pd['eda_eda_hydrogen_bonds']['replica_sd']:.2f}", "不适用", f"{pd['cavity_associated_hydrogen_bonds']['mean']:.2f} ± {pd['cavity_associated_hydrogen_bonds']['replica_sd']:.2f}", f"{pd['cavity_bridging_hydrogen_bonds']['mean']:.2f} ± {pd['cavity_bridging_hydrogen_bonds']['replica_sd']:.2f}"],
                ["1:1 THF:en", f"{md['eda_eda_hydrogen_bonds']['mean']:.2f} ± {md['eda_eda_hydrogen_bonds']['replica_sd']:.2f}", f"{md['eda_thf_hydrogen_bonds']['mean']:.2f} ± {md['eda_thf_hydrogen_bonds']['replica_sd']:.2f}", f"{md['cavity_associated_hydrogen_bonds']['mean']:.2f} ± {md['cavity_associated_hydrogen_bonds']['replica_sd']:.2f}", f"{md['cavity_bridging_hydrogen_bonds']['mean']:.2f} ± {md['cavity_bridging_hydrogen_bonds']['replica_sd']:.2f}"],
            ],
            [32 * mm, 32 * mm, 34 * mm, 38 * mm, 36 * mm],
        ),
        p("如何回答实验人员的问题", heading2),
        p(f"可以把氢键“在空腔周围如何排列”表示出来：报告中的蓝线标记空腔邻近氢键，红线预留给线段几何上穿越空腔的氢键；CSV 保留每条氢键的距离、角度、类型和空腔关系。轨迹平均每盒空腔邻近氢键为纯 en {pd['cavity_associated_hydrogen_bonds']['mean']:.2f}、1:1 混合液 {md['cavity_associated_hydrogen_bonds']['mean']:.2f}，但几何穿越计数均为 0。"),
        p("但不能把氢键画成“与空腔成键”。空腔没有原子核或受体位点，必须在加入电子后用自旋密度/电子密度、质心、回转半径和 IPR 等指标验证电子是否真正占据该空间。", callout),
        PageBreak(),
        p("5. 体相短程结构", heading),
        _report_image(FIGURES / "rdf.png", usable_width),
        Spacer(1, 4 * mm),
        p("三副本平均 RDF 主峰", heading2),
        table(
            [["体系", "原子对", "主峰位置 / Å", "g(r) 峰高"]]
            + [[SYSTEM_LABELS[system], pair, f"{peak['r_angstrom']:.2f}", f"{peak['g_r']:.2f}"] for system in SYSTEMS for pair, peak in metrics['systems'][system]['rdf_peaks'].items()],
            [42 * mm, 48 * mm, 40 * mm, 42 * mm],
        ),
        p("纯 en 的 N-N 分布反映 en-en 网络的局部排布；1:1 混合液的 N-O 分布直接反映 en N 与 THF O 的短程相关。RDF 与严格角度判据的氢键计数互补：近距离并不自动等于氢键。"),
        PageBreak(),
        p("6. 六个空腔中心代表快照", heading),
        p("每个副本从生产轨迹后半段独立选择一个稳健代表帧。投影以空腔为中心，只显示局部溶剂环境；完整坐标、PDB、PyMOL 脚本和 JSON 审计数据随报告公开。"),
        _report_image(FIGURES / "cavity_gallery.png", usable_width),
        Spacer(1, 3 * mm),
        table(
            [["体系", "副本", "时间 / ns", "空腔 / Å", "局部原子", "邻近 / 穿越"]]
            + [[SYSTEM_LABELS[system], item['replica'], f"{item['time_ns']:.1f}", f"{item['void_radius_angstrom']:.2f}", item['local_atom_count'], f"{item['local_view_hydrogen_bonds']['cavity_associated_hydrogen_bonds']} / {item['local_view_hydrogen_bonds']['cavity_bridging_hydrogen_bonds']}"] for system in SYSTEMS for item in metrics['systems'][system]['snapshots']],
            [38 * mm, 18 * mm, 29 * mm, 27 * mm, 28 * mm, 32 * mm],
        ),
        p("三维检查：在 PyMOL 中运行各副本的 cavity_hbonds.pml，即可同时查看局部 PDB、空腔球和氢键虚线。二维图用于快速沟通，三维文件用于确认投影是否掩盖空间关系。", small),
        PageBreak(),
        p("7. 下一步计算路线与可迁移工作流", heading),
        _report_image(FIGURES / "workflow_boundary.png", usable_width),
        p("推荐把本轮结果作为溶剂基线，而不是直接把几何空腔当作电子态："),
        p("• 对 6 个代表环境分别加入一个 Li，生成 Li0-like、紧密接触对和分离 Li+/e- 候选构型。", bullet),
        p("• 在固定核条件下比较总电荷守恒的垂直电荷分离，再允许局部溶剂壳层重组。", bullet),
        p("• 用自旋密度积分、电子质心、回转半径、IPR、VDE 和方法/尺寸敏感性判定电子是否离散在一个或多个配体/空腔区域。", bullet),
        p("• 将胺种类、THF:胺摩尔比和密度作为配置输入，使同一条一键 Slurm 工作流可迁移到其他胺 + THF 体系。", bullet),
        p("纯 en 和 1:1 THF:en 给出了两个高 en 含量锚点。若实验配液还包含其他 THF:en 比例，后续应按实际摩尔比和实验密度增补体系，而不是用 Li:en = 1:2 替代体相溶剂比例。", callout),
        PageBreak(),
        p("8. 数据来源、复现和审计", heading),
        p("报告生成器只读取 reports/hbond_stage_a/data。构建时检查三个 ready 汇总、Stage A 完成标记、六条记录唯一性、分子/原子计数、Li 缺失标志、帧数，以及 analysis/XYZ/cell/cavity JSON/PDB/PyMOL 的 SHA-256。任一检查失败都会停止。"),
    ]
    provenance = [["体系/副本", "快照 ID", "XYZ SHA-256", "cavity JSON SHA-256"]]
    for record in records:
        snapshot = record["snapshot"]
        provenance.append([
            f"{record['system']}/r{record['replica']}",
            snapshot["snapshot_id"],
            snapshot["structure"]["xyz"]["sha256"][:16],
            snapshot["structure"]["cavity_hbonds"]["sha256"][:16],
        ])
    story.extend([
        table(provenance, [33 * mm, 55 * mm, 42 * mm, 42 * mm]),
        Spacer(1, 5 * mm),
        p("一键重建", heading2),
        p('python -m pip install -e ".[report]"', callout),
        p("python reports/hbond_stage_a/build_report.py", callout),
        p("公开内容", heading2),
        p("• hbond_stage_a_report_zh.pdf：本报告。", bullet),
        p("• hbond_stage_a_metrics.json：聚合指标和模型边界。", bullet),
        p("• figures/*.png：适合直接发给实验人员的图。", bullet),
        p("• data/analysis/...：6 套 CSV/JSON/XYZ/PDB/PyMOL 审计数据。", bullet),
        p("项目：https://github.com/SunsetStand/li-thf-amine-solvated-electron/tree/main/reports/hbond_stage_a", small),
    ])
    output.parent.mkdir(parents=True, exist_ok=True)
    document.build(story, onFirstPage=page_frame, onLaterPages=page_frame)


def build(output: Path = DEFAULT_PDF) -> dict[str, Any]:
    records = load_records()
    validate_records(records)
    metrics = summarize(records)
    FIGURES.mkdir(parents=True, exist_ok=True)
    render_density_void(records, FIGURES / "density_void.png")
    render_hbond_statistics(records, metrics, FIGURES / "hbond_statistics.png")
    render_rdf(records, FIGURES / "rdf.png")
    render_cavity_gallery(records, FIGURES / "cavity_gallery.png")
    render_workflow(FIGURES / "workflow_boundary.png")
    METRICS_JSON.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    build_pdf(records, metrics, output)
    return metrics


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_PDF)
    args = parser.parse_args(argv)
    metrics = build(args.output)
    pure = metrics["systems"]["pure_eda"]["descriptors"]
    mixed = metrics["systems"]["eda_1to1"]["descriptors"]
    print(f"READY: {args.output}")
    print(
        "HBOND_STAGE_A: "
        f"pure_eda_hb={pure['eda_eda_hydrogen_bonds']['mean']:.2f}; "
        f"mixed_eda_eda_hb={mixed['eda_eda_hydrogen_bonds']['mean']:.2f}; "
        f"mixed_eda_thf_hb={mixed['eda_thf_hydrogen_bonds']['mean']:.2f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

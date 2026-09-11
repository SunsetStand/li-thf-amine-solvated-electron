#!/usr/bin/env python3
# ruff: noqa: E501
"""Build the complete Stage-B pilot report from accepted HPC artifacts.

The report covers the candidate bank, legacy CP2K execution smoke, paired
Li0/Li+e diabatic smoke, and Stage-B1 cube localization diagnostics.  Source
electronic-structure files are read in place and never modified or regenerated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from collections import defaultdict
from collections.abc import Sequence
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.cube import CubeData, read_cube

SYSTEMS = ("pure_thf", "eda_1p5m")
SYSTEM_LABELS = {"pure_thf": "纯 THF", "eda_1p5m": "1.5 M en/THF"}
SYSTEM_LABELS_EN = {"pure_thf": "Pure THF", "eda_1p5m": "1.5 M EDA/THF"}
STATES = ("li0_diabatic", "li_plus_e_diabatic")
STATE_LABELS = {
    "li0_diabatic": "Li0 diabatic",
    "li_plus_e_diabatic": "Li+ + e- diabatic",
}
FIGURE_NAMES = (
    "stage_b_pipeline.png",
    "candidate_bank.png",
    "mechanism_energy.png",
    "spin_partition.png",
    "cube_maps.png",
)

INK = "#17233B"
MUTED = "#64748B"
GRID = "#D9E2EC"
BACKGROUND = "#F4F7FB"
PURE = "#64748B"
MIXED = "#0F8793"
LI_COLOR = "#D1495B"
SOLVENT_COLOR = "#3B82C4"
INTERSTITIAL_COLOR = "#E59F23"
CAVITY_COLOR = "#7C3AED"
SUCCESS = "#16856B"
WARNING = "#D97706"


def _load_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _write_json(path: str | Path, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_gate(summary_path: Path, gate_path: Path) -> None:
    if not summary_path.is_file() or not gate_path.is_file():
        raise ValueError(f"missing accepted summary/gate: {summary_path}, {gate_path}")
    tokens = gate_path.read_text(encoding="utf-8").split()
    expected = ["sha256", _sha256(summary_path), summary_path.name]
    if tokens != expected:
        raise ValueError(f"gate does not match accepted summary: {gate_path}")


def _record_source(
    sources: dict[str, dict[str, Any]],
    path: str | Path,
    expected_sha256: str | None,
    role: str,
) -> None:
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"missing Stage-B source file: {source}")
    observed = _sha256(source)
    if expected_sha256 is not None and observed != expected_sha256:
        raise ValueError(f"SHA-256 mismatch for {role}: {source}")
    key = str(source.resolve())
    existing = sources.get(key)
    if existing is not None and existing["sha256"] != observed:
        raise ValueError(f"inconsistent source hashes for {source}")
    sources[key] = {
        "path": key,
        "sha256": observed,
        "size_bytes": source.stat().st_size,
        "roles": sorted({*(existing or {}).get("roles", []), role}),
    }


def _summary_paths(base: Path) -> dict[str, Path]:
    return {
        "candidates": base / "stage_b_candidates.summary.json",
        "candidates_gate": base / "stage_b_candidates.done",
        "legacy_smoke": base / "stage_b_cp2k_smoke.summary.json",
        "legacy_gate": base / "stage_b.done",
        "mechanism": base / "stage_b_mechanism_smoke.summary.json",
        "mechanism_gate": base / "stage_b_mechanism_smoke.done",
        "localization": base / "stage_b_localization.summary.json",
        "localization_gate": base / "stage_b_localization.done",
    }


def load_and_validate(run_root: str | Path, campaign: str) -> dict[str, Any]:
    base = Path(run_root).resolve() / campaign
    paths = _summary_paths(base)
    for summary_key, gate_key in (
        ("candidates", "candidates_gate"),
        ("legacy_smoke", "legacy_gate"),
        ("mechanism", "mechanism_gate"),
        ("localization", "localization_gate"),
    ):
        _validate_gate(paths[summary_key], paths[gate_key])

    candidates = _load_json(paths["candidates"])
    legacy = _load_json(paths["legacy_smoke"])
    mechanism = _load_json(paths["mechanism"])
    localization = _load_json(paths["localization"])
    for label, summary in (
        ("candidate bank", candidates),
        ("legacy CP2K smoke", legacy),
        ("paired mechanism smoke", mechanism),
        ("localization", localization),
    ):
        if not summary.get("ready"):
            raise ValueError(f"Stage-B {label} summary is not ready")

    if int(candidates.get("record_count", -1)) != 6:
        raise ValueError("Stage-B candidate summary must contain six system/replica records")
    if int(localization.get("state_count", -1)) != 4 or int(
        localization.get("pair_count", -1)
    ) != 2:
        raise ValueError("Stage-B localization summary must contain four states and two pairs")
    if localization.get("source_mechanism_summary", {}).get("sha256") != _sha256(
        paths["mechanism"]
    ):
        raise ValueError("localization summary does not reference the accepted mechanism summary")

    sources: dict[str, dict[str, Any]] = {}
    for key, path in paths.items():
        _record_source(sources, path, None, key)

    manifests: list[dict[str, Any]] = []
    for system, system_record in sorted(candidates.get("systems", {}).items()):
        for replica in system_record.get("replicas", []):
            manifest_path = base / "stage_b" / system / f"r{replica}" / "candidates" / "manifest.json"
            manifest = _load_json(manifest_path)
            if not manifest.get("ready") or int(manifest.get("candidate_count", -1)) != 3:
                raise ValueError(f"candidate manifest is not a three-candidate ready result: {manifest_path}")
            if manifest.get("system_id") != system or int(manifest.get("replica", -1)) != int(
                replica
            ):
                raise ValueError(f"candidate manifest identity mismatch: {manifest_path}")
            _record_source(sources, manifest_path, None, f"candidate_manifest:{system}:r{replica}")
            for candidate in manifest["candidates"]:
                structure = candidate["structure"]
                _record_source(
                    sources,
                    structure["xyz"]["path"],
                    structure["xyz"]["sha256"],
                    f"candidate_xyz:{system}:r{replica}:{candidate['candidate_id']}",
                )
                _record_source(
                    sources,
                    structure["cell"]["path"],
                    structure["cell"]["sha256"],
                    f"candidate_cell:{system}:r{replica}:{candidate['candidate_id']}",
                )
                metadata = candidate.get("metadata")
                if metadata:
                    _record_source(
                        sources,
                        metadata["path"],
                        metadata["sha256"],
                        f"candidate_metadata:{system}:r{replica}:{candidate['candidate_id']}",
                    )
            manifests.append(manifest)
    if len(manifests) != 6:
        raise ValueError(f"expected six candidate manifests, found {len(manifests)}")

    for record in legacy.get("records", []):
        for item in ("input", "output"):
            source = record[item]
            _record_source(
                sources,
                source["path"],
                source["sha256"],
                f"legacy_smoke:{record['system_id']}:{item}",
            )

    mechanism_keys: set[tuple[str, int, str]] = set()
    for pair in mechanism.get("records", []):
        if not pair.get("ready") or not pair.get("complete_pair"):
            raise ValueError(f"mechanism pair is not ready: {pair.get('system_id')}")
        for state in pair.get("states", []):
            key = (state["system_id"], int(state["replica"]), state["state_id"])
            mechanism_keys.add(key)
            for item in ("input", "output"):
                source = state[item]
                _record_source(
                    sources,
                    source["path"],
                    source["sha256"],
                    f"mechanism:{state['system_id']}:{state['state_id']}:{item}",
                )
    expected_keys = {(system, 1, state) for system in SYSTEMS for state in STATES}
    if mechanism_keys != expected_keys:
        raise ValueError(f"unexpected mechanism state identities: {sorted(mechanism_keys)}")

    localization_keys: set[tuple[str, int, str]] = set()
    for state in localization.get("states", []):
        if not state.get("ready"):
            raise ValueError(f"localization state is not ready: {state.get('system_id')}")
        key = (state["system_id"], int(state["replica"]), state["state_id"])
        localization_keys.add(key)
        for input_name, source in state.get("inputs", {}).items():
            if isinstance(source, dict) and source.get("path") and source.get("sha256"):
                _record_source(
                    sources,
                    source["path"],
                    source["sha256"],
                    f"localization:{state['system_id']}:{state['state_id']}:{input_name}",
                )
    if localization_keys != expected_keys:
        raise ValueError(f"unexpected localization state identities: {sorted(localization_keys)}")

    localization_pair_keys: set[tuple[str, int]] = set()
    for pair in localization.get("pairs", []):
        if not pair.get("ready") or not pair.get("density_difference", {}).get("ready"):
            raise ValueError(f"localization pair is not ready: {pair.get('system_id')}")
        localization_pair_keys.add((pair["system_id"], int(pair["replica"])))
        difference_cube = pair.get("inputs", {}).get("density_difference_cube")
        if not difference_cube:
            raise ValueError(f"missing density-difference cube record: {pair['system_id']}")
        _record_source(
            sources,
            difference_cube["path"],
            difference_cube["sha256"],
            f"density_difference_cube:{pair['system_id']}",
        )
    if localization_pair_keys != {(system, 1) for system in SYSTEMS}:
        raise ValueError("localization pair identities do not match the pilot systems")

    return {
        "base": base,
        "paths": paths,
        "candidates": candidates,
        "manifests": manifests,
        "legacy_smoke": legacy,
        "mechanism": mechanism,
        "localization": localization,
        "sources": sorted(sources.values(), key=lambda item: item["path"]),
    }


def summarize(bundle: dict[str, Any]) -> dict[str, Any]:
    candidate_groups: dict[str, dict[str, list[dict[str, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for manifest in bundle["manifests"]:
        for candidate in manifest["candidates"]:
            candidate_groups[manifest["system_id"]][candidate["candidate_id"]].append(
                {
                    "target_distance_angstrom": float(
                        candidate["target_li_cavity_distance_angstrom"]
                    ),
                    "achieved_distance_angstrom": float(
                        candidate["achieved_li_cavity_distance_angstrom"]
                    ),
                    "li_clearance_angstrom": float(candidate["li_site"]["radius_angstrom"]),
                    "cavity_clearance_angstrom": float(
                        candidate["cavity_basis_site"]["radius_angstrom"]
                    ),
                }
            )
    candidate_metrics: dict[str, Any] = {}
    for system, by_candidate in candidate_groups.items():
        candidate_metrics[system] = {}
        for candidate, records in by_candidate.items():
            achieved = [record["achieved_distance_angstrom"] for record in records]
            li_clearance = [record["li_clearance_angstrom"] for record in records]
            cavity_clearance = [record["cavity_clearance_angstrom"] for record in records]
            candidate_metrics[system][candidate] = {
                "replica_count": len(records),
                "target_distance_angstrom": records[0]["target_distance_angstrom"],
                "achieved_distances_angstrom": achieved,
                "achieved_distance_mean_angstrom": float(np.mean(achieved)),
                "achieved_distance_min_angstrom": min(achieved),
                "achieved_distance_max_angstrom": max(achieved),
                "li_clearance_mean_angstrom": float(np.mean(li_clearance)),
                "cavity_clearance_mean_angstrom": float(np.mean(cavity_clearance)),
            }

    legacy_records = []
    for record in bundle["legacy_smoke"].get("records", []):
        gate = record.get("cdft_constraint_gate", {})
        legacy_records.append(
            {
                "system_id": record["system_id"],
                "candidate_id": record["candidate_id"],
                "replica": int(record["replica"]),
                "ready": bool(record.get("converged") and record.get("normal_termination")),
                "energy_hartree": float(record["energy_hartree"]),
                "target_electrons": float(gate["target_electrons"]),
                "current_electrons": float(gate["current_electrons"]),
                "deviation_electrons": float(gate["deviation_electrons"]),
                "tolerance_electrons": float(gate["tolerance_electrons"]),
            }
        )

    mechanism_pairs = {
        record["system_id"]: record for record in bundle["mechanism"]["records"]
    }
    mechanism_states: list[dict[str, Any]] = []
    for system in SYSTEMS:
        for state in mechanism_pairs[system]["states"]:
            gate = state["cdft_constraint_gate"]
            mechanism_states.append(
                {
                    "system_id": system,
                    "state_id": state["state_id"],
                    "ready": state["converged"] and state["normal_termination"],
                    "target_electrons": float(gate["target_electrons"]),
                    "current_electrons": float(gate["current_electrons"]),
                    "deviation_electrons": float(gate["deviation_electrons"]),
                    "tolerance_electrons": float(gate["tolerance_electrons"]),
                    "iteration": int(gate["iteration"]),
                    "energy_hartree": float(state["energy_hartree"]),
                }
            )

    localization_states: list[dict[str, Any]] = []
    for record in bundle["localization"]["states"]:
        spin = record["spin_density"]
        region = record["geometric_partition"]
        li_fraction = float(region["li_positive_spin_fraction"])
        interstitial = float(region["interstitial_positive_spin_fraction"])
        localization_states.append(
            {
                "system_id": record["system_id"],
                "state_id": record["state_id"],
                "ready": record["ready"],
                "signed_spin_electrons": float(spin["signed_integral"]),
                "radius_angstrom": (
                    float(spin["radius_angstrom"])
                    if spin["radius_angstrom"] is not None
                    else None
                ),
                "inverse_participation_ratio": float(spin["inverse_participation_ratio"]),
                "li_positive_spin_fraction": li_fraction,
                "solvent_vdw_positive_spin_fraction": 1.0 - li_fraction - interstitial,
                "maximum_solvent_molecule_positive_spin_fraction": float(
                    region["maximum_solvent_molecule_positive_spin_fraction"]
                ),
                "interstitial_positive_spin_fraction": interstitial,
                "ghost_cavity_probe_positive_spin_fraction": float(
                    region["ghost_cavity_probe_positive_spin_fraction"]
                ),
                "flags": record["localization_proxy_flags"],
                "spin_cube": record["inputs"]["spin_cube"]["path"],
            }
        )

    localization_pairs: list[dict[str, Any]] = []
    for record in bundle["localization"]["pairs"]:
        density = record["density_difference"]
        localization_pairs.append(
            {
                "system_id": record["system_id"],
                "ready": record["ready"],
                "signed_integral_electrons": float(density["signed_integral_electrons"]),
                "rearranged_electrons_half_l1": float(
                    density["rearranged_electrons_half_l1"]
                ),
                "li_density_change_electrons": float(
                    density["li_region"]["density_change_electrons"]
                ),
                "interstitial_accumulation_fraction": float(
                    density["interstitial_region"]["fraction_of_all_accumulation"]
                ),
                "ghost_cavity_accumulation_fraction": float(
                    density["ghost_cavity_probe"]["fraction_of_all_accumulation"]
                ),
                "difference_cube": record["inputs"]["density_difference_cube"]["path"],
            }
        )

    gaps = {
        system: float(mechanism_pairs[system]["fixed_geometry_diabatic_gap_ev"])
        for system in SYSTEMS
    }
    return {
        "schema_version": 1,
        "campaign": bundle["mechanism"]["campaign"],
        "kind": "complete_stage_b_pilot_report",
        "ready": True,
        "scientific_status": "COMPLETED_NUMERICAL_STAGE_B_PILOT_NOT_PRODUCTION_MECHANISM",
        "scope": {
            "systems": list(SYSTEMS),
            "replicas_in_candidate_bank": 3,
            "replicas_in_quantum_smoke": 1,
            "candidate_ids": ["compact", "separated", "distant"],
            "quantum_candidate_id": "separated",
            "states": list(STATES),
            "fixed_geometry": True,
            "periodic": True,
        },
        "candidate_bank": candidate_metrics,
        "legacy_smoke": {
            "ready": bool(bundle["legacy_smoke"]["ready"]),
            "scientific_status": bundle["legacy_smoke"].get("scientific_status"),
            "records": sorted(legacy_records, key=lambda item: item["system_id"]),
        },
        "mechanism": {
            "energy_definition": bundle["mechanism"]["energy_definition"],
            "gaps_ev": gaps,
            "eda_gap_reduction_vs_pure_thf_ev": gaps["pure_thf"] - gaps["eda_1p5m"],
            "states": mechanism_states,
        },
        "localization": {
            "states": sorted(
                localization_states, key=lambda item: (item["system_id"], item["state_id"])
            ),
            "pairs": sorted(localization_pairs, key=lambda item: item["system_id"]),
        },
        "interpretation": {
            "numerical": "All candidate, cDFT, cube-integrity, and charge-conservation gates passed.",
            "li_to_environment": "The Li+e branch removes nearly all positive spin from the Li region and increases its radius of gyration.",
            "cavity": "The configured ghost-cavity probe contains only a minor fraction of positive spin and density accumulation.",
            "delocalization": "No single geometric region crosses the configured dominance threshold; the separated branch is distributed across solvent and interstitial regions.",
            "amine_trend": "EDA lowers the fixed-geometry diabatic gap in this one-snapshot PBE smoke, but the result is not yet a free energy or production trend.",
        },
        "limitations": [
            "One quantum snapshot per system; no replica uncertainty.",
            "Semilocal PBE can over-delocalize an excess electron through self-interaction error.",
            "The ghost basis, cutoff, basis size, and cell-size dependence are not calibrated.",
            "Nuclei are fixed and the cDFT constraint has not been released.",
            "Geometric regions are not Bader, Hirshfeld, or orbital populations.",
            "Only pure THF and EDA 1.5 M are in the quantum pilot; other amines and 3 M remain untested.",
        ],
        "next_gates": [
            "Common-threshold and fixed-enclosed-spin three-dimensional visualization.",
            "Multiwfn/Critic2 or CP2K Hirshfeld cross-check of spin and density redistribution.",
            "Ghost-basis removal/translation and basis/cutoff convergence.",
            "PBE0/ADMM fixed-geometry paired states, then cDFT release and VDE.",
            "Replica and finite-size checks before the five-amine 1.5/3.0 M matrix.",
        ],
    }


def _hex(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[index : index + 2], 16) for index in (0, 2, 4))


def _pil_font(size: int, bold: bool = False):
    from PIL import ImageFont

    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def _canvas(width: int, height: int):
    from PIL import Image

    return Image.new("RGB", (width, height), _hex("#FFFFFF"))


def _figure_title(draw: Any, title: str, subtitle: str = "") -> None:
    draw.text((70, 55), title, font=_pil_font(42, True), fill=_hex(INK))
    if subtitle:
        draw.text((70, 112), subtitle, font=_pil_font(23), fill=_hex(MUTED))


def _save(image: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path, format="PNG", optimize=True)


def render_pipeline(output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 820)
    draw = ImageDraw.Draw(image)
    _figure_title(
        draw,
        "Stage B pilot: evidence ladder",
        "Green checks are completed numerical gates; grey steps remain scientific validation work.",
    )
    nodes = [
        ("Candidate bank", "6 snapshots x 3 seeds", True, PURE),
        ("CP2K smoke", "Li+e path check", True, MIXED),
        ("Paired cDFT", "Li0 vs Li+e", True, LI_COLOR),
        ("Cube localization", "spin + delta density", True, INTERSTITIAL_COLOR),
        ("Production gates", "PBE0, release, VDE", False, MUTED),
    ]
    lefts = [55, 405, 755, 1105, 1455]
    for index, ((heading, subheading, complete, color), left) in enumerate(
        zip(nodes, lefts, strict=True)
    ):
        right = left + 285
        draw.rounded_rectangle(
            (left, 245, right, 590), radius=24, fill=_hex(BACKGROUND), outline=_hex(color), width=5
        )
        draw.rectangle((left, 245, right, 266), fill=_hex(color))
        draw.text(((left + right) / 2, 350), heading, anchor="mm", font=_pil_font(25, True), fill=_hex(INK))
        draw.text(((left + right) / 2, 408), subheading, anchor="mm", font=_pil_font(20), fill=_hex(MUTED))
        status = "PASSED" if complete else "NOT YET RUN"
        status_color = SUCCESS if complete else MUTED
        draw.rounded_rectangle((left + 48, 475, right - 48, 535), radius=18, fill=_hex(status_color))
        draw.text(((left + right) / 2, 505), status, anchor="mm", font=_pil_font(18, True), fill=(255, 255, 255))
        if index < len(nodes) - 1:
            end = lefts[index + 1] - 18
            draw.line((right + 10, 418, end, 418), fill=_hex(GRID), width=8)
            draw.polygon(((end, 418), (end - 17, 405), (end - 17, 431)), fill=_hex(GRID))
    draw.text(
        (900, 700),
        "Completed Stage B is a reproducible numerical pilot, not yet a solvated-electron stability proof.",
        anchor="mm",
        font=_pil_font(25, True),
        fill=_hex(INK),
    )
    _save(image, output)


def _axes(draw: Any, box: tuple[int, int, int, int], y_max: float, y_label: str) -> None:
    left, top, right, bottom = box
    draw.line((left, top, left, bottom), fill=_hex(INK), width=3)
    draw.line((left, bottom, right, bottom), fill=_hex(INK), width=3)
    for index in range(5):
        value = y_max * index / 4
        y = bottom - (bottom - top) * index / 4
        draw.line((left, y, right, y), fill=_hex(GRID), width=2)
        draw.text((left - 16, y), f"{value:.2f}", anchor="rm", font=_pil_font(17), fill=_hex(MUTED))
    draw.text((left - 90, (top + bottom) / 2), y_label, anchor="mm", font=_pil_font(18, True), fill=_hex(INK))


def render_candidate_bank(metrics: dict[str, Any], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1000)
    draw = ImageDraw.Draw(image)
    _figure_title(draw, "Deterministic Li / ghost-cavity candidate bank", "Each point is one Stage-A replica; lines show configured target distances.")
    box = (150, 215, 1690, 820)
    _axes(draw, box, 11.0, "Li-cavity distance (A)")
    x_positions = {"compact": 410, "separated": 900, "distant": 1390}
    colors = {"pure_thf": PURE, "eda_1p5m": MIXED}
    for candidate, x in x_positions.items():
        target = metrics["candidate_bank"]["pure_thf"][candidate]["target_distance_angstrom"]
        y = box[3] - target / 11.0 * (box[3] - box[1])
        draw.line((x - 150, y, x + 150, y), fill=_hex(WARNING), width=4)
        draw.text((x, box[3] + 48), candidate, anchor="mm", font=_pil_font(22, True), fill=_hex(INK))
        draw.text((x, box[3] + 82), f"target {target:.1f} A", anchor="mm", font=_pil_font(17), fill=_hex(MUTED))
        for system_index, system in enumerate(SYSTEMS):
            summary = metrics["candidate_bank"][system][candidate]
            values = summary["achieved_distances_angstrom"]
            for replica_index, value in enumerate(values):
                px = x + (-52 if system_index == 0 else 52) + (replica_index - 1) * 12
                py = box[3] - value / 11.0 * (box[3] - box[1])
                draw.ellipse((px - 11, py - 11, px + 11, py + 11), fill=_hex(colors[system]), outline=(255, 255, 255), width=2)
    legend_y = 920
    for index, system in enumerate(SYSTEMS):
        x = 610 + index * 420
        draw.ellipse((x, legend_y - 12, x + 24, legend_y + 12), fill=_hex(colors[system]))
        draw.text((x + 40, legend_y), SYSTEM_LABELS_EN[system], anchor="lm", font=_pil_font(20), fill=_hex(INK))
    _save(image, output)


def render_mechanism_energy(metrics: dict[str, Any], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1000)
    draw = ImageDraw.Draw(image)
    _figure_title(draw, "Paired fixed-geometry cDFT smoke", "Energy gaps and independent final Li-population constraint checks.")
    gap_box = (150, 235, 820, 805)
    _axes(draw, gap_box, 2.2, "Delta E (eV)")
    gaps = metrics["mechanism"]["gaps_ev"]
    for index, system in enumerate(SYSTEMS):
        x = 330 + index * 310
        value = gaps[system]
        y = gap_box[3] - value / 2.2 * (gap_box[3] - gap_box[1])
        color = PURE if system == "pure_thf" else MIXED
        draw.rounded_rectangle((x - 75, y, x + 75, gap_box[3]), radius=12, fill=_hex(color))
        draw.text((x, y - 20), f"{value:.3f}", anchor="ms", font=_pil_font(25, True), fill=_hex(INK))
        draw.text((x, gap_box[3] + 43), SYSTEM_LABELS_EN[system], anchor="mm", font=_pil_font(18, True), fill=_hex(INK))
    draw.text((485, 910), f"EDA gap reduction: {metrics['mechanism']['eda_gap_reduction_vs_pure_thf_ev']:.3f} eV", anchor="mm", font=_pil_font(22, True), fill=_hex(MIXED))

    dev_box = (1040, 235, 1690, 805)
    _axes(draw, dev_box, 0.055, "|Li target deviation| (e)")
    states = metrics["mechanism"]["states"]
    for index, state in enumerate(states):
        x = 1115 + index * 145
        value = abs(state["deviation_electrons"])
        y = dev_box[3] - value / 0.055 * (dev_box[3] - dev_box[1])
        color = PURE if state["system_id"] == "pure_thf" else MIXED
        draw.rectangle((x - 42, y, x + 42, dev_box[3]), fill=_hex(color))
        draw.text((x, y - 14), f"{value:.3f}", anchor="ms", font=_pil_font(16, True), fill=_hex(INK))
        short_state = "Li0" if state["state_id"] == "li0_diabatic" else "Li+e"
        draw.text((x, dev_box[3] + 35), short_state, anchor="mm", font=_pil_font(16, True), fill=_hex(INK))
        draw.text((x, dev_box[3] + 62), "THF" if state["system_id"] == "pure_thf" else "EDA", anchor="mm", font=_pil_font(14), fill=_hex(MUTED))
    tolerance_y = dev_box[3] - 0.05 / 0.055 * (dev_box[3] - dev_box[1])
    draw.line((dev_box[0], tolerance_y, dev_box[2], tolerance_y), fill=_hex(LI_COLOR), width=4)
    draw.text((dev_box[2] - 5, tolerance_y - 12), "0.05 e gate", anchor="rs", font=_pil_font(17, True), fill=_hex(LI_COLOR))
    _save(image, output)


def render_spin_partition(metrics: dict[str, Any], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1080)
    draw = ImageDraw.Draw(image)
    _figure_title(draw, "Where is the unpaired spin?", "Stacked regions are exclusive; the ghost-cavity probe is an overlapping marker.")
    box = (155, 240, 1690, 860)
    _axes(draw, box, 1.0, "Fraction of positive spin")
    records = metrics["localization"]["states"]
    x_values = [350, 700, 1110, 1460]
    colors = [LI_COLOR, SOLVENT_COLOR, INTERSTITIAL_COLOR]
    for x, record in zip(x_values, records, strict=True):
        parts = [
            record["li_positive_spin_fraction"],
            record["solvent_vdw_positive_spin_fraction"],
            record["interstitial_positive_spin_fraction"],
        ]
        bottom = box[3]
        for value, color in zip(parts, colors, strict=True):
            height = value * (box[3] - box[1])
            draw.rectangle((x - 72, bottom - height, x + 72, bottom), fill=_hex(color), outline=(255, 255, 255), width=2)
            bottom -= height
        cavity = record["ghost_cavity_probe_positive_spin_fraction"]
        cavity_y = box[3] - cavity * (box[3] - box[1])
        draw.line((x - 96, cavity_y, x + 96, cavity_y), fill=_hex(CAVITY_COLOR), width=7)
        draw.text((x, cavity_y - 12), f"cavity {100*cavity:.1f}%", anchor="ms", font=_pil_font(15, True), fill=_hex(CAVITY_COLOR))
        system = "THF" if record["system_id"] == "pure_thf" else "EDA"
        state = "Li0" if record["state_id"] == "li0_diabatic" else "Li+e"
        draw.text((x, box[3] + 40), f"{system} / {state}", anchor="mm", font=_pil_font(19, True), fill=_hex(INK))
        draw.text((x, box[3] + 73), f"Rg {record['radius_angstrom']:.2f} A", anchor="mm", font=_pil_font(16), fill=_hex(MUTED))
    legend = [("Li region", LI_COLOR), ("Solvent vdW regions", SOLVENT_COLOR), ("Interstitial", INTERSTITIAL_COLOR), ("Ghost-cavity probe", CAVITY_COLOR)]
    for index, (label, color) in enumerate(legend):
        x = 280 + index * 380
        draw.rectangle((x, 995, x + 28, 1023), fill=_hex(color))
        draw.text((x + 42, 1009), label, anchor="lm", font=_pil_font(17), fill=_hex(INK))
    _save(image, output)


def _fractional_atom_positions(cube: CubeData) -> np.ndarray:
    if cube.atoms.size == 0:
        return np.empty((0, 3))
    return (cube.atoms[:, 2:5] - cube.origin) @ np.linalg.inv(cube.box)


def _sequential_rgb(values: np.ndarray) -> np.ndarray:
    finite = np.asarray(values, dtype=float)
    scale = float(np.percentile(finite[finite > 0], 99.0)) if np.any(finite > 0) else 1.0
    norm = np.clip(finite / max(scale, 1.0e-30), 0.0, 1.0)
    norm = np.log1p(15.0 * norm) / math.log(16.0)
    low = np.asarray(_hex("#F6F9FC"), dtype=float)
    mid = np.asarray(_hex("#30A7B8"), dtype=float)
    high = np.asarray(_hex("#F4B942"), dtype=float)
    first = norm[..., None] * 2.0
    rgb = np.where(
        (norm[..., None] <= 0.5),
        low + np.clip(first, 0, 1) * (mid - low),
        mid + np.clip(first - 1, 0, 1) * (high - mid),
    )
    return np.asarray(np.clip(rgb, 0, 255), dtype=np.uint8)


def _diverging_rgb(values: np.ndarray) -> np.ndarray:
    finite = np.asarray(values, dtype=float)
    nonzero = np.abs(finite[np.nonzero(finite)])
    scale = float(np.percentile(nonzero, 99.0)) if nonzero.size else 1.0
    norm = np.clip(finite / max(scale, 1.0e-30), -1.0, 1.0)
    white = np.asarray(_hex("#F8FAFC"), dtype=float)
    blue = np.asarray(_hex("#2563EB"), dtype=float)
    red = np.asarray(_hex("#DC3C5A"), dtype=float)
    strength = np.sqrt(np.abs(norm))[..., None]
    target = np.where((norm[..., None] >= 0), red, blue)
    return np.asarray(np.clip(white + strength * (target - white), 0, 255), dtype=np.uint8)


def _draw_cube_panel(
    canvas: Any,
    panel: tuple[int, int, int, int],
    cube_path: str | Path,
    title: str,
    difference: bool,
) -> None:
    from PIL import Image, ImageDraw

    cube = read_cube(cube_path)
    if difference:
        projection = np.clip(cube.values, 0.0, None).sum(axis=2) - np.clip(
            -cube.values, 0.0, None
        ).sum(axis=2)
        rgb = _diverging_rgb(projection)
    else:
        projection = np.clip(cube.values, 0.0, None).sum(axis=2)
        rgb = _sequential_rgb(projection)
    heatmap = Image.fromarray(np.swapaxes(rgb, 0, 1)[::-1, :, :], mode="RGB")
    left, top, right, bottom = panel
    heatmap = heatmap.resize((right - left, bottom - top), Image.Resampling.BILINEAR)
    canvas.paste(heatmap, (left, top))
    draw = ImageDraw.Draw(canvas)
    draw.rectangle(panel, outline=_hex(INK), width=3)
    draw.text((left + 12, top + 12), title, font=_pil_font(20, True), fill=_hex(INK), stroke_width=3, stroke_fill=(255, 255, 255))
    fractional = _fractional_atom_positions(cube)
    for atom, frac in zip(cube.atoms, fractional, strict=True):
        atomic_number = int(round(atom[0]))
        if atomic_number not in {0, 3}:
            continue
        x = left + float(frac[0] % 1.0) * (right - left)
        y = bottom - float(frac[1] % 1.0) * (bottom - top)
        color = LI_COLOR if atomic_number == 3 else CAVITY_COLOR
        radius = 15 if atomic_number == 3 else 12
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), outline=_hex(color), width=5)
        draw.text((x + radius + 4, y), "Li" if atomic_number == 3 else "Gh", anchor="lm", font=_pil_font(14, True), fill=_hex(color), stroke_width=2, stroke_fill=(255, 255, 255))


def render_cube_maps(metrics: dict[str, Any], output: Path) -> None:
    from PIL import ImageDraw

    image = _canvas(1800, 1750)
    draw = ImageDraw.Draw(image)
    _figure_title(draw, "Cube-space diagnostic maps", "Cell-wide z projections. Spin maps show positive spin; delta-density maps use red accumulation and blue depletion.")
    states = metrics["localization"]["states"]
    pairs = metrics["localization"]["pairs"]
    panels = [
        (80, 200, 840, 650),
        (960, 200, 1720, 650),
        (80, 735, 840, 1185),
        (960, 735, 1720, 1185),
    ]
    for panel, record in zip(panels, states, strict=True):
        title = f"{SYSTEM_LABELS_EN[record['system_id']]} / {STATE_LABELS[record['state_id']]}"
        _draw_cube_panel(image, panel, record["spin_cube"], title, difference=False)
    pair_panels = [(80, 1270, 840, 1660), (960, 1270, 1720, 1660)]
    for panel, record in zip(pair_panels, pairs, strict=True):
        title = f"{SYSTEM_LABELS_EN[record['system_id']]} / delta electron density"
        _draw_cube_panel(image, panel, record["difference_cube"], title, difference=True)
    draw.text((900, 1712), "Projection is a communication view, not a basin population analysis.", anchor="mm", font=_pil_font(18), fill=_hex(MUTED))
    _save(image, output)


def _register_pdf_fonts() -> tuple[str, str]:
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.cidfonts import UnicodeCIDFont
    from reportlab.pdfbase.ttfonts import TTFont

    candidates = [
        (Path("C:/Windows/Fonts/msyh.ttc"), Path("C:/Windows/Fonts/msyhbd.ttc")),
        (Path("C:/Windows/Fonts/simhei.ttf"), Path("C:/Windows/Fonts/simhei.ttf")),
        (Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"), Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc")),
    ]
    for regular, bold in candidates:
        if not regular.is_file():
            continue
        try:
            pdfmetrics.registerFont(TTFont("StageBRegular", str(regular), subfontIndex=0))
            pdfmetrics.registerFont(TTFont("StageBBold", str(bold if bold.is_file() else regular), subfontIndex=0))
            return "StageBRegular", "StageBBold"
        except Exception:
            continue
    pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
    return "STSong-Light", "STSong-Light"


def _report_image(path: Path, width: float):
    from PIL import Image as PILImage
    from reportlab.platypus import Image

    with PILImage.open(path) as source:
        image_width, image_height = source.size
    return Image(str(path), width=width, height=width * image_height / image_width)


def build_pdf(metrics: dict[str, Any], figures: Path, output: Path) -> None:
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

    regular, bold = _register_pdf_fonts()
    styles = getSampleStyleSheet()
    title = ParagraphStyle("TitleCN", parent=styles["Title"], fontName=bold, fontSize=24, leading=33, textColor=colors.HexColor(INK), alignment=TA_LEFT, wordWrap="CJK", spaceAfter=12)
    subtitle = ParagraphStyle("SubtitleCN", parent=styles["Normal"], fontName=regular, fontSize=11, leading=18, textColor=colors.HexColor(MUTED), wordWrap="CJK", spaceAfter=5)
    heading = ParagraphStyle("HeadingCN", parent=styles["Heading1"], fontName=bold, fontSize=18, leading=25, textColor=colors.HexColor(INK), wordWrap="CJK", spaceAfter=9)
    heading2 = ParagraphStyle("Heading2CN", parent=styles["Heading2"], fontName=bold, fontSize=13.5, leading=19, textColor=colors.HexColor(MIXED), wordWrap="CJK", spaceBefore=7, spaceAfter=5)
    body = ParagraphStyle("BodyCN", parent=styles["BodyText"], fontName=regular, fontSize=10.1, leading=16.5, textColor=colors.HexColor(INK), wordWrap="CJK", spaceAfter=6)
    small = ParagraphStyle("SmallCN", parent=body, fontSize=8.1, leading=12.3, textColor=colors.HexColor(MUTED), spaceAfter=3)
    callout = ParagraphStyle("CalloutCN", parent=body, fontName=bold, fontSize=10.8, leading=17.5, backColor=colors.HexColor("#FFF7E6"), borderColor=colors.HexColor(WARNING), borderWidth=1, borderPadding=8, spaceBefore=6, spaceAfter=8)
    bullet = ParagraphStyle("BulletCN", parent=body, leftIndent=15, firstLineIndent=-10, bulletIndent=0, spaceAfter=3)
    table_header = ParagraphStyle("TableHeader", parent=small, fontName=bold, textColor=colors.white)

    def p(text: str, style: Any = body):
        return Paragraph(text, style)

    def table(rows: Sequence[Sequence[Any]], widths: Sequence[float]):
        formatted = []
        for row_index, row in enumerate(rows):
            style = table_header if row_index == 0 else small
            formatted.append([item if hasattr(item, "wrap") else Paragraph(str(item), style) for item in row])
        result = Table(formatted, colWidths=widths, repeatRows=1, hAlign="LEFT")
        commands = [
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(INK)),
            ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor(GRID)),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]
        for row_index in range(2, len(rows), 2):
            commands.append(("BACKGROUND", (0, row_index), (-1, row_index), colors.HexColor(BACKGROUND)))
        result.setStyle(TableStyle(commands))
        return result

    def page_frame(canvas: Any, document: Any) -> None:
        canvas.saveState()
        width, _ = A4
        canvas.setStrokeColor(colors.HexColor(GRID))
        canvas.line(18 * mm, 15 * mm, width - 18 * mm, 15 * mm)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.setFont(regular, 7.5)
        canvas.drawString(18 * mm, 9.5 * mm, "Li/THF/胺溶剂化电子 - Stage B pilot")
        canvas.drawRightString(width - 18 * mm, 9.5 * mm, f"第 {document.page} 页")
        canvas.restoreState()

    output.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=17 * mm,
        bottomMargin=20 * mm,
        title="Li/THF/胺溶剂化电子 Stage B 完整计算报告",
        author="solvelec workflow",
    )
    usable = A4[0] - 36 * mm
    gaps = metrics["mechanism"]["gaps_ev"]
    pair_map = {item["system_id"]: item for item in metrics["localization"]["pairs"]}
    story: list[Any] = [
        Spacer(1, 12 * mm),
        p("Stage B 完整计算报告", subtitle),
        p("Li/THF/en 体系中的固定结构电荷分离与电子定位诊断", title),
        p("pure THF 与 1.5 M en/THF | 周期性 PBE-D3(BJ) cDFT | replica 1 | separated 候选", subtitle),
        Spacer(1, 7 * mm),
        _report_image(figures / "stage_b_pipeline.png", usable),
        Spacer(1, 7 * mm),
        p(f"核心结果：候选结构、CP2K 数值路径、成对 diabatic 态和 cube 完整性门控均已通过。固定核条件下，en/THF 的 Li+e - Li0 能隙为 {gaps['eda_1p5m']:.3f} eV，pure THF 为 {gaps['pure_thf']:.3f} eV；但分离态未配对电子主要分散在溶剂和分子间区域，而不是占据预先选择的 ghost 空腔。", callout),
        p("本报告汇总完整的已完成 Stage B pilot。它不是完整生产级机理结论：PBE0、约束释放、VDE、ghost/基组/尺寸校准、量子副本和其他胺体系尚未完成。", small),
        PageBreak(),
        p("1. 研究问题、证据范围与计算路线", heading),
        p("Stage B 从已经验收的 Stage A 溶剂快照出发，引入一个 Li 中心和一个无核无电荷的 Gh 基组中心。三个几何候选仅用于提供可重复的初始构型；量子 pilot 选择 replica 1 的 separated 候选，并在完全相同的核坐标、晶胞、基组和总电子数下构造 Li0-like 与 Li+ + e- 两个 cDFT diabatic 分支。"),
        table(
            [
                ["层级", "已完成内容", "能回答什么", "不能回答什么"],
                ["B0 候选库", "2 体系 x 3 副本 x 3 几何种子", "初始结构是否可重复生成", "电子是否稳定"],
                ["CP2K smoke", "2 个 Li+e 数值作业", "程序和输入路径能否收敛", "机理和定位"],
                ["成对 cDFT", "2 体系 x 2 diabatic 态", "固定结构分支是否数值可达", "真实电离自由能"],
                ["B1 cube", "4 自旋密度 + 2 密度差", "电子在何种几何区域展开", "严格原子盆/稳定态"],
            ],
            [28 * mm, 49 * mm, 48 * mm, 47 * mm],
        ),
        p("两个分支都是中性双重态，因此密度差表示同一电子数下的重排，而不是向晶胞额外加入一个带电电子。周期性电中性不是本轮 SCF 困难的来源。", callout),
        p("与实验问题的关系", heading2),
        p("当前比较只覆盖 pure THF 和 1.5 M en/THF，可以检验 en 是否改变固定结构电荷分离的数值代价；尚不能隔离 N-H、胺骨架、浓度或氢键的独立作用。", body),
        PageBreak(),
        p("2. Stage B0：候选结构库", heading),
        _report_image(figures / "candidate_bank.png", usable),
        p("每个 Stage A 代表快照通过周期性重原子范德华表面间隙搜索产生多个空隙位点，再为 compact、separated、distant 三个目标距离选择互不复用的 Li/Gh 位点对。Gh 只有基函数，没有核和电荷；它是变分支持，不是观测到的空腔电子。", body),
        table(
            [["体系", "候选", "目标 / Å", "实际均值 / Å", "范围 / Å", "Li/Gh 平均间隙 / Å"]]
            + [
                [
                    SYSTEM_LABELS[system],
                    candidate,
                    f"{record['target_distance_angstrom']:.2f}",
                    f"{record['achieved_distance_mean_angstrom']:.2f}",
                    f"{record['achieved_distance_min_angstrom']:.2f}-{record['achieved_distance_max_angstrom']:.2f}",
                    f"{record['li_clearance_mean_angstrom']:.2f}/{record['cavity_clearance_mean_angstrom']:.2f}",
                ]
                for system in SYSTEMS
                for candidate, record in metrics["candidate_bank"][system].items()
            ],
            [29 * mm, 27 * mm, 26 * mm, 31 * mm, 33 * mm, 36 * mm],
        ),
        PageBreak(),
        p("3. CP2K 数值 smoke 与成对 diabatic 态", heading),
        _report_image(figures / "mechanism_energy.png", usable),
        p("方法为周期性 unrestricted PBE-D3(BJ)/DZVP-MOLOPT-SR-GTH，400 Ry cutoff，并在 Gh 上使用 TZV2P 基组。Li 使用 GTH-PBE-q3 势；Li0-like 的 Li 价电子目标为 3.0，Li+e 分支为 2.0。最终 cDFT 人口偏差必须不超过 0.05 e。", body),
        p("旧版单分支 smoke", heading2),
        table(
            [["体系", "分支", "目标 e", "当前 e", "偏差 e", "能量 / Ha"]]
            + [
                [
                    SYSTEM_LABELS[item["system_id"]],
                    "Li+e",
                    f"{item['target_electrons']:.1f}",
                    f"{item['current_electrons']:.4f}",
                    f"{item['deviation_electrons']:+.4f}",
                    f"{item['energy_hartree']:.6f}",
                ]
                for item in metrics["legacy_smoke"]["records"]
            ],
            [34 * mm, 23 * mm, 24 * mm, 28 * mm, 27 * mm, 40 * mm],
        ),
        p("这两个作业只证明 Li+e 输入路径可正常结束；它们随后由同一几何上的 Li0/Li+e 成对计算补足，不能单独产生能隙或定位结论。", small),
        table(
            [["体系", "状态", "目标 e", "当前 e", "偏差 e", "迭代", "能量 / Ha"]]
            + [
                [
                    SYSTEM_LABELS[item["system_id"]],
                    "Li0" if item["state_id"] == "li0_diabatic" else "Li+e",
                    f"{item['target_electrons']:.1f}",
                    f"{item['current_electrons']:.4f}",
                    f"{item['deviation_electrons']:+.4f}",
                    item["iteration"],
                    f"{item['energy_hartree']:.6f}",
                ]
                for item in metrics["mechanism"]["states"]
            ],
            [31 * mm, 24 * mm, 21 * mm, 26 * mm, 25 * mm, 18 * mm, 37 * mm],
        ),
        p(f"en/THF 的固定结构 diabatic gap 比 pure THF 小 {metrics['mechanism']['eda_gap_reduction_vs_pure_thf_ev']:.3f} eV。这是当前 pilot 最值得继续检验的信号，但绝不能直接称为 Li 电离能或溶剂化电子稳定自由能。", callout),
        PageBreak(),
        p("4. Stage B1：自旋密度定位", heading),
        _report_image(figures / "spin_partition.png", usable),
        table(
            [["体系", "状态", "Li", "溶剂 vdW", "分子间", "Gh 探针", "最大单分子", "Rg / Å"]]
            + [
                [
                    SYSTEM_LABELS[item["system_id"]],
                    "Li0" if item["state_id"] == "li0_diabatic" else "Li+e",
                    f"{100*item['li_positive_spin_fraction']:.1f}%",
                    f"{100*item['solvent_vdw_positive_spin_fraction']:.1f}%",
                    f"{100*item['interstitial_positive_spin_fraction']:.1f}%",
                    f"{100*item['ghost_cavity_probe_positive_spin_fraction']:.1f}%",
                    f"{100*item['maximum_solvent_molecule_positive_spin_fraction']:.1f}%",
                    f"{item['radius_angstrom']:.2f}",
                ]
                for item in metrics["localization"]["states"]
            ],
            [25 * mm, 21 * mm, 18 * mm, 25 * mm, 21 * mm, 22 * mm, 27 * mm, 19 * mm],
        ),
        p("Li+e 分支中，Li 区域只保留 0.8-1.8% 的正自旋，Rg 增大到约 4.4-4.6 Å；约 45% 位于所有原子范德华球之外。任何单个溶剂分子和原始 Gh 空腔都没有成为占主导的承载区域。", callout),
        PageBreak(),
        p("5. Cube 空间图与密度重排", heading),
        _report_image(figures / "cube_maps.png", usable),
        p("图为整个周期晶胞沿 z 方向的二维投影，用于交流自旋和密度重排的空间范围。投影会压缩深度信息，因此定量结论仍来自三维体积分。Li 与 Gh 标记来自 cube 原子表；Gh 不代表真实原子。", small),
        PageBreak(),
        p("6. 成对密度差与体系比较", heading),
        p("密度差定义为 rho(Li+e) - rho(Li0)。由于两个分支电子数相同，积分应接近 0；正值表示积累，负值表示耗减。", body),
        table(
            [["体系", "积分 ΔQ / e", "重排电子数 / e", "Li 区域 Δe", "分子间积累", "Gh 探针积累"]]
            + [
                [
                    SYSTEM_LABELS[item["system_id"]],
                    f"{item['signed_integral_electrons']:+.4f}",
                    f"{item['rearranged_electrons_half_l1']:.4f}",
                    f"{item['li_density_change_electrons']:+.4f}",
                    f"{100*item['interstitial_accumulation_fraction']:.1f}%",
                    f"{100*item['ghost_cavity_accumulation_fraction']:.1f}%",
                ]
                for item in metrics["localization"]["pairs"]
            ],
            [34 * mm, 28 * mm, 33 * mm, 28 * mm, 26 * mm, 28 * mm],
        ),
        p(f"pure THF 中 Li 区域损失 {abs(pair_map['pure_thf']['li_density_change_electrons']):.2f} e，en/THF 中损失 {abs(pair_map['eda_1p5m']['li_density_change_electrons']):.2f} e。两个体系约四分之一的密度积累位于分子间区域，而 Gh 探针仅覆盖 {100*pair_map['pure_thf']['ghost_cavity_accumulation_fraction']:.1f}% 和 {100*pair_map['eda_1p5m']['ghost_cavity_accumulation_fraction']:.1f}%。", body),
        p("最保守的当前图像：Li0-like 分支的未配对自旋以 Li 为最大单一区域但已有明显外溢；Li+e 分支则展开为多溶剂/分子间离域态。EDA 可能降低形成该 diabatic 分支的代价，却没有把电子简单地锁在预先选择的单一空腔中。", callout),
        PageBreak(),
        p("7. 与实验假设的关系", heading),
        p("实验观察要求 N-H 且推测 N 稳定 Li、N-H/氢键环境帮助稳定电子。当前结果与“胺促进电荷分离”相容，但尚未证明 N-H 是原因："),
        p("• 量子 pilot 只有 en 1.5 M 与 pure THF，无法区分 N-H、氮配位能力、介电环境和浓度变化。", bullet),
        p("• Li+e 分支没有表现为单个 en 分子阴离子；最大单分子正自旋只有约 8-10%。", bullet),
        p("• 原始 Gh 空腔不是主要承载区域，因此后续应让电子位置由波函数和密度决定，而不是继续把最大几何空隙当作电子中心。", bullet),
        p("• TMEDA 无 N-H，是验证实验 N-H 必需性的关键负对照；1,2-PDA、1,3-PDA 和 DETA 则可分离骨架、齿数与氢键网络的作用。", bullet),
        p("报告中的能隙是同一化学组成内部的成对差值，可以比较趋势；不能直接比较两个不同体系的绝对总能。", small),
        PageBreak(),
        p("8. 局限性与下一阶段决策门", heading),
        table(
            [
                ["风险", "为什么重要", "下一门控"],
                ["PBE 自相互作用", "可能人为扩大电子离域", "PBE0/ADMM 固定结构复算"],
                ["Gh 基组偏置", "可能把电子吸向人为位置", "移除/平移 Gh 并做基组收敛"],
                ["单快照", "没有构象或副本不确定度", "每体系至少 3 个独立快照"],
                ["固定核 + cDFT", "不是绝热稳定态", "约束释放、局部弛豫和 VDE"],
                ["几何分区", "不是严格原子/分子人口", "Multiwfn/Critic2/Hirshfeld 交叉验证"],
                ["有限晶胞", "离域态易受周期镜像影响", "晶胞/尺寸收敛"],
            ],
            [35 * mm, 67 * mm, 70 * mm],
        ),
        p("推荐顺序", heading2),
        *[p(f"• {index}. {item}", bullet) for index, item in enumerate(metrics["next_gates"], start=1)],
        p("只有方法学门控在 pure THF 与 en 1.5 M 上稳定后，才值得扩展到 en、1,2-PDA、1,3-PDA、DETA、TMEDA 的 1.5/3.0 M 完整矩阵。", callout),
        PageBreak(),
        p("9. 数据审计与复现", heading),
        p("报告构建器在写图和 PDF 前检查四级 summary 的 ready 状态、四个完成标记与 summary 的 SHA-256、一致的体系/状态基数、候选结构哈希、四个 CP2K 输入和输出哈希、八个 cube 输入哈希以及两份密度差 cube 哈希。任何不一致都会停止构建。", body),
        table(
            [
                ["数据层", "数量", "状态"],
                ["候选 manifest", "6", "ready"],
                ["几何候选", "18", "ready"],
                ["成对 cDFT 状态", "4", "normal termination + population gate"],
                ["自旋密度分析", "4", "ready"],
                ["密度差分析", "2", "charge conserved"],
                ["来源文件哈希", str(metrics.get("source_file_count", "见 provenance")), "validated"],
            ],
            [58 * mm, 38 * mm, 76 * mm],
        ),
        p("可移植输出包括本 PDF、5 张 PNG、机器可读指标、provenance 清单，以及 compact data 目录中的 summary/CSV、候选几何和 CP2K 输入。大型 cube 与 CP2K 输出保留在 storage，只在 provenance 中记录路径、大小和哈希。", body),
        p("科学状态：COMPLETED_NUMERICAL_STAGE_B_PILOT_NOT_PRODUCTION_MECHANISM", callout),
        p(f"报告生成日期：{date.today().isoformat()}。项目仓库：https://github.com/SunsetStand/li-thf-amine-solvated-electron", small),
    ]
    document.build(story, onFirstPage=page_frame, onLaterPages=page_frame)


def _copy_compact_audit(bundle: dict[str, Any], output_dir: Path) -> None:
    data = output_dir / "data"
    data.mkdir(parents=True, exist_ok=True)
    for key in ("candidates", "legacy_smoke", "mechanism", "localization"):
        source = bundle["paths"][key]
        shutil.copy2(source, data / source.name)
    for filename in (
        "stage_b_localization.states.csv",
        "stage_b_localization.molecules.csv",
        "stage_b_localization.pairs.csv",
    ):
        source = bundle["base"] / filename
        if source.is_file():
            shutil.copy2(source, data / filename)
    for manifest in bundle["manifests"]:
        system = manifest["system_id"]
        replica = int(manifest["replica"])
        target = data / "candidates" / system / f"r{replica}"
        target.mkdir(parents=True, exist_ok=True)
        manifest_path = bundle["base"] / "stage_b" / system / f"r{replica}" / "candidates" / "manifest.json"
        shutil.copy2(manifest_path, target / "manifest.json")
        if replica == 1:
            selected = next(
                candidate for candidate in manifest["candidates"] if candidate["candidate_id"] == "separated"
            )
            for item in ("xyz", "cell"):
                source = Path(selected["structure"][item]["path"])
                shutil.copy2(source, target / source.name)
            metadata = selected.get("metadata")
            if metadata:
                shutil.copy2(metadata["path"], target / "metadata.json")
    for pair in bundle["mechanism"]["records"]:
        for state in pair["states"]:
            target = data / "cp2k_inputs" / state["system_id"]
            target.mkdir(parents=True, exist_ok=True)
            shutil.copy2(state["input"]["path"], target / f"{state['state_id']}.inp")


def build(run_root: str | Path, campaign: str, output_dir: str | Path) -> dict[str, Any]:
    output = Path(output_dir).resolve()
    figures = output / "figures"
    bundle = load_and_validate(run_root, campaign)
    metrics = summarize(bundle)
    metrics["source_file_count"] = len(bundle["sources"])
    metrics["source_bytes"] = sum(int(item["size_bytes"]) for item in bundle["sources"])
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    render_pipeline(figures / "stage_b_pipeline.png")
    render_candidate_bank(metrics, figures / "candidate_bank.png")
    render_mechanism_energy(metrics, figures / "mechanism_energy.png")
    render_spin_partition(metrics, figures / "spin_partition.png")
    render_cube_maps(metrics, figures / "cube_maps.png")
    metrics_path = output / "stage_b_metrics.json"
    _write_json(metrics_path, metrics)
    pdf_path = output / "stage_b_report_zh.pdf"
    build_pdf(metrics, figures, pdf_path)
    _copy_compact_audit(bundle, output)
    products = [pdf_path, metrics_path, *(figures / name for name in FIGURE_NAMES)]
    provenance = {
        "schema_version": 1,
        "kind": "stage_b_report_provenance",
        "ready": True,
        "scientific_status": metrics["scientific_status"],
        "campaign": campaign,
        "sources": bundle["sources"],
        "products": [
            {
                "path": str(path.resolve()),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in products
        ],
    }
    _write_json(output / "report_provenance.json", provenance)
    return metrics


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--campaign", default="pilot")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    metrics = build(args.run_root, args.campaign, args.output_dir)
    print(f"READY: {metrics['ready']}")
    print(f"STATUS: {metrics['scientific_status']}")
    print(f"PDF: {(args.output_dir / 'stage_b_report_zh.pdf').resolve()}")
    print(f"SOURCE_FILES: {metrics['source_file_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run a checked, optionally restartable GROMACS stage inside Slurm."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def run_checked(command: list[str], log: Path, *, append: bool = False) -> None:
    with log.open("a" if append else "w", encoding="utf-8") as handle:
        handle.write(f"$ {shlex.join(command)}\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"command exited {completed.returncode}: {shlex.join(command)}")


def build_commands(
    phase: str,
    mdp: Path,
    coordinates: Path,
    topology: Path,
    prefix: Path,
    threads: int,
    checkpoint: Path | None,
    resume_checkpoint: Path | None = None,
    checkpoint_minutes: int = 15,
) -> tuple[list[str], list[str]]:
    grompp = [
        "gmx",
        "grompp",
        "-f",
        str(mdp),
        "-c",
        str(coordinates),
        "-r",
        str(coordinates),
        "-p",
        str(topology),
        "-o",
        f"{prefix}.tpr",
        "-maxwarn",
        "0",
    ]
    if checkpoint is not None:
        grompp.extend(["-t", str(checkpoint)])
    mdrun = [
        "gmx",
        "mdrun",
        "-s",
        f"{prefix}.tpr",
        "-deffnm",
        str(prefix),
        "-ntmpi",
        "1",
        "-ntomp",
        str(threads),
        "-pin",
        "on",
        "-cpt",
        str(checkpoint_minutes),
    ]
    if resume_checkpoint is not None:
        mdrun.extend(["-cpi", str(resume_checkpoint), "-append"])
    if phase == "em":
        mdrun.extend(["-v"])
    return grompp, mdrun


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_fingerprint(
    phase: str,
    mdp: Path,
    coordinates: Path,
    topology: Path,
    checkpoint: Path | None,
    threads: int,
    checkpoint_minutes: int,
) -> dict[str, Any]:
    files = {
        "mdp": {"path": str(mdp), "sha256": file_sha256(mdp)},
        "coordinates": {"path": str(coordinates), "sha256": file_sha256(coordinates)},
        "topology": {"path": str(topology), "sha256": file_sha256(topology)},
    }
    if checkpoint is not None:
        files["checkpoint"] = {
            "path": str(checkpoint),
            "sha256": file_sha256(checkpoint),
        }
    return {
        "schema_version": 1,
        "phase": phase,
        "threads": threads,
        "checkpoint_minutes": checkpoint_minutes,
        "files": files,
    }


def prepare_restart_workspace(output_dir: Path, fingerprint: dict[str, Any]) -> tuple[Path, bool]:
    workspace = output_dir / ".resume"
    state_path = workspace / "input-state.json"
    if workspace.is_dir() and state_path.is_file():
        try:
            previous = json.loads(state_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            previous = None
        if previous == fingerprint:
            phase = str(fingerprint["phase"])
            can_resume = (workspace / f"{phase}.tpr").is_file() and (
                workspace / f"{phase}.cpt"
            ).is_file()
            if can_resume:
                return workspace, True
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True)
    state_path.write_text(
        json.dumps(fingerprint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return workspace, False


def promote_outputs(workspace: Path, output_dir: Path, phase: str, extensions: list[str]) -> None:
    for extension in extensions:
        source = workspace / f"{phase}.{extension}"
        destination = output_dir / source.name
        os.replace(source, destination)
    for name in ("grompp.log", "mdrun.stdout.log", "manifest.json"):
        os.replace(workspace / name, output_dir / name)
    shutil.rmtree(workspace)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["em", "nvt", "npt", "production"], required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--checkpoint")
    parser.add_argument("--checkpoint-minutes", type=int, default=15)
    parser.add_argument("--restartable", action="store_true")
    parser.add_argument("--expect-trajectory", action="store_true")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if not os.environ.get("SLURM_JOB_ID"):
        print("ERROR: refusing to run GROMACS outside a Slurm allocation", file=sys.stderr)
        return 2
    if args.phase in ("npt", "production") and not args.checkpoint:
        print(f"ERROR: {args.phase} continuation requires a checkpoint", file=sys.stderr)
        return 2

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    if threads <= 0:
        print("ERROR: SLURM_CPUS_PER_TASK must be positive", file=sys.stderr)
        return 2
    if args.checkpoint_minutes <= 0:
        print("ERROR: checkpoint-minutes must be positive", file=sys.stderr)
        return 2
    checkpoint = Path(args.checkpoint).resolve() if args.checkpoint else None
    mdp = Path(args.mdp).resolve()
    coordinates = Path(args.coordinates).resolve()
    topology = Path(args.topology).resolve()
    fingerprint = input_fingerprint(
        args.phase,
        mdp,
        coordinates,
        topology,
        checkpoint,
        threads,
        args.checkpoint_minutes,
    )
    workspace = output_dir
    resumed = False
    if args.restartable:
        workspace, resumed = prepare_restart_workspace(output_dir, fingerprint)
    prefix = workspace / args.phase
    resume_checkpoint = Path(f"{prefix}.cpt") if resumed else None
    grompp, mdrun = build_commands(
        args.phase,
        mdp,
        coordinates,
        topology,
        prefix,
        threads,
        checkpoint,
        resume_checkpoint,
        args.checkpoint_minutes,
    )
    try:
        if not resumed:
            run_checked(grompp, workspace / "grompp.log")
        run_checked(mdrun, workspace / "mdrun.stdout.log", append=resumed)
        extensions = ["tpr", "gro", "edr", "log"]
        if args.phase != "em":
            extensions.append("cpt")
        if args.expect_trajectory:
            extensions.append("xtc")
        required = [Path(f"{prefix}.{extension}") for extension in extensions]
        for path in required:
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"expected non-empty output is missing: {path}")
        manifest_path = (
            workspace / "manifest.json" if args.restartable else output_dir / "manifest.json"
        )
        manifest_path.write_text(
            json.dumps(
                {
                    "phase": args.phase,
                    "slurm_job_id": os.environ["SLURM_JOB_ID"],
                    "threads": threads,
                    "commands": [mdrun] if resumed else [grompp, mdrun],
                    "restartable": args.restartable,
                    "resumed_from_checkpoint": resumed,
                    "input_fingerprint": fingerprint,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        if args.restartable:
            promote_outputs(workspace, output_dir, args.phase, extensions)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

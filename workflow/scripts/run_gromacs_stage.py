#!/usr/bin/env python3
"""Run one checked GROMACS smoke stage inside an existing Slurm allocation."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


def run_checked(command: list[str], log: Path) -> None:
    with log.open("w", encoding="utf-8") as handle:
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
        "15",
    ]
    if phase == "em":
        mdrun.extend(["-v"])
    return grompp, mdrun


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["em", "nvt", "npt"], required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--checkpoint")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if not os.environ.get("SLURM_JOB_ID"):
        print("ERROR: refusing to run GROMACS outside a Slurm allocation", file=sys.stderr)
        return 2
    if args.phase == "npt" and not args.checkpoint:
        print("ERROR: NPT continuation requires an NVT checkpoint", file=sys.stderr)
        return 2

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / args.phase
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    if threads <= 0:
        print("ERROR: SLURM_CPUS_PER_TASK must be positive", file=sys.stderr)
        return 2
    checkpoint = Path(args.checkpoint).resolve() if args.checkpoint else None
    grompp, mdrun = build_commands(
        args.phase,
        Path(args.mdp).resolve(),
        Path(args.coordinates).resolve(),
        Path(args.topology).resolve(),
        prefix,
        threads,
        checkpoint,
    )
    try:
        run_checked(grompp, output_dir / "grompp.log")
        run_checked(mdrun, output_dir / "mdrun.stdout.log")
        required = [
            Path(f"{prefix}.tpr"),
            Path(f"{prefix}.gro"),
            Path(f"{prefix}.edr"),
            Path(f"{prefix}.log"),
        ]
        if args.phase != "em":
            required.append(Path(f"{prefix}.cpt"))
        for path in required:
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"expected non-empty output is missing: {path}")
        (output_dir / "manifest.json").write_text(
            json.dumps(
                {
                    "phase": args.phase,
                    "slurm_job_id": os.environ["SLURM_JOB_ID"],
                    "threads": threads,
                    "commands": [grompp, mdrun],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate one GAFF2/AM1-BCC molecule template with checked external tools."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path


def mol2_summary(path: Path) -> tuple[int, float]:
    atom_count = 0
    charge_sum = 0.0
    in_atoms = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if not in_atoms or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 9:
            raise ValueError(f"malformed MOL2 atom line: {line}")
        atom_count += 1
        charge_sum += float(fields[8])
    if atom_count == 0:
        raise ValueError(f"no atoms found in {path}")
    return atom_count, charge_sum


def run_checked(command: list[str], log: Path) -> None:
    with log.open("a", encoding="utf-8") as handle:
        handle.write(f"$ {shlex.join(command)}\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"command exited {completed.returncode}: {shlex.join(command)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--residue", required=True)
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.name
    smiles = output_dir / f"{stem}.smi"
    initial = output_dir / f"{stem}.sdf"
    mol2 = output_dir / f"{stem}.mol2"
    frcmod = output_dir / f"{stem}.frcmod"
    pdb = output_dir / f"{stem}.pdb"
    log = output_dir / "parameterize.log"
    manifest = output_dir / "manifest.json"

    smiles.write_text(f"{args.smiles} {stem}\n", encoding="utf-8")
    log.write_text("", encoding="utf-8")
    commands = [
        ["obabel", str(smiles), "-O", str(initial), "--gen3d", "-h"],
        [
            "antechamber",
            "-i",
            str(initial),
            "-fi",
            "sdf",
            "-o",
            str(mol2),
            "-fo",
            "mol2",
            "-at",
            "gaff2",
            "-c",
            "bcc",
            "-nc",
            str(args.charge),
            "-rn",
            args.residue,
            "-pf",
            "y",
        ],
        ["parmchk2", "-i", str(mol2), "-f", "mol2", "-o", str(frcmod), "-s", "gaff2"],
        [
            "antechamber",
            "-i",
            str(mol2),
            "-fi",
            "mol2",
            "-o",
            str(pdb),
            "-fo",
            "pdb",
            "-rn",
            args.residue,
            "-pf",
            "y",
        ],
    ]
    try:
        for command in commands:
            run_checked(command, log)
        for required in (initial, mol2, frcmod, pdb):
            if not required.is_file() or required.stat().st_size == 0:
                raise RuntimeError(f"expected non-empty output is missing: {required}")
        atom_count, charge_sum = mol2_summary(mol2)
        if abs(charge_sum - args.charge) > 0.02:
            raise RuntimeError(
                f"MOL2 charge sum {charge_sum:.6f} differs from requested charge {args.charge}"
            )
        manifest.write_text(
            json.dumps(
                {
                    "name": args.name,
                    "smiles": args.smiles,
                    "residue": args.residue,
                    "formal_charge": args.charge,
                    "mol2_charge_sum": charge_sum,
                    "atom_count": atom_count,
                    "commands": commands,
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

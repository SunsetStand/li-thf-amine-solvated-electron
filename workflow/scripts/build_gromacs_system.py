#!/usr/bin/env python3
"""Build an Amber/GAFF2 mixture with TLeap and convert it to GROMACS."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path


def render_leap_input(
    packed_pdb: Path,
    box_angstrom: float,
    molecules: list[tuple[str, Path, Path, int]],
    prmtop: Path,
    inpcrd: Path,
) -> str:
    lines = ["source leaprc.gaff2"]
    for index, (residue, mol2, frcmod, _count) in enumerate(molecules):
        lines.append(f'loadamberparams "{frcmod}"')
        lines.append(f'MOL{index} = loadmol2 "{mol2}"')
        lines.append(f"check MOL{index}")
        lines.append(f"# residue {residue}")
    lines.extend(
        [
            f'system = loadpdb "{packed_pdb}"',
            f"set system box {{ {box_angstrom:.8f} {box_angstrom:.8f} {box_angstrom:.8f} }}",
            "check system",
            f'saveamberparm system "{prmtop}" "{inpcrd}"',
            "quit",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packed-pdb", required=True)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--molecule",
        nargs=4,
        action="append",
        metavar=("RESIDUE", "MOL2", "FRCMOD", "COUNT"),
        required=True,
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    packed_pdb = Path(args.packed_pdb).resolve()
    specification = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    box_angstrom = float(specification["initial_box_angstrom"])
    if box_angstrom <= 0:
        print("ERROR: initial_box_angstrom must be positive", file=sys.stderr)
        return 2
    molecules = [
        (residue, Path(mol2).resolve(), Path(frcmod).resolve(), int(count))
        for residue, mol2, frcmod, count in args.molecule
    ]
    leap_input = output_dir / "tleap.in"
    leap_log = output_dir / "tleap.log"
    prmtop = output_dir / "system.prmtop"
    inpcrd = output_dir / "system.inpcrd"
    topology = output_dir / "topol.top"
    coordinates = output_dir / "conf.gro"
    manifest = output_dir / "manifest.json"
    leap_input.write_text(
        render_leap_input(packed_pdb, box_angstrom, molecules, prmtop, inpcrd),
        encoding="utf-8",
    )

    command = ["tleap", "-f", str(leap_input)]
    try:
        with leap_log.open("w", encoding="utf-8") as handle:
            handle.write(f"$ {shlex.join(command)}\n")
            handle.flush()
            completed = subprocess.run(
                command, cwd=output_dir, stdout=handle, stderr=subprocess.STDOUT, check=False
            )
        leap_text = leap_log.read_text(encoding="utf-8", errors="replace")
        if completed.returncode != 0 or "Errors = 0" not in leap_text or "FATAL" in leap_text:
            raise RuntimeError(f"TLeap validation failed; inspect {leap_log}")
        for required in (prmtop, inpcrd):
            if not required.is_file() or required.stat().st_size == 0:
                raise RuntimeError(f"expected non-empty output is missing: {required}")

        import parmed as pmd

        structure = pmd.load_file(str(prmtop), xyz=str(inpcrd))
        expected_residues = sum(count for _residue, _mol2, _frcmod, count in molecules)
        if len(structure.residues) != expected_residues:
            raise RuntimeError(
                f"expected {expected_residues} residues, found {len(structure.residues)}"
            )
        structure.save(str(topology), format="gromacs", overwrite=True)
        structure.save(str(coordinates), format="gro", overwrite=True)
        manifest.write_text(
            json.dumps(
                {
                    "packed_pdb": str(packed_pdb),
                    "box_angstrom": box_angstrom,
                    "atom_count": len(structure.atoms),
                    "residue_count": len(structure.residues),
                    "molecules": [
                        {
                            "residue": residue,
                            "mol2": str(mol2),
                            "frcmod": str(frcmod),
                            "count": count,
                        }
                        for residue, mol2, frcmod, count in molecules
                    ],
                    "commands": [command],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    except (ImportError, OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

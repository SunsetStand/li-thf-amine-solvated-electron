#!/usr/bin/env python3
"""Build an Amber/GAFF2 mixture with TLeap and convert it to GROMACS."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

TLEAP_IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
CHARGE_QUANTUM = Decimal("0.000001")
MAX_ROUNDING_RESIDUAL = Decimal("0.002")


def balanced_mol2(source: Path, destination: Path) -> tuple[Path, float, float]:
    """Balance only atom-charge rounding against the nearest integer charge.

    Keep the original AM1-BCC template untouched. TLeap receives a local copy
    only when its molecular charge has a small rounding residual.
    """
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    atoms: list[tuple[int, re.Match[str], Decimal]] = []
    in_atoms = False
    for index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if not in_atoms or not line.strip():
            continue
        fields = list(re.finditer(r"\S+", line))
        if len(fields) < 9:
            raise ValueError(f"malformed MOL2 atom line: {line.rstrip()}")
        try:
            charge = Decimal(fields[8].group())
        except InvalidOperation as exc:
            raise ValueError(f"invalid MOL2 atom charge: {line.rstrip()}") from exc
        atoms.append((index, fields[8], charge))
    if not atoms:
        raise ValueError(f"no atoms found in {source}")

    original = sum((charge for _index, _field, charge in atoms), Decimal(0))
    formal_charge = int(original.to_integral_value())
    residual = Decimal(formal_charge) - original
    if abs(residual) > MAX_ROUNDING_RESIDUAL:
        raise ValueError(
            f"MOL2 charge sum {original} in {source} differs from the nearest "
            "integer charge by more than rounding; inspect parameterization"
        )
    if residual == 0:
        return source, float(original), 0.0

    units = [int((charge / CHARGE_QUANTUM).to_integral_value()) for _, _, charge in atoms]
    remaining = formal_charge * 1_000_000 - sum(units)
    share, remainder = divmod(abs(remaining), len(atoms))
    direction = 1 if remaining >= 0 else -1
    corrected = [
        charge_units + direction * (share + (position < remainder))
        for position, charge_units in enumerate(units)
    ]
    maximum_change = max(
        abs(Decimal(new_units) * CHARGE_QUANTUM - old_charge)
        for new_units, (_index, _field, old_charge) in zip(corrected, atoms, strict=True)
    )
    for (index, field, _charge), new_units in zip(atoms, corrected, strict=True):
        line = lines[index]
        replacement = f"{Decimal(new_units) * CHARGE_QUANTUM:.6f}"
        lines[index] = line[: field.start()] + replacement + line[field.end() :]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines), encoding="utf-8")
    return destination, float(original), float(maximum_change)


def render_leap_input(
    packed_pdb: Path,
    box_angstrom: float,
    molecules: list[tuple[str, Path, Path, int]],
    prmtop: Path,
    inpcrd: Path,
) -> str:
    lines = ["source leaprc.gaff2"]
    registered_residues: set[str] = set()
    for residue, mol2, frcmod, _count in molecules:
        if not TLEAP_IDENTIFIER.fullmatch(residue):
            raise ValueError(f"residue is not a valid TLeap identifier: {residue!r}")
        if residue in registered_residues:
            raise ValueError(f"duplicate TLeap residue template: {residue}")
        registered_residues.add(residue)
        lines.append(f'loadamberparams "{frcmod}"')
        # loadPdb resolves a PDB residue by looking for a same-named LEaP
        # variable, so the template variable must be the configured residue.
        lines.append(f'{residue} = loadmol2 "{mol2}"')
        lines.append(f"check {residue}")
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
    command = ["tleap", "-f", str(leap_input)]
    try:
        leap_molecules = []
        charge_adjustments = []
        for residue, mol2, frcmod, count in molecules:
            used_mol2, raw_charge, maximum_change = balanced_mol2(
                mol2, output_dir / "charge_balanced" / f"{residue}.mol2"
            )
            leap_molecules.append((residue, used_mol2, frcmod, count))
            charge_adjustments.append(
                {
                    "residue": residue,
                    "source_mol2": str(mol2),
                    "tleap_mol2": str(used_mol2),
                    "raw_molecular_charge": raw_charge,
                    "maximum_atom_charge_change": maximum_change,
                }
            )
        leap_input.write_text(
            render_leap_input(packed_pdb, box_angstrom, leap_molecules, prmtop, inpcrd),
            encoding="utf-8",
        )
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
        for residue in structure.residues:
            residue_charge = sum(atom.charge for atom in residue.atoms)
            if abs(residue_charge - round(residue_charge)) > 1e-5:
                raise RuntimeError(
                    f"noninteger charge {residue_charge:.6f} remains in "
                    f"{residue.name} after TLeap; inspect {leap_log}"
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
                    "charge_adjustments": charge_adjustments,
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

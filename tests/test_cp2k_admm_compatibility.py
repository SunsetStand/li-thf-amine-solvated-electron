"""Validate generated ADMM blocks against the CP2K 2023.2 input reference."""

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from solvelec.config import load_config
from solvelec.li_relaxation import render_input, settings_from
from solvelec.rendering import render_cp2k, render_stage_b2c_method_benchmark_cp2k

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "workflow/templates/cp2k"
METHODS = load_config(ROOT / "configs/methods.yaml")

# CP2K 2023.2 official reference (independent of the input generators):
# https://manual.cp2k.org/cp2k-2023_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/
# AUXILIARY_DENSITY_MATRIX_METHOD.html
KEYWORDS_2023_2 = {
    "ADMM_PURIFICATION_METHOD": {
        "CAUCHY",
        "CAUCHY_SUBSPACE",
        "MCWEENY",
        "MO_DIAG",
        "MO_NO_DIAG",
        "NONE",
        "NONE_DM",
    },
    "BLOCK_LIST": None,
    "EPS_FILTER": None,
    "EXCH_CORRECTION_FUNC": {
        "BECKE88X",
        "BECKE88X_LIBXC",
        "DEFAULT",
        "DEFAULT_LIBXC",
        "LDA_X_LIBXC",
        "NONE",
        "OPTX",
        "OPTX_LIBXC",
        "PBEX",
        "PBEX_LIBXC",
    },
    "EXCH_SCALING_MODEL": {"MERLOT", "NONE"},
    "METHOD": {
        "BASIS_PROJECTION",
        "BLOCKED_PROJECTION",
        "BLOCKED_PROJECTION_PURIFY_FULL",
        "CHARGE_CONSTRAINED_PROJECTION",
    },
    "OPTX_A1": None,
    "OPTX_A2": None,
    "OPTX_GAMMA": None,
}


def check_admm(text: str) -> int:
    blocks = re.findall(
        r"&AUXILIARY_DENSITY_MATRIX_METHOD\s*\n(.*?)&END AUXILIARY_DENSITY_MATRIX_METHOD",
        text.upper(),
        re.S,
    )
    for block in blocks:
        for line in block.splitlines():
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            keyword = fields[0]
            if keyword not in KEYWORDS_2023_2:
                raise ValueError(f"CP2K 2023.2 unknown ADMM keyword: {keyword}")
            values = KEYWORDS_2023_2[keyword]
            if values is not None and (len(fields) != 2 or fields[1] not in values):
                raise ValueError(f"CP2K 2023.2 invalid ADMM value: {' '.join(fields)}")
    return len(blocks)


class CP2KADMMCompatibilityTests(unittest.TestCase):
    def test_rejects_observed_unknown_keyword_and_invalid_enum(self):
        for setting in ("ADMM_TYPE ADMMS", "METHOD ADMMS"):
            with self.subTest(setting=setting), self.assertRaises(ValueError):
                check_admm(
                    "&AUXILIARY_DENSITY_MATRIX_METHOD\n"
                    + setting
                    + "\n&END AUXILIARY_DENSITY_MATRIX_METHOD\n"
                )

    def test_li_states_optimizations_and_single_points(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coordinates = root / "coordinates.xyz"
            cell = root / "cell.inc"
            coordinates.write_text("4\nwater and Li\nO 1 1 1\nH 2 1 1\nH 1 2 1\nLi 5 5 5\n")
            cell.write_text("&CELL\nA 15 0 0\nB 0 15 0\nC 0 0 15\n&END CELL\n")
            settings = settings_from(METHODS)
            for variant in settings["variants"]:
                for state in ("a", "b"):
                    for optimize in (False, True):
                        with self.subTest(variant=variant, state=state, optimize=optimize):
                            output = root / "cp2k.inp"
                            render_input(
                                TEMPLATES / "li_relaxation.inp.tpl",
                                output,
                                coordinates=coordinates,
                                cell=cell,
                                settings=settings,
                                state=state,
                                variant=variant,
                                optimize=optimize,
                            )
                            expected = int(settings["variants"][variant]["xc_family"] == "PBE0")
                            self.assertEqual(check_admm(output.read_text()), expected)

    def test_other_hybrid_input_generators(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "cp2k.inp"
            render_cp2k(
                TEMPLATES / "pbe0_cdft.inp.tpl",
                output,
                "solvated_electron",
                "test",
                "coordinates.inc",
                "cell.inc",
                METHODS["cp2k"],
                li_atom_index=1,
                constrained=True,
            )
            self.assertEqual(check_admm(output.read_text()), 1)
            coordinates = root / "coordinates.xyz"
            coordinates.write_text("2\ntest\nH 1 1 1\nGh 5 5 5\n")
            cell = root / "cell.inc"
            cell.write_text("&CELL\nA 15 0 0\nB 0 15 0\nC 0 0 15\n&END CELL\n")
            base = METHODS["stage_b2c_preferential_smoke"]
            for variant in METHODS["stage_b2c_method_benchmark"]["variants"]:
                with self.subTest(variant=variant["id"]):
                    render_stage_b2c_method_benchmark_cp2k(
                        TEMPLATES / "stage_b2c_method_benchmark.inp.tpl",
                        output,
                        project="test",
                        coordinates_path=coordinates,
                        cell_path=cell,
                        method={
                            **base,
                            **variant,
                            "scientific_status": "NUMERICAL_METHOD_SENSITIVITY_BENCHMARK_ONLY",
                        },
                        state=base["states"][0],
                    )
                    self.assertEqual(
                        check_admm(output.read_text()),
                        int(variant["xc_family"] == "PBE0"),
                    )

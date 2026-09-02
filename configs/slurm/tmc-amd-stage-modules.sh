#!/usr/bin/env bash
# Source this file inside a TMC-AMD Slurm allocation. Each engine family brings
# a different OpenMPI build, so callers must select at most one engine stage.

# shellcheck source=tmc-amd-storage.sh
source "${ROOT}/configs/slurm/tmc-amd-storage.sh"

case "${1:-}" in
  molecule_generation|conformer_search)
    solvelec_activate_support_tools
    ;;
  classical_md)
    solvelec_activate_support_tools
    module load gromacs/2023
    ;;
  trajectory_analysis)
    # The repository virtual environment contains MDAnalysis/NumPy/SciPy.
    # Do not load an engine MPI family for read-only trajectory analysis.
    ;;
  cdft)
    module load cp2k/2023.2
    # The site module exposes cp2k.psmp but does not publish its data directory.
    # CP2K otherwise looks for these relative names in the calculation directory.
    : "${CP2K_DATA_DIR:=/data/softwares/cp2k/2023.2/data}"
    export CP2K_DATA_DIR
    for cp2k_data_file in BASIS_MOLOPT GTH_POTENTIALS dftd3.dat; do
      if [[ ! -r "${CP2K_DATA_DIR}/${cp2k_data_file}" ]]; then
        printf 'ERROR: required CP2K data file is not readable: %s\n' \
          "${CP2K_DATA_DIR}/${cp2k_data_file}" >&2
        return 2
      fi
    done
    ;;
  embedded_vde)
    module load orca/6.1.1
    ;;
  plane_wave_gate)
    solvelec_activate_qe
    ;;
  "")
    ;;
  *)
    printf 'ERROR: no TMC-AMD module mapping for stage %s\n' "$1" >&2
    return 2
    ;;
esac

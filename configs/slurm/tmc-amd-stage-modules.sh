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
  cdft)
    module load cp2k/2023.2
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

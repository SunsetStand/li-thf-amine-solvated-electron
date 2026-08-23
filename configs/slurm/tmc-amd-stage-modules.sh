#!/usr/bin/env bash
# Source this file inside a TMC-AMD Slurm allocation. Each engine family brings
# a different OpenMPI build, so callers must select at most one engine stage.

case "${1:-}" in
  classical_md)
    module load gromacs/2023
    ;;
  cdft)
    module load cp2k/2023.2
    ;;
  embedded_vde)
    module load orca/6.1.1
    ;;
  "")
    ;;
  *)
    printf 'ERROR: no TMC-AMD module mapping for stage %s\n' "$1" >&2
    return 2
    ;;
esac

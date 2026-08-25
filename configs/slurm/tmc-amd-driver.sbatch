#!/usr/bin/env bash
#SBATCH --job-name=solvelec-driver
#SBATCH --partition=amd
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail

ROOT="${SOLVELEC_ROOT:-${SLURM_SUBMIT_DIR:-}}"
if [[ -z "${ROOT}" || ! -f "${ROOT}/run.sh" ]]; then
  printf 'ERROR: repository root is invalid: %s\n' "${ROOT:-unset}" >&2
  printf 'SOLVELEC_ROOT=%s  SLURM_SUBMIT_DIR=%s\n' \
    "${SOLVELEC_ROOT:-unset}" "${SLURM_SUBMIT_DIR:-unset}" >&2
  exit 2
fi
ROOT="$(cd "${ROOT}" && pwd -P)"
if [[ $# -lt 1 ]]; then
  printf 'ERROR: expected a run.sh command.\n' >&2
  exit 2
fi

export MODULEPATH="/data/modulefiles/softwares:/data/modulefiles/libraries${MODULEPATH:+:${MODULEPATH}}"
source /etc/profile.d/modules.sh
module purge
if [[ ! -x "${ROOT}/.venv/bin/python" ]]; then
  module load miniconda3
fi

# Licensed programs and site-version choices belong in this untracked file.
SITE_MODULES="${SOLVELEC_SITE_MODULES:-${ROOT}/configs/profiles/private-tmc-amd/modules.sh}"
if [[ -f "${SITE_MODULES}" ]]; then
  # shellcheck source=/dev/null
  source "${SITE_MODULES}"
fi
export PYTHONUNBUFFERED=1
export SOLVELEC_DEFAULT_PROFILE="${SOLVELEC_DEFAULT_PROFILE:-tmc-amd}"
export SOLVELEC_REQUIRE_SLURM=1

load_doctor_stage_modules() {
  local selected_stage=""
  local previous=""
  local requirement=""
  for argument in "$@"; do
    requirement=""
    if [[ "${previous}" == "--require" ]]; then
      requirement="${argument}"
    elif [[ "${argument}" == --require=* ]]; then
      requirement="${argument#--require=}"
    elif [[ "${argument}" == "--strict-engines" ]]; then
      requirement="production"
    fi

    case "${requirement}" in
      molecule_generation|conformer_search|classical_md|cdft|embedded_vde|plane_wave_gate)
        if [[ -n "${selected_stage}" && "${selected_stage}" != "${requirement}" ]]; then
          printf 'ERROR: select one TMC capability stage per doctor job: %s and %s.\n' \
            "${selected_stage}" "${requirement}" >&2
          printf 'Run one doctor command per stage.\n' >&2
          return 2
        fi
        selected_stage="${requirement}"
        ;;
      production)
        printf 'ERROR: TMC production-wide doctor would mix incompatible MPI modules.\n' >&2
        printf 'Run the documented stage-specific doctor commands instead.\n' >&2
        return 2
        ;;
    esac
    previous="${argument}"
  done

  if [[ -n "${selected_stage}" ]]; then
    # shellcheck source=tmc-amd-stage-modules.sh
    source "${ROOT}/configs/slurm/tmc-amd-stage-modules.sh" "${selected_stage}"
  fi
}

cd "${ROOT}"
if [[ "$1" == "doctor" ]]; then
  load_doctor_stage_modules "${@:2}"
fi
printf '=== solvelec Slurm driver ===\n'
printf 'Date: %s  Node: %s  Job: %s\n' "$(date)" "$(hostname)" "${SLURM_JOB_ID}"
printf 'Command:'
printf ' %q' "${ROOT}/run.sh" "$@"
printf '\n\n'
exec "${ROOT}/run.sh" "$@"

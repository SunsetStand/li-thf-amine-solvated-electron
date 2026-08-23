#!/usr/bin/env bash
#SBATCH --job-name=solvelec-driver
#SBATCH --partition=amd
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
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

cd "${ROOT}"
printf '=== solvelec Slurm driver ===\n'
printf 'Date: %s  Node: %s  Job: %s\n' "$(date)" "$(hostname)" "${SLURM_JOB_ID}"
printf 'Command:'
printf ' %q' "${ROOT}/run.sh" "$@"
printf '\n\n'
exec "${ROOT}/run.sh" "$@"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
STAGE="${1:-}"
if [[ -z "${STAGE}" || "${2:-}" != "--" || $# -lt 3 ]]; then
  printf 'ERROR: usage: run_tmc_stage.sh STAGE -- COMMAND [ARG ...]\n' >&2
  exit 2
fi
shift 2

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf 'ERROR: refusing to run stage %s outside a Slurm allocation.\n' "${STAGE}" >&2
  exit 2
fi

export MODULEPATH="/data/modulefiles/softwares:/data/modulefiles/libraries${MODULEPATH:+:${MODULEPATH}}"
source /etc/profile.d/modules.sh
module purge

SITE_MODULES="${SOLVELEC_SITE_MODULES:-${ROOT}/configs/profiles/private-tmc-amd/modules.sh}"
if [[ -f "${SITE_MODULES}" ]]; then
  # shellcheck source=/dev/null
  source "${SITE_MODULES}"
fi

# shellcheck source=configs/slurm/tmc-amd-stage-modules.sh
source "${ROOT}/configs/slurm/tmc-amd-stage-modules.sh" "${STAGE}"

# Engine commands can change to a storage work directory before launch. Keep
# repository-provided site wrappers discoverable independently of that cwd.
export PATH="${ROOT}/configs/slurm:${PATH}"

printf 'Stage %s:' "${STAGE}"
printf ' %q' "$@"
printf '\n'
exec "$@"

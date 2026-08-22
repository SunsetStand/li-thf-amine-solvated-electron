#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SOURCE="${ROOT}/src"
PIP_PROJECT="${ROOT}"
if command -v cygpath >/dev/null 2>&1; then
  PYTHON_SOURCE="$(cygpath -w "${PYTHON_SOURCE}")"
  PIP_PROJECT="$(cygpath -w "${PIP_PROJECT}")"
  export PYTHONPATH="${PYTHON_SOURCE}${PYTHONPATH:+;${PYTHONPATH}}"
else
  export PYTHONPATH="${PYTHON_SOURCE}${PYTHONPATH:+:${PYTHONPATH}}"
fi

find_python() {
  if [[ -x "${ROOT}/.venv/bin/python" ]]; then
    printf '%s\n' "${ROOT}/.venv/bin/python"
  elif [[ -x "${ROOT}/.venv/Scripts/python.exe" ]]; then
    printf '%s\n' "${ROOT}/.venv/Scripts/python.exe"
  elif command -v python3 >/dev/null 2>&1; then
    command -v python3
  elif command -v python >/dev/null 2>&1; then
    command -v python
  else
    printf 'ERROR: Python 3.10+ was not found. Load a Python module first.\n' >&2
    return 1
  fi
}

find_snakemake() {
  if [[ -x "${ROOT}/.venv/bin/snakemake" ]]; then
    printf '%s\n' "${ROOT}/.venv/bin/snakemake"
  elif [[ -x "${ROOT}/.venv/Scripts/snakemake.exe" ]]; then
    printf '%s\n' "${ROOT}/.venv/Scripts/snakemake.exe"
  elif command -v snakemake >/dev/null 2>&1; then
    command -v snakemake
  else
    return 1
  fi
}

PYTHON_BIN="$(find_python)"
COMMAND="${1:-help}"
if [[ $# -gt 0 ]]; then shift; fi

campaign_from_args() {
  local campaign="pilot"
  local previous=""
  for item in "$@"; do
    if [[ "${previous}" == "--campaign" ]]; then campaign="${item}"; fi
    previous="${item}"
  done
  printf '%s\n' "${campaign}"
}

case "${COMMAND}" in
  bootstrap)
    if [[ ! -x "${ROOT}/.venv/bin/python" && ! -x "${ROOT}/.venv/Scripts/python.exe" ]]; then
      "${PYTHON_BIN}" -m venv "${ROOT}/.venv"
    fi
    PYTHON_BIN="$(find_python)"
    "${PYTHON_BIN}" -m pip install --upgrade pip
    extras="workflow,analysis,dev"
    if command -v sbatch >/dev/null 2>&1; then
      extras="${extras},slurm"
    fi
    "${PYTHON_BIN}" -m pip install -e "${PIP_PROJECT}[${extras}]"
    "${PYTHON_BIN}" -m solvelec.cli validate
    "${PYTHON_BIN}" -m solvelec.cli test
    ;;
  doctor)
    "${PYTHON_BIN}" -m solvelec.cli doctor "$@"
    ;;
  test)
    "${PYTHON_BIN}" -m solvelec.cli validate
    "${PYTHON_BIN}" -m solvelec.cli test
    ;;
  dry-run)
    campaign="$(campaign_from_args "$@")"
    if SNAKEMAKE_BIN="$(find_snakemake)"; then
      "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
        --config campaign="${campaign}" --dry-run --printshellcmds
    else
      printf 'Snakemake not installed; showing the expanded dependency matrix instead.\n' >&2
      "${PYTHON_BIN}" -m solvelec.cli matrix --campaign "${campaign}"
    fi
    ;;
  submit|resume)
    campaign="$(campaign_from_args "$@")"
    profile="slurm"
    previous=""
    for item in "$@"; do
      if [[ "${previous}" == "--profile" ]]; then profile="${item}"; fi
      previous="${item}"
    done
    if ! SNAKEMAKE_BIN="$(find_snakemake)"; then
      printf 'ERROR: snakemake is required. Run ./run.sh bootstrap or load its module.\n' >&2
      exit 2
    fi
    "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
      --config campaign="${campaign}" --profile "${ROOT}/configs/profiles/${profile}"
    ;;
  status)
    campaign="$(campaign_from_args "$@")"
    if ! SNAKEMAKE_BIN="$(find_snakemake)"; then
      printf 'ERROR: snakemake is required. Run ./run.sh bootstrap or load its module.\n' >&2
      exit 2
    fi
    "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
      --config campaign="${campaign}" --summary
    ;;
  report)
    campaign="$(campaign_from_args "$@")"
    "${PYTHON_BIN}" -m solvelec.cli report --campaign "${campaign}" \
      --output "${ROOT}/runs/${campaign}/report/README.md"
    ;;
  update)
    if [[ -n "$(git -C "${ROOT}" status --porcelain)" ]]; then
      printf 'ERROR: refusing to update a dirty worktree. Commit or stash changes first.\n' >&2
      exit 2
    fi
    git -C "${ROOT}" pull --ff-only
    exec "${ROOT}/run.sh" bootstrap
    ;;
  matrix)
    "${PYTHON_BIN}" -m solvelec.cli matrix "$@"
    ;;
  help|-h|--help)
    cat <<'EOF'
Usage: ./run.sh COMMAND [options]

Commands:
  bootstrap   Create/update .venv, install dependencies, and run tests
  doctor      Report engines; add --require STAGE to enforce a capability gate
  test        Validate configs and run dependency-light regression tests
  dry-run     Show Snakemake DAG, or the campaign matrix before bootstrap
  submit      Generate/submit the selected campaign through a profile
  resume      Resume incomplete Snakemake jobs using the selected profile
  status      Show Snakemake output status
  report      Generate a readiness report (never fabricates scientific results)
  matrix      Print the expanded composition/replica matrix
  update      Fast-forward only, sync dependencies, then rerun tests
EOF
    ;;
  *)
    printf 'ERROR: unknown command %s\n' "${COMMAND}" >&2
    exit 2
    ;;
esac

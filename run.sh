#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SOURCE="${ROOT}/src"
PIP_PROJECT="${ROOT}"
SLURM_DRIVER="${SOLVELEC_SLURM_DRIVER:-${ROOT}/configs/slurm/tmc-amd-driver.sbatch}"
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

COMMAND="${1:-help}"
if [[ $# -gt 0 ]]; then shift; fi

submit_via_slurm() {
  local action="$1"
  shift
  local runtime="${SOLVELEC_QUICK_TIME:-00:30:00}"
  local memory="${SOLVELEC_QUICK_MEM:-4G}"
  case "${action}" in
    bootstrap)
      runtime="${SOLVELEC_BOOTSTRAP_TIME:-02:00:00}"
      memory="${SOLVELEC_BOOTSTRAP_MEM:-8G}"
      ;;
    submit|resume)
      runtime="${SOLVELEC_CONTROLLER_TIME:-12:00:00}"
      memory="${SOLVELEC_CONTROLLER_MEM:-8G}"
      ;;
  esac

  if [[ ! -f "${SLURM_DRIVER}" ]]; then
    printf 'ERROR: Slurm driver not found: %s\n' "${SLURM_DRIVER}" >&2
    return 2
  fi
  mkdir -p "${ROOT}/runs/slurm"
  local submission
  submission="$(
    cd "${ROOT}"
    sbatch --parsable \
      --partition="${SOLVELEC_SLURM_PARTITION:-amd}" \
      --job-name="solvelec-${action}" \
      --ntasks=1 \
      --cpus-per-task="${SOLVELEC_SLURM_CPUS:-4}" \
      --time="${runtime}" \
      --mem="${memory}" \
      --output="${ROOT}/runs/slurm/%x-%j.out" \
      --error="${ROOT}/runs/slurm/%x-%j.err" \
      --export="ALL,SOLVELEC_REQUIRE_SLURM=1,SOLVELEC_ROOT=${ROOT}" \
      "${SLURM_DRIVER}" "${action}" "$@"
  )"
  printf 'Submitted %s as Slurm job %s.\n' "${action}" "${submission}"
  printf 'Logs: %s/runs/slurm/solvelec-%s-%s.{out,err}\n' \
    "${ROOT}" "${action}" "${submission%%;*}"
}

is_login_safe_command() {
  case "$1" in
    help|-h|--help|queue|logs|update) return 0 ;;
    *) return 1 ;;
  esac
}

is_known_command() {
  case "$1" in
    bootstrap|doctor|probe|test|dry-run|submit|resume|status|report|update|matrix|queue|logs|help|-h|--help)
      return 0
      ;;
    *) return 1 ;;
  esac
}

if ! is_known_command "${COMMAND}"; then
  printf 'ERROR: unknown command %s\n' "${COMMAND}" >&2
  exit 2
fi

# On a Slurm host, task-like commands are submission requests.  The same entry
# point still runs directly on workstations and in CI, where sbatch is absent.
if [[ -z "${SLURM_JOB_ID:-}" ]] && command -v sbatch >/dev/null 2>&1 \
  && ! is_login_safe_command "${COMMAND}"; then
  submit_via_slurm "${COMMAND}" "$@"
  exit $?
fi

if [[ "${SOLVELEC_REQUIRE_SLURM:-0}" == "1" && -z "${SLURM_JOB_ID:-}" ]] \
  && ! is_login_safe_command "${COMMAND}"; then
  printf 'ERROR: refusing to run %s outside a Slurm allocation.\n' "${COMMAND}" >&2
  exit 2
fi

PYTHON_BIN=""
case "${COMMAND}" in
  help|-h|--help|queue|logs|update|probe) ;;
  *) PYTHON_BIN="$(find_python)" ;;
esac

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
  probe)
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      printf 'ERROR: probe must run inside a Slurm allocation.\n' >&2
      exit 2
    fi
    printf '=== solvelec compute-node probe ===\n'
    printf 'Date: %s\n' "$(date --iso-8601=seconds 2>/dev/null || date)"
    printf 'Host: %s\n' "$(hostname)"
    printf 'Working directory: %s\n' "${ROOT}"
    printf 'Slurm job: %s\n' "${SLURM_JOB_ID:-missing}"
    printf 'CPUs: %s\n' "${SLURM_CPUS_ON_NODE:-unknown}"
    printf '\n=== Slurm allocation ===\n'
    scontrol show job "${SLURM_JOB_ID}" || true
    printf '\n=== Loaded modules ===\n'
    module list 2>&1 || true
    printf '\n=== Relevant available modules ===\n'
    module -t avail 2>&1 \
      | grep -Ei 'python|miniconda|apptainer|packmol|openbabel|amber|gromacs|cp2k|orca|espresso|xtb|crest|mpi' \
      || true
    printf '\n=== Candidate module activation (isolated) ===\n'
    for candidate in cp2k/2023.1_openmpi cp2k/2023.2 gromacs/2023 orca/6.1.1; do
      (
        printf '\n--- %s ---\n' "${candidate}"
        module purge
        if module load "${candidate}"; then
          module list 2>&1 || true
          for executable in cp2k.psmp cp2k.popt cp2k gmx gmx_mpi orca mpirun; do
            if command -v "${executable}" >/dev/null 2>&1; then
              printf '%-16s %s\n' "${executable}" "$(command -v "${executable}")"
            fi
          done
        else
          printf 'MODULE_LOAD_FAILED\n'
        fi
      )
    done
    printf '\n=== Executables ===\n'
    for executable in git python3 python apptainer singularity packmol obabel antechamber \
      gmx gmx_mpi cp2k.psmp cp2k.popt cp2k orca pw.x xtb crest sbatch srun sacct \
      mpirun; do
      printf '%-16s %s\n' "${executable}" "$(command -v "${executable}" || printf 'NOT_FOUND')"
    done
    printf '\n=== Filesystems ===\n'
    df -h "${ROOT}" /data/home/storage 2>/dev/null || true
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
    profile="${SOLVELEC_DEFAULT_PROFILE:-slurm}"
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
  queue)
    squeue -u "${USER}"
    ;;
  logs)
    if [[ $# -gt 1 ]]; then
      printf 'ERROR: usage: ./run.sh logs [command]\n' >&2
      exit 2
    fi
    log_filter="${1:-*}"
    case "${log_filter}" in
      \*|bootstrap|doctor|probe|test|dry-run|submit|resume|status|report|matrix) ;;
      *)
        printf 'ERROR: unknown log command %s\n' "${log_filter}" >&2
        exit 2
        ;;
    esac
    shopt -s nullglob
    if [[ "${log_filter}" == "*" ]]; then
      log_files=("${ROOT}/runs/slurm"/solvelec-*.out)
      log_files+=("${ROOT}/runs/slurm"/solvelec-*.err)
    else
      log_files=("${ROOT}/runs/slurm"/solvelec-"${log_filter}"-*.out)
      log_files+=("${ROOT}/runs/slurm"/solvelec-"${log_filter}"-*.err)
    fi
    if [[ ${#log_files[@]} -eq 0 ]]; then
      printf 'No matching logs under %s/runs/slurm\n' "${ROOT}"
    else
      ls -1t -- "${log_files[@]}"
    fi
    ;;
  matrix)
    "${PYTHON_BIN}" -m solvelec.cli matrix "$@"
    ;;
  help|-h|--help)
    cat <<'EOF'
Usage: ./run.sh COMMAND [options]

Commands:
  bootstrap   Install dependencies and test (auto-submitted on a Slurm host)
  probe       Inspect modules, executables, allocation, and filesystems on a compute node
  doctor      Check engines; add --require STAGE to enforce a capability gate
  test        Validate configs and run dependency-light regression tests
  dry-run     Build the Snakemake DAG without executing workflow jobs
  submit      Run the Snakemake controller, which submits child Slurm jobs
  resume      Resume through the same Slurm-controller mechanism
  status      Generate Snakemake output status inside a short Slurm job
  report      Generate a readiness report (never fabricates scientific results)
  matrix      Print the expanded composition/replica matrix
  update      Fast-forward locally, then submit bootstrap through Slurm
  queue       Show the current user's Slurm queue (login-node-safe)
  logs        List Slurm logs; optionally filter by command, e.g. logs probe

On a host with sbatch, every task-like command above is automatically wrapped
in configs/slurm/tmc-amd-driver.sbatch. Only help, queue, update's Git operation,
logs, and the sbatch submission itself run on the login node.
EOF
    ;;
  *)
    printf 'ERROR: unknown command %s\n' "${COMMAND}" >&2
    exit 2
    ;;
esac

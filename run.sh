#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SOURCE="${ROOT}/src"
PIP_PROJECT="${ROOT}"
SLURM_DRIVER="${SOLVELEC_SLURM_DRIVER:-${ROOT}/configs/slurm/tmc-amd-driver.sbatch}"
TMC_STORAGE_HELPER="${ROOT}/configs/slurm/tmc-amd-storage.sh"
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

profile_uses_slurm_executor() {
  local profile_directory="$1"
  local profile_config="${profile_directory}/config.v9+.yaml"
  [[ -f "${profile_config}" ]] \
    && grep -Eq '^[[:space:]]*executor:[[:space:]]*slurm([[:space:]]|$)' \
      "${profile_config}"
}

run_snakemake_controller() {
  local profile_directory="$1"
  shift

  if [[ -n "${SLURM_JOB_ID:-}" ]] && profile_uses_slurm_executor "${profile_directory}"; then
    # Snakemake's Slurm executor also removes inherited SLURM_* variables, but
    # older or mismatched plugin code can do so after its status thread starts.
    # That race submits the first jobs and then waits forever.  Clear the parent
    # allocation only in this subshell, before Snakemake can create any thread.
    # SOLVELEC_REQUIRE_SLURM remains exported, and every workflow rule receives
    # a fresh SLURM_JOB_ID from its own child allocation.
    (
      local parent_job_id="${SLURM_JOB_ID}"
      local variable
      local -a slurm_variables=()
      mapfile -t slurm_variables < <(compgen -A variable SLURM_)
      printf 'Preparing nested Slurm controller from job %s; clearing inherited SLURM_* context before Snakemake starts.\n' \
        "${parent_job_id}"
      printf 'Workflow rules remain protected and will run as Slurm child jobs.\n'
      for variable in "${slurm_variables[@]}"; do
        unset "${variable}"
      done
      "$@"
    )
  else
    "$@"
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
    tools-install)
      runtime="${SOLVELEC_TOOLS_TIME:-12:00:00}"
      memory="${SOLVELEC_TOOLS_MEM:-16G}"
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
    bootstrap|storage-init|tools-install|doctor|probe|test|dry-run|submit|resume|unlock|status|report|update|matrix|queue|logs|help|-h|--help)
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
  help|-h|--help|queue|logs|update|probe|storage-init|tools-install) ;;
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

workflow_target_from_args() {
  local target="input_bundle"
  local previous=""
  for item in "$@"; do
    if [[ "${previous}" == "--target" ]]; then target="${item}"; fi
    if [[ "${item}" == --target=* ]]; then target="${item#--target=}"; fi
    previous="${item}"
  done
  case "${target}" in
    input_bundle|classical_smoke) ;;
    *)
      printf 'ERROR: unknown workflow target %s; choose input_bundle or classical_smoke.\n' \
        "${target}" >&2
      return 2
      ;;
  esac
  printf '%s\n' "${target}"
}

workflow_storage_root() {
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # shellcheck source=configs/slurm/tmc-amd-storage.sh
    source "${TMC_STORAGE_HELPER}"
    solvelec_validated_storage_root
  else
    printf '%s\n' "${ROOT}"
  fi
}

workflow_run_root() {
  if [[ -n "${SOLVELEC_RUN_ROOT:-}" ]]; then
    printf '%s\n' "${SOLVELEC_RUN_ROOT%/}"
  else
    printf '%s/runs\n' "$(workflow_storage_root)"
  fi
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
  storage-init)
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      printf 'ERROR: storage-init must run inside a Slurm allocation.\n' >&2
      exit 2
    fi
    # shellcheck source=configs/slurm/tmc-amd-storage.sh
    source "${TMC_STORAGE_HELPER}"
    solvelec_init_storage
    ;;
  tools-install)
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      printf 'ERROR: tools-install must run inside a Slurm allocation.\n' >&2
      exit 2
    fi
    if [[ $# -gt 1 ]]; then
      printf 'ERROR: usage: ./run.sh tools-install [chem|ambertools|qe|all]\n' >&2
      exit 2
    fi
    tool_group="${1:-all}"
    case "${tool_group}" in
      chem|ambertools|qe|all) ;;
      *)
        printf 'ERROR: unknown tool group %s; choose chem, ambertools, qe, or all.\n' \
          "${tool_group}" >&2
        exit 2
        ;;
    esac
    # Load Conda only after this Bash process has started. Loading it in the
    # outer Slurm driver changes LD_LIBRARY_PATH before run.sh's shebang is
    # evaluated and makes Bash pick up Miniconda's incompatible libtinfo.
    module load miniconda3
    if ! command -v conda >/dev/null 2>&1; then
      printf 'ERROR: conda was not found; the TMC driver must load miniconda3.\n' >&2
      exit 2
    fi
    # shellcheck source=configs/slurm/tmc-amd-storage.sh
    source "${TMC_STORAGE_HELPER}"
    solvelec_init_storage
    storage_root="$(solvelec_validated_storage_root)"
    export CONDA_PKGS_DIRS="${storage_root}/cache/conda-pkgs"
    export CONDA_CHANNEL_PRIORITY=strict
    # The site Miniconda can predate the env subcommands' --yes option.
    # CONDA_ALWAYS_YES is supported by both old and current Conda releases.
    export CONDA_ALWAYS_YES=true

    install_tool_environment() {
      local environment_name="$1"
      local definition="$2"
      local prefix
      local manifest
      prefix="$(solvelec_tool_prefix "${environment_name}")"
      manifest="${storage_root}/software/manifests/${environment_name}.explicit.txt"
      printf '\n=== Installing %s ===\n' "${environment_name}"
      printf 'Definition: %s\nPrefix: %s\n' "${definition}" "${prefix}"
      if [[ -d "${prefix}/conda-meta" ]]; then
        conda env update --prefix "${prefix}" --file "${definition}"
      else
        conda env create --prefix "${prefix}" --file "${definition}"
      fi
      conda list --prefix "${prefix}" --explicit > "${manifest}.tmp"
      mv -- "${manifest}.tmp" "${manifest}"
      printf 'Manifest: %s\n' "${manifest}"
    }

    if [[ "${tool_group}" == "chem" || "${tool_group}" == "all" ]]; then
      install_tool_environment chem-tools \
        "${ROOT}/configs/environments/tmc-chem-tools.yml"
    fi
    if [[ "${tool_group}" == "ambertools" || "${tool_group}" == "all" ]]; then
      install_tool_environment ambertools \
        "${ROOT}/configs/environments/tmc-ambertools.yml"
    fi
    if [[ "${tool_group}" == "qe" || "${tool_group}" == "all" ]]; then
      install_tool_environment qe "${ROOT}/configs/environments/tmc-qe.yml"
    fi
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
    target="$(workflow_target_from_args "$@")"
    storage_root="$(workflow_storage_root)"
    run_root="$(workflow_run_root)"
    if SNAKEMAKE_BIN="$(find_snakemake)"; then
      "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
        "${target}" --config campaign="${campaign}" run_root="${run_root}" \
        storage_root="${storage_root}" --dry-run --printshellcmds
    else
      printf 'Snakemake not installed; showing the expanded dependency matrix instead.\n' >&2
      "${PYTHON_BIN}" -m solvelec.cli matrix --campaign "${campaign}"
    fi
    ;;
  submit|resume)
    campaign="$(campaign_from_args "$@")"
    target="$(workflow_target_from_args "$@")"
    storage_root="$(workflow_storage_root)"
    run_root="$(workflow_run_root)"
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
    profile_directory="${ROOT}/configs/profiles/${profile}"
    run_snakemake_controller "${profile_directory}" \
      "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
      "${target}" --config campaign="${campaign}" run_root="${run_root}" \
      storage_root="${storage_root}" --profile "${profile_directory}"
    ;;
  unlock)
    campaign="$(campaign_from_args "$@")"
    target="$(workflow_target_from_args "$@")"
    storage_root="$(workflow_storage_root)"
    run_root="$(workflow_run_root)"
    if ! SNAKEMAKE_BIN="$(find_snakemake)"; then
      printf 'ERROR: snakemake is required. Run ./run.sh bootstrap or load its module.\n' >&2
      exit 2
    fi
    if ! active_controllers="$(
      squeue -h -u "${USER}" -n solvelec-submit,solvelec-resume -o '%A %j %T'
    )"; then
      printf 'ERROR: unable to verify the Slurm controller queue; refusing to unlock.\n' >&2
      exit 2
    fi
    if [[ -n "${active_controllers}" ]]; then
      printf 'ERROR: refusing to unlock while workflow controllers are active:\n' >&2
      printf '%s\n' "${active_controllers}" >&2
      printf 'Cancel or wait for those jobs, then submit unlock again.\n' >&2
      exit 2
    fi
    "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
      "${target}" --config campaign="${campaign}" run_root="${run_root}" \
      storage_root="${storage_root}" --unlock
    ;;
  status)
    campaign="$(campaign_from_args "$@")"
    target="$(workflow_target_from_args "$@")"
    storage_root="$(workflow_storage_root)"
    run_root="$(workflow_run_root)"
    if ! SNAKEMAKE_BIN="$(find_snakemake)"; then
      printf 'ERROR: snakemake is required. Run ./run.sh bootstrap or load its module.\n' >&2
      exit 2
    fi
    "${SNAKEMAKE_BIN}" --snakefile "${ROOT}/workflow/Snakefile" --directory "${ROOT}" \
      "${target}" --config campaign="${campaign}" run_root="${run_root}" \
      storage_root="${storage_root}" --summary
    ;;
  report)
    campaign="$(campaign_from_args "$@")"
    run_root="$(workflow_run_root)"
    "${PYTHON_BIN}" -m solvelec.cli report --campaign "${campaign}" \
      --output "${run_root}/${campaign}/report/README.md"
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
      \*|bootstrap|storage-init|tools-install|doctor|probe|test|dry-run|submit|resume|unlock|status|report|matrix) ;;
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
  storage-init  Create the validated TMC project tree on /data/home/storage
  tools-install Install staged open-source engines; optional: chem|ambertools|qe|all
  probe       Inspect modules, executables, allocation, and filesystems on a compute node
  doctor      Check engines; add --require STAGE to enforce a capability gate
  test        Validate configs and run dependency-light regression tests
  dry-run     Build the Snakemake DAG; add --target classical_smoke for real-chain preview
  submit      Run the controller; --target defaults to the input_bundle
  resume      Resume through the same Slurm-controller mechanism
  unlock      Safely clear a stale Snakemake lock after checking no controller is active
  status      Generate target status inside a short Slurm job
  report      Generate a readiness report (never fabricates scientific results)
  matrix      Print the expanded composition/replica matrix
  update      Fast-forward locally, then submit bootstrap through Slurm
  queue       Show the current user's Slurm queue (login-node-safe)
  logs        List Slurm logs; optionally filter by command, e.g. logs probe

On a host with sbatch, every task-like command above is automatically wrapped
in configs/slurm/tmc-amd-driver.sbatch. Only help, queue, update's Git operation,
logs, and the sbatch submission itself run on the login node.

Large workflow outputs on TMC are written below
/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron/runs.
The repository remains the command working directory. The first executable
chemistry target is: ./run.sh submit --campaign smoke --target classical_smoke
EOF
    ;;
  *)
    printf 'ERROR: unknown command %s\n' "${COMMAND}" >&2
    exit 2
    ;;
esac

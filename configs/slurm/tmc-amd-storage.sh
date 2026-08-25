#!/usr/bin/env bash
# TMC-AMD storage and environment helpers. This file is sourced only inside
# Slurm allocations by the repository driver.

solvelec_storage_root() {
  local user_name="${USER:-}"
  if [[ -n "${SOLVELEC_STORAGE_ROOT:-}" ]]; then
    printf '%s\n' "${SOLVELEC_STORAGE_ROOT}"
  elif [[ -n "${user_name}" ]]; then
    printf '/data/home/storage/Backup_Data/%s/li-thf-amine-solvated-electron\n' \
      "${user_name}"
  else
    printf 'ERROR: USER is unset; set SOLVELEC_STORAGE_ROOT explicitly.\n' >&2
    return 2
  fi
}

solvelec_validated_storage_root() {
  local root
  local user_name="${USER:-}"
  root="$(solvelec_storage_root)" || return $?
  if [[ -z "${user_name}" ]]; then
    printf 'ERROR: USER is unset; cannot validate storage ownership.\n' >&2
    return 2
  fi
  case "/${root}/" in
    *"/../"*|*"/./"*)
      printf 'ERROR: storage root must not contain dot path components: %s\n' "${root}" >&2
      return 2
      ;;
  esac
  case "${root}" in
    "/data/home/storage/Backup_Data/${user_name}/"*) ;;
    *)
      printf 'ERROR: storage root must be below the current user directory: %s\n' \
        "/data/home/storage/Backup_Data/${user_name}" >&2
      printf 'Rejected storage root: %s\n' "${root}" >&2
      return 2
      ;;
  esac
  printf '%s\n' "${root%/}"
}

solvelec_init_storage() {
  local root
  root="$(solvelec_validated_storage_root)" || return $?
  mkdir -p -- \
    "${root}/runs" \
    "${root}/software/conda/envs" \
    "${root}/software/manifests" \
    "${root}/cache/conda-pkgs"
  printf 'Storage root: %s\n' "${root}"
  printf 'Production data: %s\n' "${root}/runs"
  printf 'Tool environments: %s\n' "${root}/software/conda/envs"
}

solvelec_tool_prefix() {
  local environment_name="${1:?environment name is required}"
  local root
  root="$(solvelec_validated_storage_root)" || return $?
  printf '%s/software/conda/envs/%s\n' "${root}" "${environment_name}"
}

solvelec_prepend_existing_bin() {
  local prefix="${1:?environment prefix is required}"
  if [[ -d "${prefix}/bin" ]]; then
    export PATH="${prefix}/bin:${PATH}"
  fi
}

solvelec_activate_support_tools() {
  local ambertools_prefix
  local chem_tools_prefix
  ambertools_prefix="$(solvelec_tool_prefix ambertools)" || return $?
  chem_tools_prefix="$(solvelec_tool_prefix chem-tools)" || return $?
  solvelec_prepend_existing_bin "${ambertools_prefix}"
  solvelec_prepend_existing_bin "${chem_tools_prefix}"
  if [[ -d "${ambertools_prefix}" ]]; then
    export AMBERHOME="${ambertools_prefix}"
  fi
}

solvelec_activate_qe() {
  local qe_prefix
  qe_prefix="$(solvelec_tool_prefix qe)" || return $?
  solvelec_prepend_existing_bin "${qe_prefix}"
  if [[ -d "${qe_prefix}/lib" ]]; then
    export LD_LIBRARY_PATH="${qe_prefix}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    export OPAL_PREFIX="${qe_prefix}"
  fi
}

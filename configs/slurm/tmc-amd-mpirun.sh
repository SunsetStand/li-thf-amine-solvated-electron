#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 -n RANKS COMMAND [ARG ...]" >&2
}

die() {
    echo "ERROR: TMC MPI launcher: $*" >&2
    exit 2
}

[[ -n "${SLURM_JOB_ID:-}" ]] || die "refusing to run outside a Slurm allocation"
command -v scontrol >/dev/null 2>&1 || die "scontrol is unavailable"
command -v mpirun >/dev/null 2>&1 || die "mpirun is unavailable"

case "${1:-}" in
    -n|--np)
        [[ $# -ge 3 ]] || {
            usage
            exit 2
        }
        ranks="$2"
        shift 2
        ;;
    *)
        usage
        exit 2
        ;;
esac

[[ "$ranks" =~ ^[1-9][0-9]*$ ]] || die "rank count must be a positive integer"

# Snakemake's slurm-jobstep executor deliberately removes allocation-size
# variables before running MPI rules. Query Slurm itself instead of trusting
# inherited environment values, then enforce the allocation before launching.
job_record="$(scontrol show job --oneliner "$SLURM_JOB_ID")" \
    || die "cannot inspect Slurm job ${SLURM_JOB_ID}"

job_state=""
num_nodes=""
num_cpus=""
batch_host=""
node_list=""
for field in $job_record; do
    case "$field" in
        JobState=*) job_state="${field#*=}" ;;
        NumNodes=*) num_nodes="${field#*=}" ;;
        NumCPUs=*) num_cpus="${field#*=}" ;;
        BatchHost=*) batch_host="${field#*=}" ;;
        NodeList=*) node_list="${field#*=}" ;;
    esac
done

[[ "$job_state" == "RUNNING" ]] || die "Slurm job is not RUNNING (state=${job_state:-unknown})"
[[ "$num_nodes" == "1" ]] || die "only single-node allocations are supported (nodes=${num_nodes:-unknown})"
[[ "$num_cpus" =~ ^[1-9][0-9]*$ ]] || die "Slurm did not report a valid CPU allocation"
(( ranks <= num_cpus )) \
    || die "requested ${ranks} ranks but Slurm allocated only ${num_cpus} CPUs"

local_host="$(hostname -s)"
allocated_host="${batch_host:-$node_list}"
[[ -n "$allocated_host" && "$allocated_host" != "(null)" ]] \
    || die "Slurm did not report the allocated host"
[[ "$local_host" == "$allocated_host" ]] \
    || die "current host ${local_host} is not allocated host ${allocated_host}"

echo "TMC MPI launcher: job=${SLURM_JOB_ID} host=${local_host} ranks=${ranks}/${num_cpus}" >&2

# The site's OpenMPI 4.1.5 lacks the PMI interface required for direct `srun`
# launch.  The allocation is single-node, so let mpirun fork locally inside the
# existing batch cgroup.  Removing Slurm discovery variables prevents OpenMPI
# from mistaking Snakemake's one-process wrapper context for the full allocation.
# The explicit host slot count and --nooversubscribe retain the Slurm CPU cap.
while IFS='=' read -r variable _value; do
    if [[ "$variable" == SLURM_* ]]; then
        unset "$variable"
    fi
done < <(env)

exec mpirun \
    --host "${local_host}:${ranks}" \
    --map-by slot \
    --bind-to core \
    --nooversubscribe \
    -n "$ranks" \
    "$@"

#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --ntasks=1
#SBATCH --partition=xeon
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2%2
#SBATCH --mem=10G

set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------
chunks=(4 6)
studies=(7 7)

n_jobs=${#chunks[@]}

idx=$((SLURM_ARRAY_TASK_ID - 1))

if [ "$idx" -ge "$n_jobs" ]; then
  echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
  exit 1
fi

chunk=${chunks[$idx]}
study=${studies[$idx]}

# -----------------------------
# Log file (custom)
# -----------------------------
mkdir -p logs
logfile="logs/coloc_chunk${chunk}_study${study}.log"

# Reset log
: > "$logfile"

# -----------------------------
# Job metadata
# -----------------------------
echo "----------------------------------------" >> "$logfile"
echo "MR JOB STARTED" >> "$logfile"
echo "Date: $(date)" >> "$logfile"
echo "SLURM_JOB_ID: $SLURM_JOB_ID" >> "$logfile"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID" >> "$logfile"
echo "Running chunk index: $chunk" >> "$logfile"
echo "Running GWAS: $study" >> "$logfile"
echo "Node: $(hostname)" >> "$logfile"
echo "----------------------------------------" >> "$logfile"

# -----------------------------
# Bind mounts
# -----------------------------
export APPTAINER_BIND="/dev/shm:/dev/shm,/mnt/proteome:/mnt/proteome,/tmp:/tmp"

# -----------------------------
# Container path (USE ABSOLUTE PATH)
# -----------------------------
CONTAINER_PATH="/mnt/nfs_home/smartinez/images/mr_pipeline_250603.sif"

# -----------------------------
# Script path
# -----------------------------
SCRIPT_PATH="/mnt/nfs_home/smartinez/Bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_colocalization.R"

# -----------------------------
# Run
# -----------------------------
apptainer exec \
  --bind /mnt/proteome:/mnt/proteome \
  "$CONTAINER_PATH" \
  Rscript "$SCRIPT_PATH" "$chunk" "$study" \
  > "$logfile" 2>&1

echo "=== END $(date) ==="
#!/bin/bash
#SBATCH --job-name=pMR_est
#SBATCH --partition=xeon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6%6
#SBATCH --mem=10G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# -----------------------------
# Create log directory
# -----------------------------
mkdir -p logs

# -----------------------------
# Container
# -----------------------------
CONTAINER_PATH=/mnt/nfs_home/smartinez/images/mr_pipeline_250603.sif

# -----------------------------
# Bind mounts
# -----------------------------
export APPTAINER_BIND="/dev/shm:/dev/shm,/mnt/proteome:/mnt/proteome,/tmp:/tmp"

# -----------------------------
# Study × instrument set mapping
# -----------------------------
param1_array=(1 1 1 7 7 7)
param2_array=(1 2 3 1 2 3)

param1=${param1_array[$SLURM_ARRAY_TASK_ID-1]}
param2=${param2_array[$SLURM_ARRAY_TASK_ID-1]}

# -----------------------------
# Log file
# -----------------------------
logfile="logs/pMR_${param1}_${param2}.log"

# erase previous logfile
: > "$logfile"

# -----------------------------
# Job metadata
# -----------------------------
echo "----------------------------------------" >> "$logfile"
echo "MR JOB STARTED" >> "$logfile"
echo "Date: $(date)" >> "$logfile"
echo "SLURM_JOB_ID: $SLURM_JOB_ID" >> "$logfile"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID" >> "$logfile"
echo "Study index: $param1" >> "$logfile"
echo "Instrument set: $param2" >> "$logfile"
echo "Node: $(hostname)" >> "$logfile"
echo "----------------------------------------" >> "$logfile"

# -----------------------------
# Run MR pipeline
# -----------------------------
apptainer exec \
  --bind /mnt/proteome:/mnt/proteome \
  $CONTAINER_PATH \
  Rscript \
  ~/Bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/05_estimation.R \
  $param1 $param2 \
  &>> "$logfile" 2>&1

# -----------------------------
# Finish timestamp
# -----------------------------
echo "----------------------------------------" >> "$logfile"
echo "MR JOB FINISHED" >> "$logfile"
echo "Date: $(date)" >> "$logfile"
echo "----------------------------------------" >> "$logfile"
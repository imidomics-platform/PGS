#!/bin/bash
#SBATCH --job-name=MR_est
#SBATCH --partition=xeon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Mode from launcher
MODE=${MODE:-pQTL}

# Paths
CONTAINER_PATH=/mnt/nfs_home/smartinez/images/mr_pipeline_250603.sif
SCRIPT_PATH=~/Bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/MR_pipeline/05_estimation.R

# Bind mounts
export APPTAINER_BIND="/dev/shm:/dev/shm,/mnt/proteome:/mnt/proteome,/tmp:/tmp"

# Logs
mkdir -p logs

# Number of IS
n_is=3

# Get number of studies again (safe inside worker)
n_studies=$(apptainer exec $CONTAINER_PATH Rscript -e "
source('~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R');
cfg <- load_mr_config('${MODE}');
cat(length(cfg\$gwas))
")

# Mapping
idx=$((SLURM_ARRAY_TASK_ID - 1))

study_index=$(( idx / n_is + 1 ))
is_index=$(( idx % n_is + 1 ))

# Log file
logfile="logs/${MODE}_MR_s${study_index}_IS${is_index}.log"
: > "$logfile"

# Metadata
echo "----------------------------------------" >> "$logfile"
echo "MR JOB STARTED" >> "$logfile"
echo "Mode: $MODE" >> "$logfile"
echo "Date: $(date)" >> "$logfile"
echo "SLURM_JOB_ID: $SLURM_JOB_ID" >> "$logfile"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID" >> "$logfile"
echo "Study index: $study_index" >> "$logfile"
echo "Instrument set: $is_index" >> "$logfile"
echo "Node: $(hostname)" >> "$logfile"
echo "----------------------------------------" >> "$logfile"

# Run pipeline
apptainer exec \
  --bind /mnt/proteome:/mnt/proteome \
  $CONTAINER_PATH \
  Rscript \
  "$SCRIPT_PATH" \
  "$study_index" "$is_index" "$MODE"\
  &>> "$logfile"

# Finish
echo "----------------------------------------" >> "$logfile"
echo "MR JOB FINISHED" >> "$logfile"
echo "Date: $(date)" >> "$logfile"
echo "----------------------------------------" >> "$logfile"
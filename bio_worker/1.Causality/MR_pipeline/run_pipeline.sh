#!/bin/bash

set -euo pipefail

# Arguments ----

MODE=${1:-pQTL}
PROFILE=${2:-slurm,apptainer}

# Paths ----

PIPELINE_DIR=~/Bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/MR_pipeline
MAIN_NF=${PIPELINE_DIR}/main.nf
CONFIG_FILE=${PIPELINE_DIR}/nextflow.config

# Validation ----

if [[ "$MODE" != "pQTL" && "$MODE" != "eQTL" ]]; then
  echo "Invalid mode: $MODE"
  echo "Allowed values: pQTL | eQTL"
  exit 1
fi

if [[ ! -f "$MAIN_NF" ]]; then
  echo "Cannot find main.nf at: $MAIN_NF"
  exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Cannot find nextflow.config at: $CONFIG_FILE"
  exit 1
fi

# Logging ----

mkdir -p "${PIPELINE_DIR}/logs"

RUNSTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="${PIPELINE_DIR}/logs/run_pipeline_${MODE}_${RUNSTAMP}.log"

: > "$LOGFILE"

echo "----------------------------------------" >> "$LOGFILE"
echo "NEXTFLOW PIPELINE STARTED" >> "$LOGFILE"
echo "Mode: $MODE" >> "$LOGFILE"
echo "Profile: $PROFILE" >> "$LOGFILE"
echo "Date: $(date)" >> "$LOGFILE"
echo "Host: $(hostname)" >> "$LOGFILE"
echo "Working dir: $(pwd)" >> "$LOGFILE"
echo "Pipeline dir: $PIPELINE_DIR" >> "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"

# Run pipeline ----

nextflow run "$MAIN_NF" \
  -config "$CONFIG_FILE" \
  -profile "$PROFILE" \
  --mode "$MODE" \
  &>> "$LOGFILE"

# Finish ----

echo "----------------------------------------" >> "$LOGFILE"
echo "NEXTFLOW PIPELINE FINISHED" >> "$LOGFILE"
echo "Date: $(date)" >> "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"

echo "Pipeline submitted/executed successfully."
echo "Log: $LOGFILE"
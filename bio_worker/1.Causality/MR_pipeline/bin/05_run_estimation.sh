#!/bin/bash

set -euo pipefail

MODE=${1:-pQTL}

CONTAINER_PATH=/mnt/nfs_home/smartinez/images/mr_pipeline_250603.sif

# Get number of studies from config
n_studies=$(apptainer exec $CONTAINER_PATH Rscript -e "
source('~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R');
cfg <- load_mr_config('${MODE}');
cat(length(cfg\$gwas))
")

n_is=3
total_jobs=$(( n_studies * n_is ))

echo "Mode: $MODE"
echo "Studies: $n_studies"
echo "Total jobs: $total_jobs"

# Submit SLURM job with dynamic array
sbatch \
  --array=1-${total_jobs}%6 \
  --export=MODE=${MODE} \
  05_run_estimation_worker.sh
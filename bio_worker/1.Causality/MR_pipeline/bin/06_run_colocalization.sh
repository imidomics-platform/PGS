#!/bin/bash

set -euo pipefail

MODE=${1:-pQTL}

CONTAINER_PATH=/mnt/nfs_home/smartinez/images/mr_pipeline_250603.sif

n_studies=$(apptainer exec "$CONTAINER_PATH" Rscript -e "
source('~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R');
cfg <- load_mr_config('${MODE}');
cat(length(cfg\$gwas))
")

n_chunks=$(apptainer exec "$CONTAINER_PATH" Rscript -e "
source('~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R');
cfg <- load_mr_config('${MODE}');
cat(cfg\$n_chunks)
")

total_jobs=$(( n_studies * n_chunks ))

echo "Mode: $MODE"
echo "Studies: $n_studies"
echo "Chunks: $n_chunks"
echo "Total jobs: $total_jobs"

sbatch \
--array=1-"${total_jobs}"%6 \
--export=MODE="${MODE}" \
06_run_colocalization_worker.sh
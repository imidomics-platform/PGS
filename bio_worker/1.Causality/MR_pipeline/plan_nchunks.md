# Plan: Add nchunks to estimate step

**TL;DR:** Adapt the ESTIMATE process to use nchunks like COLOC does—partitioning genes across parallel tasks while keeping n_is iterations. This creates studies × chunks × is_idx jobs, each processing a subset of genes. Mirror the chunking pattern from coloc: split gene list by `cut()`, select chunk-specific subset, and loop over it.

## Discovery Notes

**Current Structure:**
- COLOC process: Input `tuple val(chunk), val(study)`, passes chunk + n_chunks to 06_colocalization.R
- ESTIMATE process: Input `tuple val(study), val(is_idx)`, no chunking; loops `for (gene in genes) { ... }`
- Channels: `coloc_jobs_ch = studies × 1..nchunks`, `estimate_jobs_ch = studies × 1..n_is`

**User Intent:** Keep both chunking AND n_is iterations (studies × chunks × is_idx), with genes partitioned by chunk for parallelization.

## Implementation Plan

### Phase 1: Update main.nf workflow
1. Modify `estimate_jobs_ch` channel to include chunks: `studies × 1..nchunks × 1..n_is`
2. Adjust tuple mapping to `tuple(chunk, study, is_idx)`

### Phase 2: Update ESTIMATE process module
1. Change input tuple from `tuple val(study), val(is_idx)` → `tuple val(chunk), val(study), val(is_idx)`
2. Update tag to include chunk: `"${study}_chunk${chunk}_IS${is_idx}"`
3. Update log output tag similarly
4. Pass chunk and n_chunks to R script: add `${chunk}` and `${params.n_chunks}` arguments

### Phase 3: Adapt 05_estimation.R
1. Parse new command-line args: [chunk, study_index, is_index, mode, n_chunks]
2. Add arg validation for chunk (similar to coloc: lines 65-66)
3. After loading genes (line 89), add chunking logic:
   ```R
   chunks <- split(genes, cut(seq_along(genes), n_chunks, labels = FALSE))
   chosen_genes <- chunks[[chunk]]
   ```
4. Update loop: `for (gene in chosen_genes)` instead of `for (gene in genes)`
5. Update logging to show chunk info (like coloc line 77)

## Relevant Files
- `Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/MR_pipeline/main.nf` — update workflow and estimate_jobs_ch channel
- `Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/MR_pipeline/modules/local/estimate.nf` — update ESTIMATE process input/tag/script
- `Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/MR_pipeline/bin/05_estimation.R` — add chunk parsing, validation, and gene partitioning logic

## Verification
1. Test with mock data: ensure chunk partitioning produces same gene count distribution as coloc
2. Run pipeline with n_chunks=2 on single study → verify 2 chunks × 3 IS = 6 estimate jobs per study
3. Confirm all genes across chunks are processed exactly once per IS iteration
4. Check log files confirm each job shows correct chunk/n_chunks info
5. Ensure output file naming remains consistent with current convention (adjust if needed)

## Dependencies
- None — all changes are independent (can be done in parallel or sequentially)

## Status
- [ ] Phase 1 complete
- [ ] Phase 2 complete
- [ ] Phase 3 complete
- [ ] Verification complete

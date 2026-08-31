#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(stringr)
  library(jsonlite)
})

set.seed(1)

# Project settings ----
proj_dir <- "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_MR_pQTL/"

# Load shared configuration ----
source("../../../Shared/imidomics/R/readme_generator.R")
source("../../../Shared/imidomics/R/config.R")

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir  <- imidomics_config$raw_data_dir

# Define paths ----
meta_dir <- file.path(
  root_data_dir,
  "/Projects/Internal_Projects/2021_ExternalDatasets_Mining/_MultiOmics/01_crossTalk/pQTL/MetaAnalysis/Results/"
)

out_base <- file.path(root_data_dir, proj_dir, "01_pQTLs_Preprocessed")
dir_clean <- file.path(out_base, "01_MetaAnalyzed_pQTLs_Clean")
dir_formatted <- file.path(out_base, "02_Filtered_Formatted")

dir.create(dir_clean, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_formatted, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_base, "log.txt")

tic <- Sys.time()

# Load required objects ----
featureTable <- readRDS(file.path(
  root_data_dir,
  "/Data/Molecular/Proteomics/Plasma/2021_IMX_P4V4_SOMAscan41_IA/Metadata/featureTable.rds"
))

# Identify meta-analysis files ----
files <- list.files(meta_dir, pattern = "tbl$", full.names = FALSE)
processed_files <- character(0)

message("Found ", length(files), " meta-analysis files.")
 
# Main loop ----
for (file in files) {
  
  prot_Apt <- gsub("1.tbl", "", file)
  protApt  <- gsub("_", ".", paste0("seq_", prot_Apt))
  protSymb <- featureTable$EntrezGeneSymbol[featureTable$AptName == protApt]
  
  RES <- tryCatch(
    fread(file.path(meta_dir, paste0(prot_Apt, "1.tbl"))),
    error = function(e) {
      write.table(
        cbind(prot_Apt, " failed to read"),
        file = log_file, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE, append = TRUE
      )
      return(NULL)
    }
  )
  
  if (is.null(RES) || nrow(RES) == 0) {
    write.table(
      cbind(protSymb, " has 0 rows after METAL"),
      file = log_file, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE, append = TRUE
    )
    next()
  }
  
  # Direction consistency ----
  RES$NPos <- str_count(RES$Direction, "\\+")
  RES$NNeg <- str_count(RES$Direction, "-")
  RES$NConsistent <- RES$NPos + RES$NNeg
  
  RES <- RES[which(RES$NConsistent != 0),]
  
  # Save CLEAN version (reverse MR usage)
  fwrite(RES, file = file.path(dir_clean, paste0(prot_Apt, ".txt")))
  
  # Instrument candidate filtering ----
  idx <- which(RES$NPos > 0 & RES$NNeg > 0)
  if (length(idx) > 0) RES <- RES[-idx, ]
  
  if (max(RES$NConsistent) > 1) {
    RES <- RES[which(RES$HetPVal > 5E-8 | RES$NConsistent >= 3),]
  }
  
  # P-value handling ----
  pvals <- RES$`P-value`
  RES$`P-value` <- as.numeric(RES$`P-value`)
  RES$LP <- -log10(RES$`P-value`)
  RES$LP[RES$LP == Inf] <- as.numeric(gsub(".*e-", "", pvals[RES$LP == Inf]))
  
  RES <- RES[which(RES$`P-value` < 0.05),]
  
  if (nrow(RES) == 0) next()
  
  if (sum(RES$`P-value` < 10^-315) > 0) {
    idx <- which(RES$`P-value` < 10^-315)
    lpvalues <- RES$LP[idx]
    idx_ordered <- order(lpvalues, decreasing = FALSE)
    newvalues <- 1 / (1:length(lpvalues)) * 10^-315
    RES$`P-value`[idx[idx_ordered]] <- newvalues
  }
  
  # Convert to TwoSampleMR format ----
  RES$exposure <- prot_Apt
  
  exp.dat <- format_data(
    dat = RES,
    snp_col = "MarkerName",
    beta_col = "Effect",
    se_col = "StdErr",
    eaf_col = "Freq1",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    pval_col = "P-value",
    samplesize_col = "SampleSize",
    chr_col = "Chromosome",
    pos_col = "Position",
    log_pval = FALSE,
    min_pval = 1e-320,
    phenotype_col = "exposure"
  )
  
  exp.dat$NConsistent <- RES$NConsistent
  exp.dat$exposure.Direction <- RES$Direction
  exp.dat$LP <- RES$LP
  
  saveRDS(
    exp.dat,
    file = file.path(dir_formatted, paste0(prot_Apt, ".rds"))
  )
  
  processed_files <- c(processed_files, prot_Apt)
}

# Save file list for clumping ----
saveRDS(
  processed_files,
  file = file.path(out_base, "files_for_clumping.rds")
)

# README (human-readable provenance) ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 1,
  titol = "pQTL data preprocessing",
  autors = "Sergio",
  lloc = "01_pQTL_preprocessing.R",
  temps = temps,
  descripcio = "First we simply remove variants that have beta equal to 0 avoid numeric problems. This is kept in 01_pQTLs_Preprocessed/01_MetaAnalyzed_pQTLs_Clean. The second step is to format adequately the P-values to avoid 0s and preserve the significance ranking, keep the pQTLs that were direction-concordant in the meta-analysis, filter by heterogeneity p-value and keep only the nominally significant results to reduce the size of the files. The data is stored in TwoSampleMR format in the folder 01_pQTLs_Preprocessed/02_Filtered_Formatted.",
  altres = ""
)

# JSON metadata (machine-readable) ----
metadata <- list(
  script = "01_pQTL_preprocessing.R",
  timestamp = as.character(Sys.time()),
  n_files_input = length(files),
  n_files_processed = length(processed_files),
  R_version = R.version.string
)

write_json(
  metadata,
  path = file.path(out_base, "run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Step 01 completed successfully.")
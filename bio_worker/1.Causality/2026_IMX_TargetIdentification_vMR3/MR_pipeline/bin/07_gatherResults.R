#!/usr/bin/env Rscript

# Arguments ----

args <- commandArgs(trailingOnly = TRUE)

mode <- ifelse(length(args) >= 1, args[1], "pQTL")
study_filter <- ifelse(length(args) >= 2, args[2], NA_character_)

# Libraries ----

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Logging ----

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

# Config ----

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
source("~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R")

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir

cfg <- load_mr_config(mode)

proj_dir <- cfg$proj_dir
studies <- names(cfg$gwas)

if (!is.na(study_filter)) {
  if (!study_filter %in% studies) {
    stop("study_filter not found in cfg$gwas: ", study_filter)
  }
  studies <- study_filter
}

tic <- Sys.time()

log_msg("=== START GATHER RESULTS ===")
log_msg("Mode: ", mode)
log_msg("Studies: ", paste(studies, collapse = ", "))

# Paths ----

mr_pergene_dir <- file.path(root_data_dir, proj_dir, "04_MREstimation", "perGene")
coloc_pergene_dir <- file.path(root_data_dir, proj_dir, "05_Colocalization", "perGene")

mr_out_dir <- file.path(root_data_dir, proj_dir, "04_MREstimation")
coloc_out_dir <- file.path(root_data_dir, proj_dir, "05_Colocalization")

# Helper functions ----

safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) NULL)
}

strip_suffix <- function(x, suffix_pattern) {
  sub(suffix_pattern, "", basename(x))
}

extract_gene_name <- function(file, study, is, suffix_pattern) {
  x <- strip_suffix(file, suffix_pattern)
  sub(paste0("_", study, "_", is, "$"), "", x)
}

# Gather MR results ----

for (study in studies) {
  
  log_msg("Gathering MR results for ", study)
  
  for (is in c("IS1", "IS2", "IS3")) {
    
    pattern_main <- paste0("_", study, "_", is, "_Main\\.rds$")
    pattern_het  <- paste0("_", study, "_", is, "_Heterogeneity\\.rds$")
    pattern_ple  <- paste0("_", study, "_", is, "_Pleiotropy\\.rds$")
    pattern_loo  <- paste0("_", study, "_", is, "_LOO\\.rds$")
    pattern_snp  <- paste0("_", study, "_", is, "_SNPs\\.rds$")
    pattern_ste  <- paste0("_", study, "_", is, "_Steiger\\.rds$")
    
    files_main <- list.files(mr_pergene_dir, pattern = pattern_main, full.names = TRUE)
    files_het  <- list.files(mr_pergene_dir, pattern = pattern_het,  full.names = TRUE)
    files_ple  <- list.files(mr_pergene_dir, pattern = pattern_ple,  full.names = TRUE)
    files_loo  <- list.files(mr_pergene_dir, pattern = pattern_loo,  full.names = TRUE)
    files_snp  <- list.files(mr_pergene_dir, pattern = pattern_snp,  full.names = TRUE)
    files_ste  <- list.files(mr_pergene_dir, pattern = pattern_ste,  full.names = TRUE)
    
    RES <- list()
    HET <- list()
    PLE <- list()
    LOO <- list()
    SNP <- list()
    STE <- list()
    
    for (f in files_main) {
      gene <- extract_gene_name(f, study, is, "_Main\\.rds$")
      RES[[gene]] <- safe_read_rds(f)
    }
    
    for (f in files_het) {
      gene <- extract_gene_name(f, study, is, "_Heterogeneity\\.rds$")
      HET[[gene]] <- safe_read_rds(f)
    }
    
    for (f in files_ple) {
      gene <- extract_gene_name(f, study, is, "_Pleiotropy\\.rds$")
      PLE[[gene]] <- safe_read_rds(f)
    }
    
    for (f in files_loo) {
      gene <- extract_gene_name(f, study, is, "_LOO\\.rds$")
      LOO[[gene]] <- safe_read_rds(f)
    }
    
    for (f in files_snp) {
      gene <- extract_gene_name(f, study, is, "_SNPs\\.rds$")
      SNP[[gene]] <- safe_read_rds(f)
    }
    
    for (f in files_ste) {
      gene <- extract_gene_name(f, study, is, "_Steiger\\.rds$")
      STE[[gene]] <- safe_read_rds(f)
    }
    
    saveRDS(RES, file = file.path(mr_out_dir, paste0(study, "_", is, "_Main.rds")))
    saveRDS(HET, file = file.path(mr_out_dir, paste0(study, "_", is, "_Heterogeneity.rds")))
    saveRDS(PLE, file = file.path(mr_out_dir, paste0(study, "_", is, "_Pleiotropy.rds")))
    saveRDS(LOO, file = file.path(mr_out_dir, paste0(study, "_", is, "_LOO.rds")))
    saveRDS(SNP, file = file.path(mr_out_dir, paste0(study, "_", is, "_SNPs.rds")))
    saveRDS(STE, file = file.path(mr_out_dir, paste0(study, "_", is, "_Steiger.rds")))
    
    log_msg("Saved MR gathered objects for ", study, " ", is)
  }
}

# Gather colocalization results ----

for (study in studies) {
  
  log_msg("Gathering colocalization results for ", study)
  
  files_coloc <- list.files(
    coloc_pergene_dir,
    pattern = paste0("_", study, "_coloc\\.rds$"),
    full.names = TRUE
  )
  
  files_ldold <- list.files(
    coloc_pergene_dir,
    pattern = paste0("_", study, "_LDcheckOld\\.rds$"),
    full.names = TRUE
  )
  
  files_ldnew <- list.files(
    coloc_pergene_dir,
    pattern = paste0("_", study, "_LDcheck\\.rds$"),
    full.names = TRUE
  )
  
  COL <- list()
  COL2 <- list()
  COL3 <- list()
  
  for (f in files_coloc) {
    gene <- sub(paste0("_", study, "_coloc\\.rds$"), "", basename(f))
    COL[[gene]] <- safe_read_rds(f)
  }
  
  for (f in files_ldold) {
    gene <- sub(paste0("_", study, "_LDcheckOld\\.rds$"), "", basename(f))
    COL2[[gene]] <- safe_read_rds(f)
  }
  
  for (f in files_ldnew) {
    gene <- sub(paste0("_", study, "_LDcheck\\.rds$"), "", basename(f))
    COL3[[gene]] <- safe_read_rds(f)
  }
  
  saveRDS(COL,  file = file.path(coloc_out_dir, paste0(study, "_coloc.rds")))
  saveRDS(COL2, file = file.path(coloc_out_dir, paste0(study, "_LDcheckOld.rds")))
  saveRDS(COL3, file = file.path(coloc_out_dir, paste0(study, "_LDcheck.rds")))
  
  log_msg("Saved colocalization gathered objects for ", study)
}

# Metadata ----

metadata <- list(
  script = "07_gatherResults.R",
  mode = mode,
  study_filter = ifelse(is.na(study_filter), NULL, study_filter),
  timestamp = as.character(Sys.time()),
  n_studies = length(studies),
  R_version = R.version.string
)

write_json(
  metadata,
  path = file.path(root_data_dir, proj_dir, paste0("run_metadata_07", ifelse(is.na(study_filter), "", paste0("_", study_filter)), ".json")),
  pretty = TRUE,
  auto_unbox = TRUE
)

# README ----

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 7,
  titol = "Gather per-gene results",
  autors = "Sergio",
  lloc = "07_gatherResults.R",
  temps = temps,
  descripcio = "Per-gene MR and colocalization results are gathered into study-level objects for downstream scoring.",
  altres = ""
)

log_msg("=== FINISHED GATHER RESULTS ===")
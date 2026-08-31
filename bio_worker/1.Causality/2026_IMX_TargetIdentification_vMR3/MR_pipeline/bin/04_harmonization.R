#!/usr/bin/env Rscript

# Arguments ----

args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "pQTL")
study_filter <- ifelse(length(args) >= 2, args[2], NA_character_)

# Libraries ----
.libPaths(c("/renv", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(stringr)
  library(gwasvcf)
  library(gwasglue)
  library(jsonlite)
})

set.seed(1)

# Load config ----

bioinf_root <- Sys.getenv("BIOINF_ROOT", unset = "~/Bioinformatics")
bioinf_root <- path.expand(bioinf_root)

source(file.path(bioinf_root, "Shared/imidomics/R/config.R"))
source(file.path(bioinf_root, "Shared/imidomics/R/readme_generator.R"))
source(file.path(bioinf_root, "Shared/imidomics/R/mr_pipeline_config.R"))

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir  <- imidomics_config$raw_data_dir

cfg <- load_mr_config(mode)

proj_dir  <- cfg$proj_dir
proj_dir0 <- cfg$proj_dir0

message("Running harmonization for mode: ", mode)
if (!is.na(study_filter)) {
  message("Study filter: ", study_filter)
}

tic <- Sys.time()

# Paths ----

vcf_dir <- file.path(root_data_dir, proj_dir0, "tmp_files")
dbfile  <- file.path(vcf_dir, "mydb.sqlite")

harm_dir <- file.path(root_data_dir, proj_dir, cfg$harmonized)
dir.create(harm_dir, recursive = TRUE, showWarnings = FALSE)

# Tools and reference ----

plink_ref <- file.path(root_data_dir, cfg$reference$plink_ref)
set_bcftools(cfg$tools$bcftools)

# Load instruments ----

INS <- readRDS(file.path(root_data_dir, proj_dir, cfg$instruments))

if (!"rsid" %in% names(INS[[1]])) {
  stop("A column rsid is required in instruments object")
}

# Helper functions ----

extract_outcome_vcf <- function(vcf_path, ins, cfg, dbfile) {
  
  vcf <- query_gwas(
    vcf = vcf_path,
    rsid = ins$rsid,
    proxies = "yes",
    dbfile = dbfile,
    tag_r2 = cfg$ld$tag_r2
  )
  
  out.dat <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, "outcome")
  out.dat <- out.dat[out.dat$beta.outcome != 0, ]
  
  out.dat
}

harmonise_wrapper <- function(ins, out) {
  
  har <- harmonise_data(exposure_dat = ins, outcome_dat = out)
  har <- har[har$palindromic == FALSE | har$ambiguous == FALSE, ]
  
  har
}

# GWAS list ----

gwas_list <- cfg$gwas

if (!is.na(study_filter)) {
  if (!study_filter %in% names(gwas_list)) {
    stop("study_filter not found in cfg$gwas: ", study_filter)
  }
  gwas_list <- gwas_list[study_filter]
}

# Main loop ----

for (g in names(gwas_list)) {
  
  message("Processing GWAS: ", g)
  
  HAR <- vector("list", 3)
  
  for (i in 1:3) {
    
    ins <- INS[[i]]
    
    if (gwas_list[[g]]$type == "vcf") {
      out <- extract_outcome_vcf(
        vcf_path = file.path(root_data_dir, cfg$gwas_dir, gwas_list[[g]]$file),
        ins = ins,
        cfg = cfg,
        dbfile = dbfile
      )
    } else {
      stop("Only type 'vcf' is currently supported")
    }
    
    HAR[[i]] <- harmonise_wrapper(ins, out)
  }
  
  saveRDS(HAR, file = file.path(harm_dir, paste0(g, ".rds")))
}

# Metadata ----

metadata <- list(
  script = "04_harmonization.R",
  mode = mode,
  study_filter = ifelse(is.na(study_filter), NULL, study_filter),
  timestamp = as.character(Sys.time()),
  n_gwas = length(gwas_list),
  gwas_names = names(gwas_list),
  R_version = R.version.string
)

write_json(
  metadata,
  file.path(
    harm_dir,
    ifelse(
      is.na(study_filter),
      "run_metadata_04.json",
      paste0("run_metadata_04_", study_filter, ".json")
    )
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

# README ----

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 4,
  titol = "Outcome data preparation & harmonization",
  autors = "Sergio",
  lloc = "04_harmonization.R",
  temps = temps,
  descripcio = "Selected instruments are extracted from GWAS datasets. Missing SNPs are replaced by LD proxies (r²≥threshold). Harmonisation is performed using TwoSampleMR.",
  altres = ""
)

message("Step 04 completed successfully.")
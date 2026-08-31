#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))

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

# Configuration ----
source("~/.imiducb/config.R")
source("~/.imiducb/mr_pipeline_config.R")
# Loads: root_data_dir, raw_data_dir, volume_dir, analyses_dir, reference_source_dir, tmp_dir
# Sources: readme_generator.R and register_analysis.R

if (!mode %in% c("pQTL", "eQTL")) {
  stop("Unsupported mode: ", mode)
}

cfg <- load_mr_config(mode)
proj_dir <- file.path(analyses_dir, cfg$proj_dir)
proj_dir0 <- file.path(analyses_dir, cfg$proj_dir0)

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

plink_ref <- cfg$reference$plink_ref
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
      vcf_path <- file.path(reference_source_dir, gwas_list[[g]]$file)
      out <- extract_outcome_vcf(
        vcf_path = vcf_path,
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

# Register analysis ----
register_analysis(
  analysis_type = "harmonization",
  analysis_subtype = mode,
  analysis_method = "custom",
  analysis_design = ifelse(is.na(study_filter), "all", paste0("study_filter=", study_filter)),
  parameters = list(mode = mode, study_filter = ifelse(is.na(study_filter), NA_character_, study_filter)),
  datasets = list(),
  reference_sources = as.list(unique(c(unname(unlist(lapply(cfg$gwas, function(x) x$reference_id))), cfg$reference$plink_ref_id))),
  input_files = list(),
  analysis_path = harm_dir,
  analysis_code = script_path,
  status = "done",
  computational_running_time = temps,
  description = "Selected instruments are extracted from GWAS datasets and harmonised with TwoSampleMR for downstream MR estimation."
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
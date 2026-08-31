#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(stringr)
  library(gwasvcf)
  library(gwasglue)
  library(jsonlite)
})

set.seed(1)

# Arguments ----

args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "pQTL")

# Load config ----

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
source("~/Bioinformatics/Shared/imidomics/R/mr_pipeline_config.R")

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir  <- imidomics_config$raw_data_dir

cfg <- load_mr_config(mode)

proj_dir  <- cfg$proj_dir
proj_dir0 <- cfg$proj_dir0
gwas_dir  <- cfg$gwas_dir

message("Running harmonization for mode: ", mode)

tic <- Sys.time()

# Paths ----

vcf_dir <- file.path(root_data_dir, proj_dir0, "tmp_files")
dbfile  <- file.path(vcf_dir, "mydb.sqlite")

harm_dir <- file.path(root_data_dir, proj_dir, cfg$harmonized)
dir.create(harm_dir, recursive = TRUE, showWarnings = FALSE)

# Tools & reference ----

plink_ref <- file.path(root_data_dir, cfg$reference$plink_ref)
set_bcftools(cfg$tools$bcftools)

# Optional: set plink if needed later
# set_plink(cfg$tools$plink_bin)

# Load instruments ----

INS <- readRDS(file.path(root_data_dir, proj_dir, cfg$instruments))

if(!"rsid" %in% names(INS[[1]])){
  stop("A column rsid is required in instruments object")
}

# Helper functions ----

## Extract outcome from VCF (native proxy handling) ----

extract_outcome_vcf <- function(vcf_path, ins, cfg){
  
  vcf <- query_gwas(
    vcf = vcf_path,
    rsid = ins$rsid,
    proxies = "yes",
    dbfile = dbfile,
    tag_r2 = cfg$ld$tag_r2
  )
  
  out.dat <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, "outcome")
  out.dat <- out.dat[out.dat$beta.outcome != 0,]
  
  return(out.dat)
}

## Harmonisation wrapper ----

harmonise_wrapper <- function(ins, out){
  
  har <- harmonise_data(exposure_dat = ins, outcome_dat = out)
  
  # Remove problematic palindromic SNPs
  har <- har[har$palindromic == FALSE | har$ambiguous == FALSE,]
  
  return(har)
}

# GWAS list ----

gwas_list <- cfg$gwas

# Main loop ----

for(g in names(gwas_list)){
  
  message("Processing GWAS: ", g)
  
  HAR <- vector("list", 3)
  
  for(i in 1:3){
    
    ins <- INS[[i]]
    
    if(g$type == "vcf"){
      out <- extract_outcome_vcf(file.path(cfg$gwas_dir, gwas_list[[g]]$file), ins, cfg)
    } else{
      stop("Only type vcf supported")
    }
    
    HAR[[i]] <- harmonise_wrapper(ins, out)
  }
  
  saveRDS(HAR, file = file.path(harm_dir, paste0(g, ".rds")))
}

# Metadata ----

metadata <- list(
  script = "04_harmonization.R",
  mode = mode,
  timestamp = as.character(Sys.time()),
  n_gwas = length(gwas_list),
  R_version = R.version.string
)

write_json(
  metadata,
  file.path(harm_dir, "run_metadata_04.json"),
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
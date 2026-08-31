#!/usr/bin/env Rscript

tic_global <- Sys.time()

# Config ----
source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir0 <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/"
proj_dir <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/"

studies <- c(
  "RA_Ishigaki2022",
  "PSA_Soomro2022",
  # "CD_deLange2017",
  # "UC_deLange2017",
  # "SLE_Bentham2015",
  "PS_Dand2025",
  "SS_Sakaue2021",
  "AD_BuduAggrey2023"
)

# Gather results ----
STE <- COL3 <- COL2 <- COL <- SNP <- HET <- PLE <- LOO <- RES <- list()
for(name in studies){
  STE[[name]] <- COL3[[name]] <- COL2[[name]] <- COL[[name]] <- SNP[[name]] <- HET[[name]] <- PLE[[name]] <- LOO[[name]] <- RES[[name]] <- list()
  
  # All files 05
  all.files <- list.files(paste0(root_data_dir, proj_dir, "/05_Colocalization/perGene/"), full.names = T)
  all.files <- grep(name, all.files, value = T)

  # Coloc results
  files <- grep("_coloc.rds", all.files, value = T)
  for(file in files){
    apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
    COL[[name]][[apt]] <- readRDS(file)
  }
  saveRDS(COL[[name]], paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_coloc.rds"))

  # LDcheck results (old version)
  files <- grep("_LDcheckOld.rds", all.files, value = T)
  for(file in files){
    apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
    COL2[[name]][[apt]] <- readRDS(file)
  }
  saveRDS(COL2[[name]], paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_LDcheckOld.rds"))

  # LDcheck results
  files <- grep("_LDcheck.rds", all.files, value = T)
  for(file in files){
    apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
    COL3[[name]][[apt]] <- readRDS(file)
  }
  saveRDS(COL3[[name]], paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_LDcheck.rds"))

  # All files 04
  all.files <- list.files(paste0(root_data_dir, proj_dir, "/04_MREstimation/perGene/"), full.names = T)
  all.files <- grep(name, all.files, value = T)
  for(is in c("IS1", "IS2", "IS3")){

    # Main results
    files <- grep(is, all.files, value = T)
    files <- grep("_Main.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      RES[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(RES[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Main.rds"))

    # LOO results
    files <- grep(is, all.files, value = T)
    files <- grep("LOO.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      LOO[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(LOO[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_LOO.rds"))

    # PLE results
    files <- grep(is, all.files, value = T)
    files <- grep("_Pleiotropy.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      PLE[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(PLE[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Pleiotropy.rds"))

    # HET results
    files <- grep(is, all.files, value = T)
    files <- grep("_Heterogeneity.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      HET[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(HET[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Heterogeneity.rds"))

    # SNP results
    files <- grep(is, all.files, value = T)
    files <- grep("_SNPs.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      SNP[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(SNP[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_SNP.rds"))

    # SNP results
    files <- grep(is, all.files, value = T)
    files <- grep("_Steiger.rds", files, value = T)
    for(file in files){
      apt <- gsub(paste0("(.*/)|(_", name, ".*)"), "", file)
      STE[[name]][[is]][[apt]] <- readRDS(file)
    }
    saveRDS(STE[[name]][[is]], paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Steiger.rds"))
  } #is
}# name

# README + metadata ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic_global, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 7,
  titol = "Put together per-gene results",
  autors = "Sergio",
  lloc = "07_gatherResults.R",
  temps = temps,
  descripcio = "This step creates an object for each study and each instrument selection set, gathering the results for all genes/proteins created in both the estimatino and colocalization steps.",
  altres = ""
)

message("07 completed successfully.")
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

study_index <- as.numeric(args[1])
is_index    <- as.numeric(args[2])

tic <- Sys.time()

#### Libraries ####
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(jsonlite)
})

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/"
ld_dir <- file.path(root_data_dir, proj_dir, "LD_reference")
source(paste0("~/Bioinformatics/", proj_dir, "/tweaked_functions_estimation.R"))

#### Parameters ####
# Studies ----
studies <- c(
  "RA_Ishigaki2022",
  "PSA_Soomro2022",
  "CD_deLange2017",
  "UC_deLange2017",
  "SLE_Bentham2015",
  "PS_Dand2025",
  "SS_Sakaue2021",
  "AD_BuduAggrey2023"
)

study <- studies[study_index]
is <- paste0("IS",is_index)

message("Running study=", study, " instrument_set=", is)

# Prevalence ----
prevs <- c(1, 0.5, 0.3, 0.5, 0.24, 2, 0.5, 9.3)/100
names(prevs) <- studies

# Load case/sample size ----
ao <- readRDS(file.path(root_data_dir, proj_dir, "opengwas_availableoutcomes.rds"))

cases <- c(
  RA_Ishigaki2022 = 22350,
  PSA_Soomro2022 = 5065,
  CD_deLange2017 = ao$ncase[ao$id=="ebi-a-GCST004132"],
  UC_deLange2017 = ao$ncase[ao$id=="ebi-a-GCST004133"],
  SLE_Bentham2015 = ao$ncase[ao$id=="ebi-a-GCST003156"],
  PS_Dand2025 = 36466,
  SS_Sakaue2021 = 1599,
  AD_BuduAggrey2023 = 60653
)

ss <- c(
  RA_Ishigaki2022 = 22350+74823,
  PSA_Soomro2022 = 5065+21286,
  CD_deLange2017 = ao$sample_size[ao$id=="ebi-a-GCST004132"],
  UC_deLange2017 = ao$sample_size[ao$id=="ebi-a-GCST004133"],
  SLE_Bentham2015 = ao$sample_size[ao$id=="ebi-a-GCST003156"],
  PS_Dand2025 = 36466+458078,
  SS_Sakaue2021 = 659915,
  AD_BuduAggrey2023 = 60653+804329
)

# Load MoE model ----
load(file.path(root_data_dir, proj_dir, "MoE_Model/rf.rdata"))

# Set Plink ----
if(root_data_dir == "/mnt/proteome/bioinformatics"){
  # Binary
  plink_bin <- "/mnt/nfs_home/smartinez/tools/plink/plink"
  # set_plink(path = "/mnt/nfs_home/smartinez/tools/plink/plink")
  Sys.setenv(PATH=paste("/home/sergio/bin/", Sys.getenv("PATH"), sep=":"))
  
  # Reference
  plink_ref <- paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR")
  # sudo mkdir /tmp/ramdisk
  # sudo mount -t tmpfs -o size=1G tmpfs /tmp/ramdisk
  # plink_ref <- "/dev/shm/tmp/EUR"
} else{
  # Pending
  # Sys.setenv(PATH=paste("/home/sergio/bcftools/", Sys.getenv("PATH"), sep=":"))
}

## Reference panel for LD
refPanel <- fread(paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR.bim"))

# Load harmonized data ----
harm_data <- readRDS(file.path(root_data_dir, proj_dir, "03_Harmonized", paste0(study, ".rds")))[[is_index]]
if(nrow(harm_data)==0) quit(save="no")
harm_data$exposure <- gsub(".*-","",harm_data$exposure)
genes <- unique(harm_data$exposure)

# Output folder ----
out_dir <- file.path(root_data_dir,proj_dir,"04_MREstimation/perGene")
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

# Comprovar genes already done
files <- list.files(paste0(root_data_dir, proj_dir, "/04_MREstimation/perGene/"), full.names = T)

# MR estimations ----
for(gene in genes){
  
  # Check pending files
  # is.done <- sum(grepl(gene, files) & grepl(study, files) & grepl(is, files))==6
  # if(is.done) next()
  
  message("Processing gene ", gene)
  
  prefix <- file.path(out_dir, paste(gene, study, is, sep = "_"))
  
  # main_file <- paste0(prefix,"_Main.rds")
  # if(file.exists(main_file)) next
  
  aux <- harm_data[harm_data$exposure==gene,]
  aux <- aux[aux$mr_keep==TRUE,]
  
  if(nrow(aux)==0) next
      
  tic2 <- Sys.time()
      
  # Units and prevalence
  aux$units.outcome <- "log odds"
  aux$prevalence.outcome <- prevs[study]
  
  # The following data must be available
  aux$eaf.outcome <- aux$eaf.exposure
  # aux$eaf.exposure
  
  # Cases
  aux$ncase.outcome <- cases[study]
  
  # Controls
  aux$samplesize.outcome <- ss[study]
  aux$ncontrol.outcome <- aux$samplesize.outcome-aux$ncase.outcome
  
  # When total sample size is lower, compute cases by ratio with respect to max sample size
  idx <- which(aux$samplesize.outcome < max(aux$samplesize.outcome))
  aux$ncase.outcome[idx] <- round(cases[study]/max(aux$samplesize.outcome)*aux$samplesize.outcome[idx])
  
  # Fix controls as well
  aux$ncontrol.outcome[idx] <- round(aux$samplesize.outcome[idx]-aux$ncase.outcome[idx])
  
  # Steiger test
  aux$r.exposure <- get_r_from_bsen(b = aux$beta.exposure, se = aux$se.exposure, n = aux$samplesize.exposure)
  aux$r.outcome <- get_r_from_lor(lor = aux$beta.outcome, af = aux$eaf.outcome, ncase = aux$ncase.outcome, ncontrol = aux$ncontrol.outcome, prevalence = aux$prevalence.outcome)
  suppressWarnings({res.dir <- directionality_test(aux)})
  saveRDS(res.dir, file = paste0(prefix, "_Steiger.rds"))
  
  # Transform to MRInput format
  if("rsid" %in% names(aux)) aux$SNP <- aux$rsid
  df <- dat_to_MRInput(aux, get_correlation=TRUE)
  if(length(df[[1]]@snps)==0) next()
  
  # Analyze with multiple instrument methods (>5 instruments but can be >2)
  res <- lapply(df, function(x) if(length(x@betaX)>2) mr_allmethods(x, method = "main")@Values else NA)
  names(res) <- gsub("\\..*", "", names(res))
  # names(res) <- gsub(".*-", "", names(res))
  
  # Analyze when < 6 instruments
  res2 <- mr(aux)
  res.het <- list()
  res.snps <- list()
  res.egger <- list()
  res.loo <- list()
  for(id in names(res)){
    
    # Get assumption tests
    if(unique(res2$nsnp[res2$exposure==id])>1){
      if(is.null(res.het[[id]])) res.het[[id]] <- mr_heterogeneity(aux[aux$exposure==id,])
      if(is.null(res.egger[[id]])) res.egger[[id]] <- mr_pleiotropy_test(aux[aux$exposure==id,])
      res.loo[[id]] <- mr_leaveoneout(aux[aux$exposure==id,])
      if(nrow(aux[aux$exposure==id,])==2) res.loo[[id]] <- rbind(mr_singlesnp(aux[aux$exposure==id,])[1:2,], res.loo[[id]])
    }
    res.snps[[id]] <- aux$SNP[aux$exposure==id]
    
    # If only 1 or 2 instruments, get estimates
    if(!is.data.frame(res[[id]]) & id %in% res2$exposure){
      res[[id]] <- data.frame(res2[res2$exposure==id,c("method", "b", "se")], LL=res2[res2$exposure==id, "b"]-1.96*res2[res2$exposure==id, "se"], UL=res2[res2$exposure==id, "b"]+1.96*res2[res2$exposure==id, "se"], Pvalue=res2[res2$exposure==id,"pval"])
      names(res[[id]]) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value")
      rownames(res[[id]]) <- NULL
    }
  }
  
  #### Keep heterogeneity statistics as well!!!
  
  # Analyze MR-RAPS & MR-MoE if >5 instruments
  for(id in names(res)){
    if(nrow(aux[aux$exposure==id,]) > 5){#is.data.frame(res[[id]])
      
      #### MR-RAPS ####
      res2 <- mr(dat = aux[aux$exposure==id,], method_list = c("mr_raps"))
      res2$method <- "RAPS"
      aux2 <- data.frame(res2[res2$exposure==id,c("method", "b", "se")], LL=res2[res2$exposure==id, "b"]-1.96*res2[res2$exposure==id, "se"], UL=res2[res2$exposure==id, "b"]+1.96*res2[res2$exposure==id, "se"], Pvalue=res2[res2$exposure==id,"pval"])
      names(aux2) <- names(res[[id]])
      res[[id]] <- rbind(res[[id]], aux2)
      rownames(res[[id]]) <- NULL
      
      #### MR-MoE ####
      res2 <- mr_wrapper(dat = aux[aux$exposure==id,])
      res2 <- mr_moe(res2, rf)
      
      # Get heterogeneity that corresponds to top model 
      est <- as.data.frame(res2[[1]]$estimates)
      het.aux <- as.data.frame(res2[[1]]$heterogeneity)
      het.aux <- het.aux[het.aux$steiger_filtered==est$steiger_filtered[1] & 
                           het.aux$outlier_filtered==est$outlier_filtered[1] &
                           het.aux$method=="IVW",]
      res.het[[id]] <- het.aux 
      
      # Get pleiotropy test that corresponds to top model 
      ple.method <- ifelse(grepl("FE", est$method[1]), "FE Egger intercept", "RE Egger intercept")
      ple.aux <- as.data.frame(res2[[1]]$directional_pleiotropy)
      ple.aux <- ple.aux[ple.aux$steiger_filtered==est$steiger_filtered[1] & 
                           ple.aux$outlier_filtered==est$outlier_filtered[1] &
                           ple.aux$method==ple.method,]
      res.egger[[id]] <- ple.aux 
      
      # Get best MoE estimate
      res2 <- as.data.frame(res2[[1]]$estimates)
      
      # Reformat to add to main results
      aux2 <- data.frame(res2[1,c("method", "b", "se")], LL=res2[1, "ci_low"], UL=res2[1, "ci_upp"], Pvalue=res2[1,"pval"])
      
      # Add the results to the main object
      names(aux2) <- names(res[[id]])
      res[[id]] <- rbind(res[[id]], aux2)
      rownames(res[[id]]) <- NULL
    }
  }
  
  saveRDS(res, file = paste0(prefix, "_Main.rds"))
  saveRDS(res.het, file = paste0(prefix, "_Heterogeneity.rds"))
  saveRDS(res.snps, file = paste0(prefix, "_SNPs.rds"))
  saveRDS(res.egger, file = paste0(prefix, "_Pleiotropy.rds"))
  saveRDS(res.loo, file = paste0(prefix, "_LOO.rds"))
  
  toc2 <- Sys.time()
  print(paste0("Time for ", gene, ": ", as.numeric(difftime(toc2, tic2, units='mins'))))
  print(paste0("Total time so far: ", as.numeric(difftime(toc2, tic, units='mins'))))
  
}# gene

metadata <- list(
  script="05_estimation.R",
  study=study,
  instrument_set=is,
  timestamp=as.character(Sys.time()),
  R_version=R.version.string
)

write_json(
  metadata,
  file.path(root_data_dir, proj_dir, "04_MREstimation", paste0("run_metadata_",study,"_",is,".json")),
  pretty=TRUE,
  auto_unbox=TRUE
)

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units='mins'))
readme(path = paste0(root_data_dir, proj_dir, "/README"),
       subanalisi = 5,
       titol = "MR analysis using the harmonised data",
       autors = "Sergio",
       lloc = "05_estimation.R",
       temps = temps,
       descripcio = "For each disease and each instrument selection method, the MR analysis pipeline is applied. The pipeline considers two cases depending of the number of instrument variables. For more than 5 instruments, 3 methods are applied: MendelianRandomization::mr_allmethods (the main method being random effects IVW accounting for correlation), MR-RAPS and MR-MoE. For 5 or less instruments, the methods included in TwoSampleMR::mr are applied. The results are stored in a per-protein list as /04_MREstimation/perGene/IMID_STUDY_ISX_Main.rds. In both cases, pleiotropy and heterogeneity tests are performed in addition to the MR analysis, and the Leave-One-Out test is performed in the case of <=5 instruments. These 3 sensitivity analyses are stored in per-protein lists with the sufixes _Pleiotropy.rds, _Heterogeneity.rds and _LOO.rds.",
       altres = "")

message("Step 05 finished.")






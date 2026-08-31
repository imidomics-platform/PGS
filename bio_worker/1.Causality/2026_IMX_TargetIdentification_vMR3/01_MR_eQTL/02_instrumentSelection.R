#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(MRInstruments)
library(data.table)
library(stringr)

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir <- "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_MR_eQTL/"
jwtkey <- "" # this needs to be filled according to the jwt api password

tic <- Sys.time()

# Available outcomes from TwoSampleMR
ao <- readRDS(paste0(root_data_dir, "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR/tmp_files/opengwas_availableoutcomes.rds"))

# Available genes in Vosa
ids <- ao$id[ao$author=="Vosa U"]
genes <- gsub(".*-", "", ids)

# Keep only those that have TSS
ensembl.tss <- readRDS(paste0(root_data_dir, proj_dir, "/Annotations/All_ENS_SYMB_TSS_v37.rds"))
genes <- genes[genes %in% ensembl.tss$ensembl_gene_id]
ids <- ids[gsub(".*-", "", ids) %in% genes]

#### Analysis #### 
tic <- Sys.time()

lds <- c(0.001, 0.1, 0.2)
for(i in args[1]){
  i <- as.numeric(i)
  INS <- data.frame()
  failed_genes <- c()
  for(gene in ids){
    
    # Extract significant instruments and clump
    ins <- try({extract_instruments(outcomes = gene, kb = 250000, r2 = lds[i], opengwas_jwt = jwtkey)})
    if(class(ins)=="try-error"){
      failed_genes <- c(failed_genes, gene)
      next()
    }
    
    # Select only cis if possible
    if(!is.null(ins)){
      ins$cis <- NA
      ins$symbol <- NA
      id <- unique(ins$exposure)
      ens <- gsub(".*-", "", id)
      symb <- unique(ensembl.tss$hgnc_symbol[ensembl.tss$ensembl_gene_id==ens])
      ins$symb <- ifelse(length(symb)==1, symb, NA)
      
      # Check that ensID was found in biomart
      if(ens %in% ensembl.tss$ensembl_gene_id){
        cis.chr <- unique(ensembl.tss$chromosome_name[ensembl.tss$ensembl_gene_id == ens])
        
        # Mark the trans SNPs if any
        if(any(ins$chr.exposure!=cis.chr)) ins$cis[which(ins$chr.exposure!=cis.chr)] <- F 
        
        # Mark cis SNPs and, if any, remove the trans
        if(any(ins$chr.exposure==cis.chr)){
          
          # Get TSS as the average position if there are many
          tss <- mean(ensembl.tss$transcription_start_site[ensembl.tss$ensembl_gene_id == ens])
          
          # Classify CIS/TRANS based on position
          ins$cis[which(ins$chr.exposure==cis.chr)] <- abs(ins$pos.exposure[which(ins$chr.exposure==cis.chr)]-tss)<1E6
          
          # If there are both CIS and TRANS, keep only CIS
          if(sum(ins$cis)>0 & sum(!ins$cis)>0){
            ins <- ins[-which(ins$cis==F),]
          }# both cis and trans
        }# any cis
      }# in ensembl.tss
      INS <- rbind(INS, ins)
      toc2 <- Sys.time()
      temps <- as.numeric(difftime(toc2, tic, units='mins'))
      print(paste0(which(ids==gene), " genes in ", temps))
    }# not null
  }# gene
  
  saveRDS(INS, file = paste0(paste0(root_data_dir, proj_dir, "/01_eQTL_Instruments/IS", i, "_full.rds")))
  saveRDS(failed_genes, file = paste0(paste0(root_data_dir, proj_dir, "/01_eQTL_Instruments/failed_genes_IS", i, ".rds")))
}# is
toc <- Sys.time()
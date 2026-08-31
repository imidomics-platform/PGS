#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
chunk <- as.numeric(args[1])
do2 <- as.numeric(args[2])

library(dplyr)
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(coloc)
library(gwasvcf)

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir0 <- "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/"
proj_dir <- "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_MR_pQTL/"
source(paste0("~/Bioinformatics/", proj_dir0, "/HelperFunctions/tweaked_functions_coloc.R"))

#### Set binaries ####
if(root_data_dir == "/mnt/proteome/bioinformatics"){
  Sys.setenv(PATH=paste("/mnt/nfs_home/smartinez/tools/", Sys.getenv("PATH"), sep=":"))
  set_plink(path = "/mnt/nfs_home/smartinez/tools/plink/plink")
  set_bcftools(path = "/mnt/nfs_home/smartinez/tools/bcftools/bcftools")
} else{
  Sys.setenv(PATH=paste("/home/sergio/bcftools/", Sys.getenv("PATH"), sep=":"))
}

tic3 <- Sys.time()

#### Parameters ####
studies <- c("RA_Okada2014", "PSA_Finn2021", "CD_deLange2017", "UC_deLange2017", "SLE_Bentham2015", "PS_Stuart2022", "SS_Sakaue2021", "AD_Oliva2025")
codes <- c("ieu-a-833", "finn-b-L12_PSORI_ARTHRO", "ebi-a-GCST004132", "ebi-a-GCST004133", "ebi-a-GCST003156")
ao <- readRDS(paste0(root_data_dir, proj_dir0, "/tmp_files/opengwas_availableoutcomes.rds"))

ss <- cases <- numeric(length(studies))
names(ss) <- names(cases) <- studies
for(i in 1:length(codes)){
  code <- codes[i]
  cases[i] <- ao$ncase[ao$id==code]
  ss[i] <- ao$sample_size[ao$id==code]
}
ss["PSA_Finn2021"] <- 212242+1637
cases["PSA_Finn2021"] <- 1637
ss["PS_Stuart2022"] <- 15967+28194
cases["PS_Stuart2022"] <- 15967
ss["SS_Sakaue2021"] <- 659915
cases["SS_Sakaue2021"] <- 1599
ss["AD_Oliva2025"] <- 408472+42963
cases["AD_Oliva2025"] <- 42963

#### Data ####
ensembl.tss <- readRDS(paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/refTSS/SOMAscanProts_FinalTSSMapping_lift37.rds"))
ensembl.tss <- ensembl.tss[ensembl.tss$CHR %in% c(1:22, "X", "Y"),]

# refPanel <- fread("/tmp/ramdisk/EUR.bim")
refPanel <- fread(paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR.bim"))
plink_ref <- paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR")

#### Analysis ####
already.done <- list.files(paste0(root_data_dir, proj_dir, "/05_Colocalization/perGene/"), full.names = T)
failed.genes <- RES <- RES2 <- RES3 <- vector(mode = "list", length = length(studies));  names(failed.genes) <- names(RES) <- names(RES2) <- names(RES3) <- studies
tic <- Sys.time()

# Files to read during the analysis
# fn <- tempfile()
fn <- paste0(root_data_dir, proj_dir, "/05_Colocalization/input_temp_", chunk, "_", do2)
vcf_dir <- paste0(root_data_dir, proj_dir0, "/tmp_files/")
dbfile <- paste0(vcf_dir, "dbfile")

# We keep track of the failed aptamers by using the log file, but we need to initialize it to distinguish between runs
write.table(cbind("This is a new run"), file = paste0(root_data_dir, proj_dir, "/05_Colocalization/log", chunk, "_", do2, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F, append = T)

for(name in studies[do2]){
  
  # Read summary stats outcome
  gwas <- readRDS(paste0(root_data_dir, "/Projects/Internal_Projects/2021_ExternalDatasets_Mining/SNPs/01_Features/01_Clinical_Association/01_caseCtrl/", name, "/", gsub("_.*", "", name), ".rds"))
  gwas$pos <- as.numeric(gwas$pos)
  gwas <- gwas[!is.na(gwas$se.outcome),]
  gwas <- gwas[!is.na(gwas$beta.outcome),]
  
  # Proteins to be analyzed
  # Get cis instruments to test for colocalization
  allIns <- readRDS(paste0(root_data_dir, proj_dir, "03_Harmonized/", name, ".rds"))[[1]]
  allIns <- allIns[allIns$Tier %in% c(1,2,4),]
  prots <- unique(allIns$exposure)
  
  # Do it in chunks
  s <- floor(seq(from=1, to=length(prots), length.out=8+1))
  if(!chunk %in% 1:8) stop("chunk argument must be in 1:8")
  if(chunk == 1){
    chosen.prots <- prots[1:s[2]]
  } else{
    chosen.prots <- prots[(s[chunk]+1):s[chunk+1]]
  }
  
  for(prot_Apt in chosen.prots){
    
    # Check pending files
    # is.done <- sum(grepl(prot_Apt, already.done) & grepl("coloc", already.done))==6
    # if(is.done) next()
    
    cat("Starting ", prot_Apt, " ... \n")
    tic2 <- Sys.time()
    
    # Get protein instruments
    ins <- allIns[which(allIns$exposure == prot_Apt),]
    
    if(length(ins$rsid)>0){
      
      # Get exposure data around cis region
      cisC <- unique(ins$cisChromosome)
      U <- unique(ins$TSS)+2E6
      L <- unique(ins$TSS)-2E6
      system(paste0("awk -F ',' '{ if(NR==1 || ($1 == ", cisC, " && $2 < ", U, " && $2 > ", L, ")) { print }}' ", root_data_dir, proj_dir, "/01_pQTLs_Preprocessed/01_MetaAnalyzed_pQTLs_Clean/", prot_Apt, ".txt > ", fn))
      exposure_dat <- fread(fn, data.table = F)
      exposure_dat[c("Effect", "StdErr", "Freq1", "SampleSize", "P-value")] <- lapply(exposure_dat[c("Effect", "StdErr", "Freq1", "SampleSize", "P-value")], as.numeric)
      
      if(nrow(exposure_dat)<500) print(paste0(prot_Apt, " has less than 500 cis pQTLs"))
      for(i in 1:length(ins$rsid)){
        
        # Get location
        rs <- ins$rsid[i]
        chr <- ins$chr.exposure[i]
        pos <- ins$pos.exposure[i]
        chrompos <- paste0(chr, ":", pos)
        
        # Get the SNPs in the region from the LD matrix
        refPanel2 <- refPanel[which(refPanel$V1==chr),]
        idx <- which(abs(refPanel2$V4 - pos)<25E4)
        rsids <- refPanel2$V2[idx]
        chrompos <- paste0(refPanel2$V1[idx], ":", refPanel2$V4[idx])
        
        # Get exposure data on the region around the snp
        exp.dat <- exposure_dat[match(chrompos, exposure_dat$MarkerName),]
        exp.dat$rsid <- rsids
        exp.dat <- exp.dat[which(!is.na(exp.dat$Chromosome)),]
        if(nrow(exp.dat)<400){
          next()
        }
        
        # Format exposure data
        exp.dat <- format_data(dat = exp.dat,
                               snp_col = "rsid",
                               beta_col = "Effect",
                               se_col = "StdErr",
                               eaf_col = "Freq1",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               samplesize_col = "SampleSize",
                               chr_col = "Chromosome",
                               pos_col = "Position",
                               log_pval = F, 
                               min_pval = 1e-320
        )
        exp.dat$exposure <- prot_Apt
        
        # Get outcome data
        outcome_dat <- gwas[gwas$SNP %in% exp.dat$SNP,]
        outcome_dat <- outcome_dat[!duplicated(outcome_dat$SNP),]
        
        # Harmonise the exposure and outcome data
        tic <- Sys.time()
        dat <- harmonise_data(exp.dat, outcome_dat)
        toc <- Sys.time()
        print(as.numeric(difftime(toc, tic, units='mins')))
        
        # Make sure all statistics are available
        # if(is.na(dat$eaf.outcome[1])) dat$eaf.outcome <- dat$eaf.exposure
        dat <- dat[which(!is.na(dat$eaf.exposure)),]
        if(!rs %in% dat$SNP) next()
        if(is.null(dat$ncase.outcome)) dat$ncase.outcome <- cases[name]
        if(is.na(dat$ncase.outcome[1])) dat$ncase.outcome <- cases[name]
        if(is.null(dat$ncontrol.outcome)) dat$ncontrol.outcome <- dat$samplesize.outcome-dat$ncase.outcome
        if(is.na(dat$ncontrol.outcome[1])) dat$ncontrol.outcome <- dat$samplesize.outcome-dat$ncase.outcome
        idx <- which(dat$samplesize.outcome < max(dat$samplesize.outcome))
        if(length(idx)>0){
          dat$ncase.outcome[idx] <- round(cases[name]/max(dat$samplesize.outcome)*dat$samplesize.outcome[idx])
          dat$ncontrol.outcome[idx] <- round(dat$samplesize.outcome[idx]-dat$ncase.outcome[idx])
        }
        
        # Get LD matrix
        tic <- Sys.time()
        # print(dim(dat))
        # print(dat$SNP[1:5])
        # saveRDS(dat, file = paste0(root_data_dir, proj_dir, "/test", chunk, ".rds"))
        M <- try({my_ld_matrix_local(dat$SNP, refPanel2, 
                                     plink_bin = "/mnt/nfs_home/smartinez/tools/plink/plink",
                                     bfile = plink_ref)})
        rownames(M) <- gsub("_.*", "", rownames(M))
        colnames(M) <- gsub("_.*", "", colnames(M))
        if("try-error" %in% class(M)){
          cat(M[1], "\n")
          next()
        }
        if(is.null(M)){
          cat("M is null")
          next()
        }
        toc <- Sys.time()
        print(as.numeric(difftime(toc, tic, units='mins')))
        
        # Harmonise LD matrix
        hdat <- dat[match(gsub("_.*", "", colnames(M)), dat$SNP),]
        # tic <- Sys.time()
        # hdat <- harmonise_ld_dat(hdat, M)
        # toc <- Sys.time()
        # print(as.numeric(difftime(toc, tic, units='mins')))
        
        # Check minimum number of variants for coloc
        # if(nrow(hdat$ld) < 350){
        #   write.table(cbind(prot_Apt, " too few variants after harmonization"), file = paste0(root_data_dir, proj_dir, "/05_Colocalization/log", chunk, "_", do2, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F, append = T)
        #   next()
        # }
        
        # Check there are enough variants
        # if(all(class(hdat$ld)=="numeric")) next()

        # Exposure data
        exp <- list()
        exp[["beta"]] <- hdat$beta.exposure
        exp[["varbeta"]] <- hdat$se.exposure^2
        exp[["snp"]] <- hdat$SNP
        exp[["position"]] <- hdat$pos.exposure
        exp[["type"]] <- "quant"
        exp[["N"]] <- hdat$samplesize.exposure
        exp[["MAF"]] <- hdat$eaf.exposure
        exp[["MAF"]][exp[["MAF"]]>0.5] <- 1-exp[["MAF"]][exp[["MAF"]]>0.5]
        
        # Outcome data
        out <- list()
        out[["beta"]] <- hdat$beta.outcome
        out[["varbeta"]] <- hdat$se.outcome^2
        out[["snp"]] <- hdat$SNP
        out[["position"]] <- hdat$pos.outcome
        out[["type"]] <- "cc"
        out[["N"]] <- hdat$samplesize.outcome
        out[["s"]] <- hdat$ncase.outcome/out[["N"]]
        
        # Both have the same LD matrix
        exp$LD <- out$LD <- hdat$ld
        
        # Compute coloc
        tic <- Sys.time()
        RES[[name]][[prot_Apt]][[rs]] <- my.coloc.abf(dataset1=exp, dataset2=out)
        toc <- Sys.time()
        print(as.numeric(difftime(toc, tic, units='mins')))
        
        # LD check (old version)
        tic <- Sys.time()
        top.exposure <- rs
        top.outcome <- hdat$SNP[order(hdat$pval.outcome, decreasing=F)][1:30]
        if(top.exposure %in% rownames(M)){
          check <- any(M[top.exposure, top.outcome]>0.8)
        } else{
          check <- F
        }
        RES2[[name]][[prot_Apt]][[rs]] <- check
        
        # LD check
        tic <- Sys.time()
        top.exposure <- rs
        tops.outcome <- hdat$SNP[order(hdat$pval.outcome, decreasing=F)][1:30]
        top.outcome <- tops.outcome[1]
        if(top.exposure %in% rownames(M)){
          check <- M[top.exposure, top.outcome]^2>0.8 | any(top.exposure %in% tops.outcome)
        } else{
          check <- F
        }
        RES3[[name]][[prot_Apt]][[rs]] <- check
        toc <- Sys.time()
        print(as.numeric(difftime(toc, tic, units='mins')))
        
      }#snp
    }#there are instruments
    
    # Write results
    if(!is.null(RES[[name]][[prot_Apt]])){
      
      saveRDS(RES[[name]][[prot_Apt]], file = paste0(root_data_dir, proj_dir, "/05_Colocalization/perGene/", prot_Apt, "_", name, "_coloc.rds"))
      saveRDS(RES2[[name]][[prot_Apt]], file = paste0(root_data_dir, proj_dir, "/05_Colocalization/perGene/", prot_Apt, "_", name, "_LDcheckOld.rds"))
      saveRDS(RES3[[name]][[prot_Apt]], file = paste0(root_data_dir, proj_dir, "/05_Colocalization/perGene/", prot_Apt, "_", name, "_LDcheck.rds"))
      
    } else{ #write if not ok
      write.table(cbind(prot_Apt, " " , name,  " failed for all instruments"), file = paste0(root_data_dir, proj_dir, "/05_Colocalization/log", chunk, "_", do2, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F, append = T)        
      # saveRDS(failed.genes, file = paste0(root_data_dir, proj_dir, "/05_Colocalization/failed_genes.rds"))
    }# write if not ok
    
    # Time for prot (all diseases)
    toc2 <- Sys.time()
    temps <- as.numeric(difftime(toc2, tic2, units='mins'))
    print(paste0(prot_Apt, " took ", temps))
    toc3 <- Sys.time()
    temps <- as.numeric(difftime(toc3, tic3, units='mins'))
    print(paste0(which(chosen.prots==prot_Apt), " in ", temps))
  }#prot
}#study

print("\nEnded")

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units='mins'))
readme(path = paste0(root_data_dir, proj_dir, "/README"),
       subanalisi = 6,
       titol = "Colocalization analysis",
       autors = "Sergio",
       lloc = "06_colocalization.R",
       temps = temps,
       descripcio = "For each disease and each aptamer, the colocalization is applied. The pipeline considers ... The results are stored in a per-protein list as /05_Colocalization/perGene/aptName_IMID_STUDY_coloc.rds. Also the alternative methods as _SuSiEcoloc.rds and _LDcheck.rds",
       altres = "")
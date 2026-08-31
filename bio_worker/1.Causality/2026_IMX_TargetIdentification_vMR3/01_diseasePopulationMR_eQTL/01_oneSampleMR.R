#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
source("~/Bioinformatics/Shared/imidomics/R/CFMR.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir <- "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_diseasePopulationMR_eQTL/"

library(data.table)
library(ivreg)
tic <- Sys.time()

if(root_data_dir == "/mnt/proteome/bioinformatics"){
  Sys.setenv(R_ENVIRON_USER = "~/.Renviron")
  Sys.setenv(PATH=paste("/mnt/nfs_home/smartinez/tools/bcftools/", Sys.getenv("PATH"), sep=":"))
  plink_bin <- "/mnt/nfs_home/smartinez/tools/plink/plink"
  plink_ref <- "/dev/shm/tmp_files/EUR"
  refPanel <- fread("/dev/shm/tmp_files/EUR.bim")
} else{
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/home/sergio/bin/", sep=":"))
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/home/sergio/bcftools/", sep=":"))
  plink_ref <- "/tmp/ramdisk/"
  refPanel <- fread("/tmp/ramdisk/EUR.bim")
}

#### Parameters ####
diseases <- c("RA", "PS", "PSA", "CD", "UC", "SLE")
da.vars <- c("AR_H_02", "PS_G_0302", "AP_I_020307", "EC_F_02", "CU_F_13", "LES_F_12")
names(da.vars) <- diseases

in.path0 <- paste0(raw_data_dir, "/RAW_S3_Bucket_Synchronized/2022_IMXGRR_CSP4004_IA/03_Imputed/VCF/")
in.path <- paste0(root_data_dir, "/Data/Molecular/Genotype/2022_IMXGRR_CSP4004_IA/01_Preprocessed/pQTL_BCF/")

#### Data ####
# Gene expression
m1 <- readRDS(paste0(root_data_dir, "/Data/Molecular/GeneExpression/WholeBlood/2024_IMX_P4V4_NovaSeqKPE_IA/02_Preprocessed/Batch_Adjusted_VST/Global/Adjusted_VST.rds"))

# Expression metadata
eSampleTable <- readRDS(paste0(root_data_dir, "/Data/Molecular/GeneExpression/WholeBlood/2024_IMX_P4V4_NovaSeqKPE_IA/Metadata/sampleTable.rds"))
eFeatureTable <- readRDS(paste0(root_data_dir, "/Data/Molecular/GeneExpression/WholeBlood/2021_IMX_002_NovaSeqKPE_IA/Metadata/featureTable_afterQC.rds"))

# Processed BD
DC <- readRDS(paste0(root_data_dir, "/Data/Clinical/2008_IMIDKIT_IA/04_subsetPreprocessed/IAP4V4/clinicalData_forPredictors.rds"))
DC$imid <- DC$pathologyId 
levels(DC$imid)[3] <- "CTRL"
rownames(DC) <- DC$indivId
DF <- readRDS(paste0(root_data_dir, "/Data/Clinical/2008_IMIDKIT_IA/04_subsetPreprocessed/IAP4V4/varInfo_forPredictors.rds"))

# Genotype metadata
gSampleTable <- readRDS(paste0(root_data_dir, "/Data/Molecular/Genotype/2022_IMXGRR_CSP4004_IA/Metadata/sampleTable_afterQC.rds"))

# Genotype components
g.pca.scores <- readRDS(paste0(root_data_dir, "/Data/Molecular/Genotype/2022_IMXGRR_CSP4004_IA/DerivedFeatures/Quantitative/Principal_Components/Results/Genetic_PCs.rds"))
rownames(g.pca.scores) <- g.pca.scores$phenotypeId

### IDs ####
# Intersect of expression and genetic samples
ids <- intersect(rownames(eSampleTable), rownames(g.pca.scores))

# Dictionary for genotype id translation
dic <- readRDS(paste0(raw_data_dir, "/RAW_S3_Bucket_Synchronized/2022_IMXGRR_CSP4004_IA/03_Imputed/VCF/Sample_Corresp.rds"))

# Expressed genes
Fids <- eFeatureTable$featureID[eFeatureTable$ExpressedPercentageCategory=="Expressed"]

# Auxiliar dataframe with all data
aux.df <- data.frame(eSampleTable[ids,], DC[ids,], row.names = ids)
# aux.df[ids,paste0("E_PC", 1:60)] <- e.pca.scores[match(ids, rownames(e.pca.scores)),1:60]
aux.df[ids,paste0("G_PC", 1:10)] <- g.pca.scores[match(ids, rownames(g.pca.scores)),1:10]

#### Analysis ####
cis <- fread(paste0(root_data_dir, "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_diseasePopulationMR_eQTL/2019-12-11-cis-eQTLsFDR.txt.gz"))

RES <- lapply(diseases, function(x) list()); names(RES) <- diseases

genes <- unique(cis$Gene)

# Do it in chunks
s <- floor(seq(from=1, to=length(genes), length.out=20+1))
chunk <- as.numeric(args[2])
if(!chunk %in% 1:20) stop("chunk argument must be in 1:20")
if(chunk == 1){
  chosen.genes <- genes[1:s[2]]
} else{
  chosen.genes <- genes[(s[chunk]+1):s[chunk+1]]
}

for(geneid in chosen.genes){
  gc()
  tic2 <- Sys.time()
  
  # Gene info
  # gene <- gsub(".*a-", "", geneid)
  gene <- geneid
  # if(!symb %in% ensembl.tss$hgnc_symbol) next()
  chr.gene <- unique(cis$SNPChr[which(cis$Gene==geneid)])
  if(chr.gene %in% c("X", "Y")) next()
  
  # Reference results for the protein
  fn <- tempfile()
  svars <- cis[which(cis$Gene==gene), c("SNPChr", "SNPPos")]
  if(nrow(svars)<5) next() 
  fwrite(svars, file = fn, quote = F, row.names = F, col.names = F, sep = "\t")
  cat("ok1")
  
  # Extract significant SNPs from chr1 file and read in with R
  # system(paste0("bcftools index ", in.path, "chr", chr.prot, ".dose.bcf"))
  # system(paste0("bcftools view -R ", fn, " -Oz -o ", fn, " ", in.path, "chr", chr.gene, ".dose.bcf"))
  system(paste0("bcftools view -R ", fn, " -Oz -o ", fn, " ", in.path0, "chr", chr.gene, ".dose.vcf.gz"))
  cat("ok2")
  a <- vcfR::read.vcfR(file = fn)
  cat("ok3")
  b <- t(vcfR::extract.gt(x = a, element = "DS", as.numeric = T))
  cat("ok4")
  b <- as.data.frame(b)
  rownames(b) <- dic$phenotypeId[match(gsub("_Sample.*", "", rownames(b)), dic$V2)]
  # rownames(b) <- gsub("0_", "", rownames(b))
  colnames(b) <- gsub(":", "_", colnames(b))
  gvars <- make.names(colnames(b))
  cat("ok5")
  
  if(!gene %in% rownames(m1)) next()
  if(ncol(b)<5) next()
  cat("ok6")
  
  # Get names in the intersection of genotype and expression data
  ids <- intersect(colnames(m1), rownames(b))
  
  for(dis in diseases[as.numeric(args[1])]){
    
    # Merge data from genotype, protein and clinical
    df <- data.frame(b[ids,, drop=F], Gene=t(m1)[ids, gene], DA=DC[ids, da.vars[dis]])
    df <- df[!is.na(df$DA),]
    cat("ok7")
    
    # CFI estimation
    # CFI <- try(suppressWarnings({build_IV_sub_sample(k=5, G=df[,gvars,drop=F], X=df[,"Gene"])}), silent = T)
    
    # Check if an error occurred
    CFI <- tryCatch({
      build_IV_sub_sample(k=5, G=df[,gvars,drop=F], X=df[,"Gene"])
    },
    warning = function(w) {
      if (grepl("Convergence for .* lambda value not reached", w$message)) {
        cat("Warning: Convergence issue detected, continuing execution.\n")
        invokeRestart("muffleWarning")
      }
      return(NULL)  # You can return NULL or take other action if you want
    },
    error = function(e) {
      cat("An error occurred: ", e$message, "\n")
      return(NULL)
    })
    if (!is.null(CFI)) {
      print("Model fitted successfully despite warnings.")
    } else {
      print("An error occurred or convergence was not reached.")
      next()
    }
    
    # CFI <- build_IV_sub_sample(k=5, G=df[,gvars,drop=F], X=df[,"Gene"])
    
    df[rownames(CFI), "CFI"] <- CFI
    
    #cross fitted IV regression
    Y <- as.numeric(df$DA)
    X <- df$Gene
    CFI <- df$CFI
    res <- ivreg(Y~X|CFI)
    RES[[dis]][[gene]] <- res
    cat("ok8")
  }# dis
  
  cat("Gene ", which(unique(cis$Gene)==gene), " done\n")
  toc2 <- Sys.time()
  print(paste0("Time for ", gene, ": ", as.numeric(difftime(toc2, tic2, units='mins'))))
  print(paste0("Total time so far: ", as.numeric(difftime(toc2, tic, units='mins'))))
}# prot
# for(dis in diseases) saveRDS(RES[[dis]], file = paste0(root_data_dir, "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_diseasePopulationMR_eQTL/01_OneSampleMR_CFMR/", dis, ".rds"))
saveRDS(RES[[dis]], file = paste0(root_data_dir, "/Projects/Internal_Projects/2021_IMX_TargetIdentification_vMR2/01_diseasePopulationMR_eQTL/01_OneSampleMR_CFMR/", dis, "_chunk", chunk, ".rds"))

cat("Finished successfully")
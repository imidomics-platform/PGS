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

# Load config ----
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
source("~/Bioinformatics/Shared/imidomics/R/config.R")

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir

proj_dir0 <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/"
proj_dir  <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/"
gwas_dir  <- "/Projects/Internal_Projects/2021_ExternalDatasets_Mining/SNPs/01_Features/01_Clinical_Association/01_caseCtrl/"

vcf_dir <- file.path(root_data_dir, proj_dir0, "tmp_files")
dbfile  <- file.path(vcf_dir, "mydb.sqlite")

tic <- Sys.time()

# Set Plink ----
if(root_data_dir == "/mnt/proteome/bioinformatics"){
  # Binary
  plink_bin <- "/mnt/nfs_home/smartinez/tools/plink/plink"
  set_bcftools("/mnt/nfs_home/smartinez/tools/bcftools/bcftools")
  # set_plink(path = "/mnt/nfs_home/smartinez/tools/plink/plink")
  Sys.setenv(PATH=paste("/home/sergio/bin/", Sys.getenv("PATH"), sep=":"))
  
  # Reference
  plink_ref <- paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR")
  # sudo mkdir /tmp/ramdisk
  # sudo mount -t tmpfs -o size=1G tmpfs /tmp/ramdisk
  # plink_ref <- "/dev/shm/tmp/EUR"
} else{
  # Pending
}

# Load instruments ----
INS <- readRDS(file.path(root_data_dir, proj_dir, "02_Instruments/ALL.rds"))

# PLINK reference ----
plink_ref <- file.path(root_data_dir,"Data/Molecular/Biological_Annotations/1000G/EUR")

# Helper functions ----
extract_outcome_rds <- function(out.dat, ins){
  
  ins$SNP <- ins$rsid
  
  out <- out.dat[out.dat$SNP %in% ins$SNP,]
  out <- out[!duplicated(out$SNP),]
  
  if(nrow(out) < nrow(ins)){
    
    rsids <- ins$SNP[!ins$SNP %in% out$SNP]
    
    aux <- gwasvcf::get_ld_proxies(
      unique(rsids),
      bfile = plink_ref,
      tag_r2 = 0.8,
      threads = 5
    )
    
    aux$prox.snp <- aux$SNP_B
    aux <- aux[order(aux$R, decreasing=TRUE),]
    
    aux <- aux[aux$prox.snp %in% out.dat$SNP,]
    aux <- aux[!duplicated(aux$SNP_A),]
    
    if(nrow(aux)>0){
      
      newlines <- out.dat[match(aux$prox.snp,out.dat$SNP),]
      
      A1matchesB1 <- newlines$effect_allele.outcome == aux$B1
      A1matchesB2 <- !A1matchesB1
      
      newlines$effect_allele.outcome[A1matchesB1] <- aux$A1[A1matchesB1]
      newlines$other_allele.outcome[A1matchesB1]  <- aux$A2[A1matchesB1]
      
      newlines$effect_allele.outcome[A1matchesB2] <- aux$A2[A1matchesB2]
      newlines$other_allele.outcome[A1matchesB2]  <- aux$A1[A1matchesB2]
      
      newlines$proxy.outcome <- TRUE
      newlines$target_snp.outcome <- aux$SNP_A
      newlines$proxy_snp.outcome  <- aux$SNP_B
      
      add.names <- setdiff(names(newlines),names(out))
      out[,add.names] <- NA
      
      out <- rbind(out,newlines)
    }
  }
  
  return(out)
}

extract_outcome_vcf <- function(vcf_path, ins){
  
  ins$SNP <- ins$rsid
  
  vcf <- query_gwas(
    vcf = vcf_path,
    rsid = ins$SNP,
    proxies = "yes",
    dbfile = dbfile,
    tag_r2 = 0.8
  )
  
  out.dat <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, "outcome")
  out.dat <- out.dat[out.dat$beta.outcome!=0,]
  
  return(out.dat)
}

harmonise_wrapper <- function(ins, out){
  
  har <- harmonise_data(exposure_dat=ins, outcome_dat=out)
  
  har <- har[har$palindromic==FALSE | har$ambiguous==FALSE,]
  
  return(har)
}

# GWAS list ----
gwas_list <- list(
  
  list(name="RA_Ishigaki2022",
       type="vcf",
       path=file.path(root_data_dir,gwas_dir,"RA_Ishigaki2022/RA.vcf.bgz")),
  
  # list(name="RA_Okada2014",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"RA_Okada2014/ieu-a-833.vcf.gz")),
   
  # list(name="PS_Dand2025",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"PS_Dand2025/PS.vcf.bgz")),
  
  # list(name="PS_Stuart2022",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"PS_Stuart2022/ebi-a-GCST90019017.vcf.gz")),
  
  # list(name="CD_deLange2017",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"CD_deLange2017/ebi-a-GCST004132.vcf.gz")),
  
  # list(name="UC_deLange2017",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"UC_deLange2017/ebi-a-GCST004133.vcf.gz")),
  
  # list(name="SLE_Bentham2015",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"SLE_Bentham2015/ebi-a-GCST003156.vcf.gz")),
  
  # list(name="PSA_Soomro2022",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"PSA_Soomro2022/PSA.vcf.bgz")),
  
  # list(name="PSA_Finn2021",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"PSA_Finn2021/finn-b-L12_PSORI_ARTHRO.vcf.gz")),
  
  # list(name="AD_BuduAggrey2023",
  #      type="vcf",
  #      path=file.path(root_data_dir,gwas_dir,"AD_BuduAggrey2023/AD.vcf.bgz"))
  
  # list(name="AD_Oliva2025",
  #      type="rds",
  #      path=file.path(root_data_dir,gwas_dir,"AD_Oliva2025/AD.rds")),
  
  list(name="SS_Sakaue2021",
       type="vcf",
       path=file.path(root_data_dir,gwas_dir,"SS_Sakaue2021/SS.vcf.bgz"))
)

# Main loop ----
for(g in gwas_list){
  
  message("Processing ",g$name)
  HAR <- vector("list", 3)
  if(g$type=="rds"){
    out.dat <- readRDS(g$path)
  }

  for(i in 1:3){
    ins <- INS[[i]]
    ins$SNP <- ins$rsid
    
    if(g$type=="rds"){
      out <- extract_outcome_rds(out.dat, ins)
    }else{
      out <- extract_outcome_vcf(g$path, ins)
    }
    HAR[[i]] <- harmonise_wrapper(ins, out)
  }
  
  saveRDS(HAR, file=file.path(root_data_dir, proj_dir, "03_Harmonized", paste0(g$name,".rds")))
}

# README ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc,tic,units="mins"))

readme(
  path=file.path(root_data_dir,proj_dir,"README"),
  subanalisi=4,
  titol="Outcome data preparation & harmonization",
  autors="Sergio",
  lloc="04_harmonization.R",
  temps=temps,
  descripcio="Selected instruments are extracted from IMID GWAS datasets. If the exact SNP is absent, LD proxies (r²≥0.8) from 1000G EUR are used. Harmonisation is performed using TwoSampleMR.",
  altres=""
)

metadata <- list(
  script="04_harmonization.R",
  timestamp=as.character(Sys.time()),
  n_gwas=length(gwas_list),
  R_version=R.version.string
)

write_json(
  metadata,
  file.path(root_data_dir,proj_dir,"03_Harmonized/run_metadata_04.json"),
  pretty=TRUE,
  auto_unbox=TRUE
)

message("Step 04 completed successfully.")
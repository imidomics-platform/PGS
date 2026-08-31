#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(args[1])
study_index <- as.integer(args[2])

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

log_msg("=== START COLOCALIZATION ===")

if (length(chunk) != 1 || is.na(chunk) || !chunk %in% 1:8) {
  stop("chunk must be an integer in 1:8")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(coloc)
  library(jsonlite)
  library(gwasvcf)
  library(TwoSampleMR)
})

set.seed(1)

# Config ----
source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")
imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir
proj_dir0 <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/"
proj_dir <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/"
source(paste0("~/Bioinformatics/", proj_dir, "/tweaked_functions_colocalization.R"))

tic_global <- Sys.time()

# Binaries ----
if (root_data_dir == "/mnt/proteome/bioinformatics") {
  Sys.setenv(PATH = paste("/mnt/nfs_home/smartinez/tools/", Sys.getenv("PATH"), sep = ":"))
  set_plink(path = "/mnt/nfs_home/smartinez/tools/plink/plink")
  set_bcftools(path = "/mnt/nfs_home/smartinez/tools/bcftools/bcftools")
  plink_bin <- "/mnt/nfs_home/smartinez/tools/plink/plink"
} else {
  Sys.setenv(PATH = paste("/home/sergio/bcftools/", Sys.getenv("PATH"), sep = ":"))
  plink_bin <- "/home/sergio/bin/plink"
}

# Parameters ----
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

if (length(study_index) != 1 || is.na(study_index) || !study_index %in% seq_along(studies)) {
  stop("study_index must be between 1 and ", length(studies))
}

study <- studies[study_index]

log_msg("Chunk: ", chunk, " | Study index: ", study_index, " (", study, ")")

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

# Inputs ----
refPanel <- fread(file.path(root_data_dir, "Data/Molecular/Biological_Annotations/1000G/EUR.bim"), header = FALSE)

plink_ref <- file.path(root_data_dir, "Data/Molecular/Biological_Annotations/1000G/EUR")

cis_dir <- file.path(root_data_dir,proj_dir,"01_pQTLs_Preprocessed/03_CisRegions")

harm_file <- file.path(root_data_dir, proj_dir, "03_Harmonized", paste0(study, ".rds"))

coloc_dir <- file.path(root_data_dir, proj_dir, "05_Colocalization/perGene")
dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)

fn <- file.path(root_data_dir, proj_dir, "05_Colocalization", paste0("input_temp_", chunk, "_", study_index))


# Load outcome GWAS once ----
gwas_path <- file.path(root_data_dir, "Projects/Internal_Projects/2021_ExternalDatasets_Mining/SNPs/01_Features/01_Clinical_Association/01_caseCtrl/", study, paste0(gsub("_.*", "", study), ".rds"))

if(!file.exists(gwas_path)){
  stop("Missing GWAS file: ", gwas_path)
}

gwas <- readRDS(gwas_path)
setDT(gwas)

log_msg("Loaded GWAS: ", study, " | N SNPs: ", nrow(gwas))

if("pos.outcome" %in% names(gwas) && !"pos" %in% names(gwas)){
  gwas[, pos := as.numeric(pos.outcome)]
} else if ("pos" %in% names(gwas)) {
  gwas[, pos := as.numeric(pos)]
}

gwas <- gwas[!is.na(se.outcome) & !is.na(beta.outcome)]
setkey(gwas, SNP)

# Load cis instruments (IS1) ----
allIns <- readRDS(harm_file)[[1]]
allIns <- allIns[allIns$Tier %in% c(1, 2, 4), ]

if(nrow(allIns) == 0){
  quit(save = "no", status = 0)
}

prots <- unique(allIns$exposure)
chunks <- split(prots, cut(seq_along(prots), 8, labels = FALSE))
chosen.prots <- chunks[[chunk]]

if(length(chosen.prots) == 0){
  quit(save = "no", status = 0)
}

# Main loop ----
for(prot_Apt in chosen.prots){

  tic_prot <- Sys.time()
  log_msg("---- Protein ", prot_Apt, " ----")
  log_msg("Index: ", which(chosen.prots == prot_Apt), "/", length(chosen.prots))
  
  coloc_outfile <- file.path(coloc_dir, paste0(prot_Apt, "_", study, "_coloc.rds"))
  ldold_outfile <- file.path(coloc_dir, paste0(prot_Apt, "_", study, "_LDcheckOld.rds"))
  ldnew_outfile <- file.path(coloc_dir, paste0(prot_Apt, "_", study, "_LDcheck.rds"))
  
  # restartable
  if(file.exists(coloc_outfile) && file.exists(ldold_outfile) && file.exists(ldnew_outfile)){
    next()
  }
  
  prot_res <- list()
  prot_ldold <- list()
  prot_ldnew <- list()
  
  ins <- allIns[allIns$exposure == prot_Apt, ]
  if (nrow(ins) == 0 || all(is.na(ins$rsid))) {
    log_msg(prot_Apt, " has no usable rsids")
    next()
  }
  
  # Get exposure data around cis region
  cisC <- unique(ins$cisChromosome)
  U <- unique(ins$TSS)+2E6
  L <- unique(ins$TSS)-2E6
  system(paste0("awk -F ',' '{ if(NR==1 || ($1 == ", cisC, " && $2 < ", U, " && $2 > ", L, ")) { print }}' ", root_data_dir, proj_dir, "/01_pQTLs_Preprocessed/01_MetaAnalyzed_pQTLs_Clean/", prot_Apt, ".txt > ", fn))
  exposure_dat <- fread(fn, data.table = F)
  exposure_dat[c("Effect", "StdErr", "Freq1", "SampleSize", "P-value")] <- lapply(exposure_dat[c("Effect", "StdErr", "Freq1", "SampleSize", "P-value")], as.numeric)
  
  setDT(exposure_dat)
  
  log_msg("Loaded cis-region data...")
  log_msg("Cis variants: ", nrow(exposure_dat))
  
  if (nrow(exposure_dat) < 500) {
    log_msg(prot_Apt, " has less than 500 cis pQTLs")
  }
  
  for(i in seq_len(nrow(ins))){
    
    rs <- ins$rsid[i]
    chr <- ins$chr.exposure[i]
    pos <- ins$pos.exposure[i]
    
    log_msg("  SNP: ", rs, " | chr: ", chr, " pos: ", pos)
    
    if (is.na(rs) || is.na(chr) || is.na(pos)){
      log_msg(prot_Apt, " instrument number ", i, ", has missing rs/chr/pos")
      next
    }
    
    # Variants in +/-250kb around instrument from reference panel
    refPanel2 <- refPanel[V1 == chr]
    idx <- which(abs(refPanel2$V4 - pos) < 25e4)
    if(length(idx) < 400){
      log_msg(prot_Apt, " instrument number ", i, ", has less than 100 variants in the 250k window in the reference panel")
      next
    }
    if(length(idx) > 5000){
      log_msg(prot_Apt, " instrument number ", i, ", has more than 5k variants in the 250k window in the reference panel, reducing window")
      distances <- abs(refPanel2$V4[idx] - pos)
      idx <- idx[order(distances)][1:5000]
    }
    rsids_region <- refPanel2$V2[idx]
    chrompos_region <- paste0(refPanel2$V1[idx], ":", refPanel2$V4[idx])
    
    # Exposure region using cis object
    exp_region <- exposure_dat[MarkerName %in% chrompos_region]
    exp_region$rsid <- rsids_region[match(exp_region$MarkerName, chrompos_region)]
    exp_region <- exp_region[!is.na(rsid)]
    if (nrow(exp_region) < 400){
      log_msg(prot_Apt, " instrument number ", i, ", has less than 400 variants in the 250k window")
      next
    }
    log_msg("    Exposure region SNPs: ", nrow(exp_region))
    
    # Format exposure
    exp.dat <- format_data(
      dat = data.frame(exp_region),
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
      log_pval = FALSE,
      min_pval = 1e-320
    )
    exp.dat$exposure <- prot_Apt
    
    # Outcome subset by keyed join
    outcome_dat <- unique(gwas[J(exp.dat$SNP), nomatch = 0], by = "SNP")
    if (nrow(outcome_dat) < 400){
      log_msg(prot_Apt, " instrument number ", i, ", has less than 400 variants in the 250k window in the outcome dataset")
      next
    }
    
    # Harmonise
    dat <- harmonise_data(exp.dat, outcome_dat)
    dat <- dat[!is.na(dat$eaf.exposure),]
    log_msg("    Harmonised SNPs: ", nrow(dat))
    if (!rs %in% dat$SNP){
      log_msg(prot_Apt, " instrument number ", i, ", does not appear after harmonization")
      next
    }
    
    
    # Case/control info
    if (is.null(dat$ncase.outcome) || all(is.na(dat$ncase.outcome))) {
      dat$ncase.outcome <- cases[study]
    }
    if (is.null(dat$ncontrol.outcome) || all(is.na(dat$ncontrol.outcome))) {
      dat$ncontrol.outcome <- dat$samplesize.outcome - dat$ncase.outcome
    }
    
    idx_small <- which(dat$samplesize.outcome < max(dat$samplesize.outcome))
    if (length(idx_small) > 0) {
      dat$ncase.outcome[idx_small] <- round(cases[study] / max(dat$samplesize.outcome) * dat$samplesize.outcome[idx_small])
      dat$ncontrol.outcome[idx_small] <- round(dat$samplesize.outcome[idx_small] - dat$ncase.outcome[idx_small])
    }
    
    # LD matrix
    log_msg("    Computing LD matrix...")
    M <- try(
      my_ld_matrix_local(
        dat$SNP,
        refPanel2,
        plink_bin = plink_bin,
        bfile = plink_ref
      ),
      silent = TRUE
    )
    
    if (inherits(M, "try-error") || is.null(M)) {
      log_msg(prot_Apt, " failed LD matrix for ", rs)
      next
    }
    log_msg("    LD matrix size: ", nrow(M), "x", ncol(M))
    
    # rownames(M) <- gsub("_.*", "", rownames(M))
    # colnames(M) <- gsub("_.*", "", colnames(M))
    
    # Harmonise LD matrix to harmonised data orientation
    hdat <- try(harmonise_ld_dat(dat, M), silent = TRUE)
    if (inherits(hdat, "try-error") || is.null(hdat)) {
      log_msg(prot_Apt, " failed harmonise_ld_dat for ", rs)
      next
    }
    if (is.null(hdat$ld)) {
      log_msg(prot_Apt, " instrument ", i, " no LD matrix after harmonisation for ", rs)
      next
    }
    
    if (nrow(hdat$x) < 350){
      log_msg(prot_Apt, " instrument ", i, " has less than 350 variants after harmonise_ld_dat ", rs)
      next
    }
    
    # Exposure dataset for coloc
    exp_list <- list(
      beta = hdat$x$beta.exposure,
      varbeta = hdat$x$se.exposure^2,
      snp = hdat$x$SNP,
      position = hdat$x$pos.exposure,
      type = "quant",
      N = hdat$x$samplesize.exposure,
      MAF = hdat$x$eaf.exposure,
      LD = hdat$ld
    )
    exp_list$MAF[exp_list$MAF > 0.5] <- 1 - exp_list$MAF[exp_list$MAF > 0.5]
    
    # Outcome dataset for coloc
    out_list <- list(
      beta = hdat$x$beta.outcome,
      varbeta = hdat$x$se.outcome^2,
      snp = hdat$x$SNP,
      position = hdat$x$pos.outcome,
      type = "cc",
      N = hdat$x$samplesize.outcome,
      s = mean(hdat$x$ncase.outcome / hdat$x$samplesize.outcome),
      LD = hdat$ld
    )
    
    # coloc
    coloc_res <- try(
      my.coloc.abf(dataset1 = exp_list, dataset2 = out_list),
      silent = TRUE
    )
    
    if (inherits(coloc_res, "try-error")) {
      log_msg(prot_Apt, " coloc failed for ", rs)
      next
    }
    log_msg("    Coloc done: PP.H4 = ", signif(coloc_res$summary["PP.H4.abf"], 3))
    
    prot_res[[rs]] <- coloc_res
    
    # LD checks
    top.exposure <- rs
    tops.outcome <- hdat$x$SNP[order(hdat$x$pval.outcome, decreasing = FALSE)][1:min(30, nrow(hdat))]
    
    if(top.exposure %in% rownames(hdat$ld)){
      prot_ldold[[rs]] <- any(hdat$ld[top.exposure, tops.outcome] > 0.8, na.rm = TRUE)
      
      top.outcome <- tops.outcome[1]
      prot_ldnew[[rs]] <- (hdat$ld[top.exposure, top.outcome]^2 > 0.8) ||
        any(top.exposure %in% tops.outcome)
    } else {
      prot_ldold[[rs]] <- FALSE
      prot_ldnew[[rs]] <- FALSE
    }
    
    rm(exp.dat, outcome_dat, dat, hdat, M, exp_list, out_list, coloc_res, exp_region)
    gc(verbose = FALSE)
  }
  
  # Save per protein
  if(length(prot_res) > 0){
    saveRDS(prot_res, coloc_outfile)
    saveRDS(prot_ldold, ldold_outfile)
    saveRDS(prot_ldnew, ldnew_outfile)
    
    log_msg("Finished protein ", prot_Apt)
    log_msg("  Instruments processed: ", length(ins$rsid))
    log_msg("  Successful coloc: ", length(prot_res))
    
  } else{
    log_msg(prot_Apt, " ", study, " failed for all instruments")
  }
  
  toc_prot <- Sys.time()
  log_msg("Time for protein: ",
          round(as.numeric(difftime(Sys.time(), tic_prot, units = "mins")), 2),
          " mins")
  
  rm(prot_res, prot_ldold, prot_ldnew, ins, exposure_dat)
  gc(verbose = FALSE)
}

# README + metadata ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic_global, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 6,
  titol = "Colocalization analysis",
  autors = "Sergio",
  lloc = "06_colocalization.R",
  temps = temps,
  descripcio = "For each disease and each protein, colocalization is applied for each independent instrument (IS1) in the cis region (otherwise not applied) using a 250k bp window. Per-protein results are stored under 05_Colocalization/perGene.",
  altres = ""
)

metadata <- list(
  script = "06_colocalization.R",
  timestamp = as.character(Sys.time()),
  study = study,
  chunk = chunk,
  n_proteins_in_chunk = length(chosen.prots),
  R_version = R.version.string
)

write_json(
  metadata,
  path = file.path(
    root_data_dir,
    proj_dir,
    "05_Colocalization",
    paste0("run_metadata_", study, "_chunk", chunk, ".json")
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("06 completed successfully.")
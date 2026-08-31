#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))

# Arguments ----
args <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(args[1])
study_index <- as.integer(args[2])
mode <- args[3]
n_chunks <- as.integer(args[4])
if (is.na(n_chunks) || n_chunks < 1) {
  stop("n_chunks must be a positive integer")
}

# Logging ----

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

# Libraries ----
.libPaths(c("/renv", .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(coloc)
  library(jsonlite)
  library(gwasvcf)
  library(gwasglue)
  library(TwoSampleMR)
})

set.seed(1)

log_msg("=== START COLOCALIZATION ===")
log_msg("Args: ", paste(args, collapse = " | "))
log_msg("Mode: ", mode)

# Configuration ----
source("~/.imiducb/config.R")
source("~/.imiducb/mr_pipeline_config.R")
# Loads: root_data_dir, raw_data_dir, volume_dir, analyses_dir, reference_source_dir, tmp_dir
# Sources: readme_generator.R and register_analysis.R

if (!mode %in% c("pQTL", "eQTL")) {
  stop("Unsupported mode: ", mode)
}

cfg <- load_mr_config(mode)
proj_dir0 <- file.path(analyses_dir, cfg$proj_dir0)
proj_dir <- file.path(analyses_dir, cfg$proj_dir)

source(cfg$tools$tweaked_functions_colocalization)

tic_global <- Sys.time()

# Validate inputs ----

if (length(chunk) != 1 || is.na(chunk) || !chunk %in% seq_len(n_chunks)) {
  stop("chunk must be an integer in 1:", n_chunks)
}

studies <- names(cfg$gwas)

if (length(study_index) != 1 || is.na(study_index) || !study_index %in% seq_along(studies)) {
  stop("study_index must be between 1 and ", length(studies))
}

study <- studies[study_index]

log_msg("Chunk: ", chunk, "/", n_chunks, " | Study index: ", study_index, " (", study, ")")

# Binaries ----

if (root_data_dir == cfg$env$cluster_root) {
  Sys.setenv(PATH = paste(dirname(cfg$tools$plink_bin), Sys.getenv("PATH"), sep = ":"))
  set_plink(path = cfg$tools$plink_bin)
  set_bcftools(path = cfg$tools$bcftools)
  plink_bin <- cfg$tools$plink_bin
} else {
  plink_bin <- cfg$tools$plink_bin
}

# Study parameters ----

cases <- sapply(cfg$gwas, function(x) x$ncase)
ss <- sapply(cfg$gwas, function(x) x$totalsamplesize)

# Inputs ----

refPanel <- fread(
  paste0(cfg$reference$plink_ref, ".bim"),
  header = FALSE
)

plink_ref <- cfg$reference$plink_ref

harm_file <- file.path(proj_dir, "03_Harmonized", paste0(study, ".rds"))

coloc_dir <- file.path(root_data_dir, proj_dir, "05_Colocalization/perGene")
dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)

tmp_file <- file.path(
  root_data_dir,
  proj_dir,
  "05_Colocalization",
  paste0("tmp_", chunk, "_", study_index, ".txt")
)

# Load outcome GWAS ----

gwas_info <- cfg$gwas[[study]]
gwas_path <- file.path(reference_source_dir, gwas_info$file)

if (!file.exists(gwas_path)) {
  stop("Missing GWAS file: ", gwas_path)
}

gwas <- readRDS(gwas_path)
setDT(gwas)

log_msg("Loaded GWAS: ", study, " | N SNPs: ", nrow(gwas))

if ("pos.outcome" %in% names(gwas) && !"pos" %in% names(gwas)) {
  gwas[, pos := as.numeric(pos.outcome)]
} else if ("pos" %in% names(gwas)) {
  gwas[, pos := as.numeric(pos)]
}

gwas <- gwas[!is.na(se.outcome) & !is.na(beta.outcome)]
setkey(gwas, SNP)

# Load harmonized instruments ----

allIns <- readRDS(harm_file)[[1]]

if ("cis" %in% names(allIns)) {
  log_msg("Using column: cis")
  allIns <- allIns[!is.na(allIns$cis) & allIns$cis == TRUE, ]
  
} else if ("Cis" %in% names(allIns)) {
  log_msg("Using column: Cis")
  allIns <- allIns[!is.na(allIns$Cis) & allIns$Cis == TRUE, ]
  
} else if ("Tier" %in% names(allIns)) {
  log_msg("Using column: Tier with cis tiers 1,2,4")
  allIns <- allIns[!is.na(allIns$Tier) & allIns$Tier %in% c(1, 2, 4), ]
  
} else {
  stop(
    "Cannot determine cis instruments: none of 'cis', 'Cis', or 'Tier' found in harmonized object. Columns are: ",
    paste(names(allIns), collapse = ", ")
  )
}

log_msg("allIns rows after cis filter: ", nrow(allIns))

if (nrow(allIns) == 0) {
  log_msg("No cis instruments available")
  quit(save = "no", status = 0)
}

prots <- unique(allIns$exposure)

# Split into chunks ----

chunks <- split(prots, cut(seq_along(prots), n_chunks, labels = FALSE))
chosen_prots <- chunks[[chunk]]

if (length(chosen_prots) == 0) {
  log_msg("Empty chunk")
  quit(save = "no", status = 0)
}

# Main loop ----

for (prot in chosen_prots) {
  
  tic_prot <- Sys.time()
  log_msg("---- Processing ", prot, " ----")
  
  coloc_outfile <- file.path(coloc_dir, paste0(prot, "_", study, "_coloc.rds"))
  ldold_outfile <- file.path(coloc_dir, paste0(prot, "_", study, "_LDcheckOld.rds"))
  ldnew_outfile <- file.path(coloc_dir, paste0(prot, "_", study, "_LDcheck.rds"))
  
  if (file.exists(coloc_outfile) && file.exists(ldold_outfile) && file.exists(ldnew_outfile)) {
    log_msg("Skipping (already done)")
    next
  }
  
  prot_res <- list()
  prot_ldold <- list()
  prot_ldnew <- list()
  
  ins <- allIns[allIns$exposure == prot, ]
  
  if (nrow(ins) == 0 || all(is.na(ins$rsid))) {
    log_msg("No usable rsids")
    next
  }
  
  cisC <- unique(ins$chr.exposure)
  U <- unique(ins$TSS) + 2e6
  L <- unique(ins$TSS) - 2e6
  
  # Exposure data ----
  
  if (mode == "pQTL") {
    
    system(paste0(
      "awk -F ',' '{ if(NR==1 || ($1 == ", cisC, " && $2 < ", U, " && $2 > ", L, ")) { print }}' ",
      root_data_dir, proj_dir,
      "/01_pQTLs_Preprocessed/01_MetaAnalyzed_pQTLs_Clean/",
      prot, ".txt > ", tmp_file
    ))
    
    exposure_dat <- fread(tmp_file)
    
  } else if (mode == "eQTL") {
    
    vcf_path <- file.path(root_data_dir, cfg$vcf_dir, paste0("eqtl-a-", prot, ".vcf.gz"))
    
    if (!file.exists(vcf_path)) {
      log_msg("Missing VCF: ", vcf_path)
      next
    }
    
    vcf <- query_gwas(vcf = vcf_path, chrompos = paste0(cisC, ":", max(1, L), "-", U))
    
    tmp <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, "exposure")
    
    exposure_dat <- data.table(
      MarkerName = paste0(tmp$chr.exposure, ":", tmp$pos.exposure),
      Chromosome = tmp$chr.exposure,
      Position = tmp$pos.exposure,
      Effect = tmp$beta.exposure,
      StdErr = tmp$se.exposure,
      Freq1 = tmp$eaf.exposure,
      Allele1 = tmp$effect_allele.exposure,
      Allele2 = tmp$other_allele.exposure,
      SampleSize = tmp$samplesize.exposure,
      `P-value` = tmp$pval.exposure
    )
    
  } else {
    stop("Unsupported mode: ", mode)
  }
  
  if (nrow(exposure_dat) < 400) {
    log_msg("Too few variants in exposure region")
    next
  }
  
  # Per-instrument loop ----
  
  for (i in seq_len(nrow(ins))) {
    
    rs <- ins$rsid[i]
    chr <- ins$chr.exposure[i]
    pos <- ins$pos.exposure[i]
    
    log_msg("Instrument ", i, " | SNP: ", rs, " | chr: ", chr, " | pos: ", pos)
    
    if (is.na(rs) || is.na(chr) || is.na(pos)) next
    
    refPanel2 <- refPanel[V1 == chr]
    
    idx <- which(abs(refPanel2$V4 - pos) < 25e4)
    
    if (length(idx) < 400) {
      log_msg("Too few variants in reference panel window")
      next
    }
    
    if (length(idx) > 5000) {
      log_msg("More than 5000 variants in reference panel window, reducing to closest 5000")
      idx <- idx[order(abs(refPanel2$V4[idx] - pos))][1:5000]
    }
    
    rsids <- refPanel2$V2[idx]
    chrompos <- paste0(refPanel2$V1[idx], ":", refPanel2$V4[idx])
    
    exp_region <- exposure_dat[MarkerName %in% chrompos]
    exp_region$rsid <- rsids[match(exp_region$MarkerName, chrompos)]
    exp_region <- exp_region[!is.na(rsid)]
    
    if (nrow(exp_region) < 400) {
      log_msg("Too few variants in exposure region after overlap")
      next
    }
    
    exp.dat <- format_data(
      dat = as.data.frame(exp_region),
      snp_col = "rsid",
      beta_col = "Effect",
      se_col = "StdErr",
      eaf_col = "Freq1",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      pval_col = "P-value",
      samplesize_col = "SampleSize",
      chr_col = "Chromosome",
      pos_col = "Position"
    )
    
    exp.dat$exposure <- prot
    
    outcome_dat <- unique(gwas[J(exp.dat$SNP), nomatch = 0], by = "SNP")
    
    if (nrow(outcome_dat) < 400) {
      log_msg("Too few outcome variants in window")
      next
    }
    
    dat <- harmonise_data(exp.dat, outcome_dat)
    dat <- dat[!is.na(dat$eaf.exposure),]
    
    if (!rs %in% dat$SNP) {
      log_msg("Lead SNP absent after harmonisation")
      next
    }
    
    dat$ncase.outcome <- cases[study]
    dat$ncontrol.outcome <- dat$samplesize.outcome - dat$ncase.outcome
    
    idx_small <- which(dat$samplesize.outcome < max(dat$samplesize.outcome))
    if (length(idx_small) > 0) {
      dat$ncase.outcome[idx_small] <- round(
        cases[study] / max(dat$samplesize.outcome) * dat$samplesize.outcome[idx_small]
      )
      dat$ncontrol.outcome[idx_small] <- round(
        dat$samplesize.outcome[idx_small] - dat$ncase.outcome[idx_small]
      )
    }
    
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
      log_msg("LD matrix failed")
      next
    }
    
    hdat <- try(harmonise_ld_dat(dat, M), silent = TRUE)
    if (inherits(hdat, "try-error") || is.null(hdat$ld)) {
      log_msg("harmonise_ld_dat failed")
      next
    }
    
    if (nrow(hdat$x) < 350) {
      log_msg("Too few variants after harmonise_ld_dat")
      next
    }
    
    exp_list <- list(
      beta = hdat$x$beta.exposure,
      varbeta = hdat$x$se.exposure^2,
      snp = hdat$x$SNP,
      position = hdat$x$pos.exposure,
      type = "quant",
      N = hdat$x$samplesize.exposure,
      MAF = pmin(hdat$x$eaf.exposure, 1 - hdat$x$eaf.exposure),
      LD = hdat$ld
    )
    
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
    
    coloc_res <- try(my.coloc.abf(exp_list, out_list), silent = TRUE)
    
    if (inherits(coloc_res, "try-error")) {
      log_msg("coloc failed")
      next
    }
    
    prot_res[[rs]] <- coloc_res
    
    top.exposure <- rs
    tops.outcome <- hdat$x$SNP[order(hdat$x$pval.outcome, decreasing = FALSE)][1:min(30, nrow(hdat$x))]
    
    if (top.exposure %in% rownames(hdat$ld)) {
      prot_ldold[[rs]] <- any(hdat$ld[top.exposure, tops.outcome] > 0.8, na.rm = TRUE)
      
      top.outcome <- tops.outcome[1]
      prot_ldnew[[rs]] <- (hdat$ld[top.exposure, top.outcome]^2 > 0.8) ||
        any(top.exposure %in% tops.outcome)
    } else {
      prot_ldold[[rs]] <- FALSE
      prot_ldnew[[rs]] <- FALSE
    }
  }
  
  # Save per protein ----
  
  if (length(prot_res) > 0) {
    saveRDS(prot_res, coloc_outfile)
    saveRDS(prot_ldold, ldold_outfile)
    saveRDS(prot_ldnew, ldnew_outfile)
    
    log_msg("Finished protein ", prot)
    log_msg("Instruments processed: ", nrow(ins))
    log_msg("Successful coloc: ", length(prot_res))
  } else {
    log_msg(prot, " failed for all instruments")
  }
  
  toc_prot <- Sys.time()
  log_msg("Time for protein: ", round(as.numeric(difftime(toc_prot, tic_prot, units = "mins")), 2), " mins")
  
  rm(prot_res, prot_ldold, prot_ldnew, ins, exposure_dat)
  gc(verbose = FALSE)
}

# Metadata ----

metadata <- list(
  script = "06_colocalization.R",
  mode = mode,
  timestamp = as.character(Sys.time()),
  study = study,
  study_index = study_index,
  chunk = chunk,
  n_chunks = n_chunks,
  n_proteins_in_chunk = length(chosen_prots),
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

# README ----

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic_global, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 6,
  titol = "Colocalization analysis",
  autors = "Sergio",
  lloc = "06_colocalization.R",
  temps = temps,
  descripcio = "For each disease and each cis target, colocalization is applied per independent instrument using a 250kb window around the lead variant. For pQTL mode, exposure data are read from clean meta-analyzed cis-region text files. For eQTL mode, exposure data are extracted from per-gene VCF files. Results are stored in 05_Colocalization/perGene.",
  altres = ""
)

log_msg("=== FINISHED COLOCALIZATION ===")
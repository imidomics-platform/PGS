#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

chunk <- as.integer(args[1])
study_index <- as.integer(args[2])
is_index <- as.integer(args[3])
mode <- ifelse(length(args) >= 4, args[4], "pQTL")
n_chunks <- ifelse(length(args) >= 5, as.integer(args[5]), 1)

.libPaths(c("/renv", .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(jsonlite)
})

# Load config ----
bioinf_root <- Sys.getenv("BIOINF_ROOT", unset = "~/Bioinformatics")
bioinf_root <- path.expand(bioinf_root)

source(file.path(bioinf_root, "Shared/imidomics/R/config.R"))
source(file.path(bioinf_root, "Shared/imidomics/R/readme_generator.R"))
source(file.path(bioinf_root, "Shared/imidomics/R/mr_pipeline_config.R"))

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir <- imidomics_config$raw_data_dir

cfg <- load_mr_config(mode)

proj_dir0 <- cfg$proj_dir0
proj_dir <- cfg$proj_dir
ld_dir <- file.path(root_data_dir, proj_dir, "LD_reference")

source(cfg$tools$tweaked_functions_estimation)

tic <- Sys.time()

# Arguments and study selection ----

gwas_names <- names(cfg$gwas)

if (is.na(study_index) || study_index < 1 || study_index > length(gwas_names)) {
  stop("Invalid study_index. Must be between 1 and ", length(gwas_names))
}

if (is.na(is_index) || !is_index %in% 1:3) {
  stop("Invalid is_index. Must be 1, 2, or 3")
}

if (is.na(n_chunks) || n_chunks < 1) {
  stop("n_chunks must be a positive integer")
}

if (length(chunk) != 1 || is.na(chunk) || !chunk %in% seq_len(n_chunks)) {
  stop("chunk must be an integer in 1:", n_chunks)
}

study <- gwas_names[study_index]
ginfo <- cfg$gwas[[study]]
is <- paste0("IS", is_index)

message("Running mode=", mode, " study=", study, " instrument_set=", is)

# Epidemiology parameters ----

prev <- ginfo$prevalence
ncase <- ginfo$ncase
samplesize <- ginfo$totalsamplesize

# Load MoE model ----

load(file.path(root_data_dir, proj_dir0, "MoE_Model/rf.rdata"))

# Tools and references ----

plink_bin <- cfg$tools$plink_bin
plink_ref <- file.path(root_data_dir, cfg$reference$plink_ref)

refPanel <- fread(paste0(root_data_dir, cfg$reference, ".bim"))

# Load harmonized data ----

harm_all <- readRDS(file.path(root_data_dir, proj_dir, cfg$harmonized, paste0(study, ".rds")))
harm_data <- harm_all[[is_index]]

if (nrow(harm_data) == 0) {
  quit(save = "no")
}

harm_data$exposure <- gsub(".*-", "", harm_data$exposure)
genes <- unique(harm_data$exposure)

# Split into chunks ----

chunks <- split(genes, cut(seq_along(genes), n_chunks, labels = FALSE))
chosen_genes <- chunks[[chunk]]

if (length(chosen_genes) == 0) {
  message("Empty chunk")
  quit(save = "no")
}

message("Chunk: ", chunk, "/", n_chunks, " | Study: ", study, " | IS: ", is)

# Output folder ----

out_dir <- file.path(root_data_dir, proj_dir, "04_MREstimation/perGene")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(out_dir, full.names = TRUE)

# MR estimation loop ----

for (gene in chosen_genes) {
  
  message("Processing gene ", gene)
  
  prefix <- file.path(out_dir, paste(gene, study, is, sep = "_"))
  
  # Check if already computed
  if (file.exists(paste0(prefix, "_Main.rds"))) {
    message("Gene ", gene, " already computed, skipping")
    next()
  }
  
  aux <- harm_data[harm_data$exposure == gene, ]
  aux <- aux[aux$mr_keep == TRUE, ]
  
  if (nrow(aux) == 0) next()
  
  tic2 <- Sys.time()
  tic3 <- Sys.time()
  
  aux$units.outcome <- "log odds"
  aux$prevalence.outcome <- prev
  
  aux$eaf.outcome <- aux$eaf.exposure
  
  aux$ncase.outcome <- ncase
  aux$samplesize.outcome <- samplesize
  aux$ncontrol.outcome <- aux$samplesize.outcome - aux$ncase.outcome
  
  idx <- which(aux$samplesize.outcome < max(aux$samplesize.outcome))
  if (length(idx) > 0) {
    aux$ncase.outcome[idx] <- round(
      ncase / max(aux$samplesize.outcome) * aux$samplesize.outcome[idx]
    )
    aux$ncontrol.outcome[idx] <- round(
      aux$samplesize.outcome[idx] - aux$ncase.outcome[idx]
    )
  }
  
  aux$r.exposure <- get_r_from_bsen(
    b = aux$beta.exposure,
    se = aux$se.exposure,
    n = aux$samplesize.exposure
  )
  
  aux$r.outcome <- get_r_from_lor(
    lor = aux$beta.outcome,
    af = aux$eaf.outcome,
    ncase = aux$ncase.outcome,
    ncontrol = aux$ncontrol.outcome,
    prevalence = aux$prevalence.outcome
  )
  
  suppressWarnings({
    res.dir <- directionality_test(aux)
  })
  saveRDS(res.dir, file = paste0(prefix, "_Steiger.rds"))
  toc3 <- Sys.time()
  print(paste0("Time for Steiger part: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
  
  # Convert to MR input with LD matrix
  tic3 <- Sys.time()
  if ("rsid" %in% names(aux)) aux$SNP <- aux$rsid
  df <- dat_to_MRInput(aux, get_correlation = TRUE)
  if (length(df[[1]]@snps) == 0) next()
  toc3 <- Sys.time()
  print(paste0("Time for LD matrix reading part: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
  
  # Compute all MR methods for >1 instrument with correlation
  tic3 <- Sys.time()
  res <- lapply(df, function(x) {
    if (length(x@betaX) > 2) {
      mr_allmethods(x, method = "main")@Values
    } else {
      NA
    }
  })
  names(res) <- gsub("\\..*", "", names(res))
  toc3 <- Sys.time()
  print(paste0("Time for mr_allmethods part: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
  
  # Compute all MR methods without correlation
  tic3 <- Sys.time()
  
  res2 <- mr(aux)
  res.het <- list()
  res.snps <- list()
  res.egger <- list()
  res.loo <- list()
  
  for (id in names(res)) {
    
    if (unique(res2$nsnp[res2$exposure == id]) > 1) {
      if (is.null(res.het[[id]])) {
        res.het[[id]] <- mr_heterogeneity(aux[aux$exposure == id, ])
      }
      if (is.null(res.egger[[id]])) {
        res.egger[[id]] <- mr_pleiotropy_test(aux[aux$exposure == id, ])
      }
      res.loo[[id]] <- mr_leaveoneout(aux[aux$exposure == id, ])
      if (nrow(aux[aux$exposure == id, ]) == 2) {
        res.loo[[id]] <- rbind(
          mr_singlesnp(aux[aux$exposure == id, ])[1:2, ],
          res.loo[[id]]
        )
      }
    }
    
    res.snps[[id]] <- aux$SNP[aux$exposure == id]
    
    if (!is.data.frame(res[[id]]) && id %in% res2$exposure) {
      res[[id]] <- data.frame(
        res2[res2$exposure == id, c("method", "b", "se")],
        LL = res2[res2$exposure == id, "b"] - 1.96 * res2[res2$exposure == id, "se"],
        UL = res2[res2$exposure == id, "b"] + 1.96 * res2[res2$exposure == id, "se"],
        Pvalue = res2[res2$exposure == id, "pval"]
      )
      names(res[[id]]) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value")
      rownames(res[[id]]) <- NULL
    }
  }
  toc3 <- Sys.time()
  print(paste0("Time for all MR methods without correlation: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
  
  tic3 <- Sys.time()
  for (id in names(res)) {
    if (nrow(aux[aux$exposure == id, ]) > 5) {
      
      # Compute MR raps
      tic3 <- Sys.time()
      res2 <- mr(dat = aux[aux$exposure == id, ], method_list = c("mr_raps"))
      res2$method <- "RAPS"
      toc3 <- Sys.time()
      print(paste0("Time for MR raps: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
      
      aux2 <- data.frame(
        res2[res2$exposure == id, c("method", "b", "se")],
        LL = res2[res2$exposure == id, "b"] - 1.96 * res2[res2$exposure == id, "se"],
        UL = res2[res2$exposure == id, "b"] + 1.96 * res2[res2$exposure == id, "se"],
        Pvalue = res2[res2$exposure == id, "pval"]
      )
      names(aux2) <- names(res[[id]])
      res[[id]] <- rbind(res[[id]], aux2)
      rownames(res[[id]]) <- NULL
      
      # Compute MR-MoE
      tic3 <- Sys.time()
      res2 <- mr_wrapper(dat = aux[aux$exposure == id, ])
      res2 <- mr_moe(res2, rf)
      toc3 <- Sys.time()
      print(paste0("Time for MR-MoE: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
      
      est <- as.data.frame(res2[[1]]$estimates)
      het.aux <- as.data.frame(res2[[1]]$heterogeneity)
      het.aux <- het.aux[
        het.aux$steiger_filtered == est$steiger_filtered[1] &
          het.aux$outlier_filtered == est$outlier_filtered[1] &
          het.aux$method == "IVW",
      ]
      res.het[[id]] <- het.aux
      
      ple.method <- ifelse(grepl("FE", est$method[1]), "FE Egger intercept", "RE Egger intercept")
      ple.aux <- as.data.frame(res2[[1]]$directional_pleiotropy)
      ple.aux <- ple.aux[
        ple.aux$steiger_filtered == est$steiger_filtered[1] &
          ple.aux$outlier_filtered == est$outlier_filtered[1] &
          ple.aux$method == ple.method,
      ]
      res.egger[[id]] <- ple.aux
      
      res2 <- as.data.frame(res2[[1]]$estimates)
      
      aux2 <- data.frame(
        res2[1, c("method", "b", "se")],
        LL = res2[1, "ci_low"],
        UL = res2[1, "ci_upp"],
        Pvalue = res2[1, "pval"]
      )
      
      names(aux2) <- names(res[[id]])
      res[[id]] <- rbind(res[[id]], aux2)
      rownames(res[[id]]) <- NULL
    }
  }
  toc3 <- Sys.time()
  print(paste0("Time for >5 instruments part: ", as.numeric(difftime(toc3, tic3, units = "mins"))))
  
  saveRDS(res, file = paste0(prefix, "_Main.rds"))
  saveRDS(res.het, file = paste0(prefix, "_Heterogeneity.rds"))
  saveRDS(res.snps, file = paste0(prefix, "_SNPs.rds"))
  saveRDS(res.egger, file = paste0(prefix, "_Pleiotropy.rds"))
  saveRDS(res.loo, file = paste0(prefix, "_LOO.rds"))
  
  toc2 <- Sys.time()
  print(paste0("Time for ", gene, ": ", as.numeric(difftime(toc2, tic2, units = "mins"))))
  print(paste0("Total time so far (", which(chosen_genes==gene), "/", length(chosen_genes), "): ", as.numeric(difftime(toc2, tic, units = "mins"))))
}

# Metadata ----

metadata <- list(
  script = "05_estimation.R",
  mode = mode,
  study = study,
  instrument_set = is,
  prevalence = prev,
  ncase = ncase,
  totalsamplesize = samplesize,
  timestamp = as.character(Sys.time()),
  R_version = R.version.string
)

write_json(
  metadata,
  file.path(root_data_dir, proj_dir, "04_MREstimation", paste0("run_metadata_", study, "_", is, ".json")),
  pretty = TRUE,
  auto_unbox = TRUE
)

# README ----

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = paste0(root_data_dir, proj_dir, "/README"),
  subanalisi = 5,
  titol = "MR analysis using the harmonised data",
  autors = "Sergio",
  lloc = "05_estimation.R",
  temps = temps,
  descripcio = "For each disease and each instrument selection method, the MR analysis pipeline is applied. For more than 5 instruments, MendelianRandomization::mr_allmethods, MR-RAPS and MR-MoE are applied. For 5 or fewer instruments, the methods included in TwoSampleMR::mr are applied. Pleiotropy, heterogeneity and leave-one-out sensitivity analyses are stored in per-protein files.",
  altres = ""
)

message("Step 05 finished.")
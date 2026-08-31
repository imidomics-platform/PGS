#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readxl)
})

# Config ----
source("~/Bioinformatics/Shared/imidomics/R/config.R")
source("~/Bioinformatics/Shared/imidomics/R/readme_generator.R")

imidomics_config <- readIMIDomicsConfig()
root_data_dir <- imidomics_config$root_data_dir
raw_data_dir  <- imidomics_config$raw_data_dir

proj_dir <- "/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/"
out_dir  <- paste0(root_data_dir, proj_dir, "/06_Scoring/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tic <- Sys.time()

# Parameters ----
studies <- c(
  "RA_Ishigaki2022",
  "PS_Dand2025",
  "CD_deLange2017",
  "UC_deLange2017",
  "SLE_Bentham2015",
  "PSA_Soomro2022",
  "SS_Sakaue2021",
  "AD_BuduAggrey2023"
)

# External annotation ----
druggable <- read_xlsx(
  paste0(
    root_data_dir,
    "/Data/Molecular/Biological_Annotations/DruggableGenome/NIHMS80906-supplement-Table_S1.xlsx"
  )
)

newlist <- readRDS(
  paste0(
    root_data_dir,
    "/Projects/Internal_Projects/2025_IMX_PaperMR/Summaries/expandedAptamers.rds"
  )
)

ple.regions <- list(
  MHC     = c(6, 29645000, 33365000),
  ABO     = c(9, 136131052, 136150605),
  CFH     = c(1, 196621173, 196716634),
  ARHGEF3 = c(3, 56761448, 57113296),
  BCHE    = c(3, 165490692, 165555211),
  VTN     = c(17, 26694305, 26697328),
  CFD     = c(19, 859664, 863641),
  APOE    = c(19, 45409053, 45412650)
)

# Load gathered results ----
# Keep original nested list format 
cat("Loading gathered results...\n")

STE <- COL3 <- COL2 <- COL <- SNP <- HET <- PLE <- LOO <- RES <- list()

for (name in studies) {
  STE[[name]]  <- list()
  COL3[[name]] <- list()
  COL2[[name]] <- list()
  COL[[name]]  <- list()
  SNP[[name]]  <- list()
  HET[[name]]  <- list()
  PLE[[name]]  <- list()
  LOO[[name]]  <- list()
  RES[[name]]  <- list()
  
  # Coloc outputs from gathered step
  COL[[name]]  <- readRDS(paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_coloc.rds"))
  COL2[[name]] <- readRDS(paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_LDcheckOld.rds"))
  COL3[[name]] <- readRDS(paste0(root_data_dir, proj_dir, "/05_Colocalization/", name, "_LDcheck.rds"))
  
  # MR outputs from gathered step
  for (is in c("IS1", "IS2", "IS3")) {
    RES[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Main.rds"))
    LOO[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_LOO.rds"))
    PLE[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Pleiotropy.rds"))
    HET[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Heterogeneity.rds"))
    SNP[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_SNP.rds"))
    STE[[name]][[is]] <- readRDS(paste0(root_data_dir, proj_dir, "/04_MREstimation/", name, "_", is, "_Steiger.rds"))
  }
}

# Harmonized objects also stay study-dependent
har <- list()
for (name in studies) {
  har[[name]] <- readRDS(paste0(root_data_dir, proj_dir, "/03_Harmonized/", name, ".rds"))
}

# Put MR results together ----
cat("Building res.list...\n")

res.list <- list()

for (name in studies) {
  res <- NULL
  
  for (is in c("IS1", "IS2", "IS3")) {
    aux <- do.call(
      rbind,
      lapply(names(RES[[name]][[is]]), function(x) {
        df <- RES[[name]][[is]][[x]][[x]]
        if (is.null(df)) return(NULL)
        names(df) <- c("Method", "Estimate", "StdErr", "LL", "UL", "PVAL")
        cbind(Prot = x, df)
      })
    )
    
    res <- rbind(res, cbind(IS = is, aux))
  }
  
  res.list[[name]] <- res
}

# Instrument tier + symbols ----
aux <- readRDS(paste0(root_data_dir, proj_dir, "/02_Instruments/ALL.rds"))
aux <- do.call(rbind, aux)

tier <- sapply(unique(aux$exposure), function(x) unique(aux[aux$exposure == x, "Tier"]))
names(tier) <- gsub(".*-", "", names(tier))

featureTable <- readRDS(paste0(root_data_dir, "/Data/Molecular/Proteomics/Plasma/2021_IMX_P4V4_SOMAscan41_IA/Metadata/featureTable.rds"))

apts <- names(tier)
apts <- gsub("_", ".", apts)
apts <- paste0("seq.", apts)

hgnc_names <- featureTable$EntrezGeneSymbol[match(apts, featureTable$AptName)]
names(hgnc_names) <- names(tier)

# Build targets.list ----
# Keep structure and logic close to original
cat("Building targets.list...\n")

targets.list <- list()

for (name in studies) {
  
  res <- res.list[[name]]
  snps <- SNP[[name]]
  
  targets <- data.frame(
    Prot = unique(res$Prot),
    IS = NA,
    NIns1 = 0, NIns2 = 0, NIns3 = 0,
    Directionality = 0,
    PrimaryAverage = 0,
    PrimaryConcordant = FALSE,
    Primary1 = NA, Primary2 = NA, Primary3 = NA,
    SecondaryAverage = 0,
    SecondaryConcordant = NA,
    Secondary1 = NA, Secondary2 = NA, Secondary3 = NA,
    RobustAverage = 0,
    RobustConcordant = NA,
    RAPS1 = NA, RAPS2 = NA, RAPS3 = NA,
    MoE1 = NA, MoE2 = NA, MoE3 = NA,
    AllDirectionConcordant = FALSE,
    ProportionConcordant = 0,
    Heterogeneity_check1 = NA, Heterogeneity_check2 = NA, Heterogeneity_check3 = NA,
    Pleiotropy_check1 = NA, Pleiotropy_check2 = NA, Pleiotropy_check3 = NA,
    LOO_check1 = NA, LOO_check2 = NA, LOO_check3 = NA,
    Steiger_check1 = FALSE, Steiger_check2 = FALSE, Steiger_check3 = FALSE,
    LDcheck_v1 = NA, LDcheck_v2 = NA,
    SimpleColoc_04 = NA, SimpleColoc_08 = NA,
    Tier = 0,
    PleiotropicRegions1 = NA, PleiotropicRegions2 = NA, PleiotropicRegions3 = NA,
    RankingScore = 0,
    row.names = unique(res$Prot)
  )
  
  for (prot in unique(res$Prot)) {
    
    aux <- res[res$Prot == prot, ]
    
    if (sum(is.nan(aux$PVAL)) > 0) aux$PVAL[is.nan(aux$PVAL)] <- 1
    if (any(aux$PVAL == 0)) aux$PVAL[aux$PVAL == 0] <- 1e-100
    
    targets$IS[targets$Prot == prot] <- length(unique(aux$IS))
    
    for (is in unique(aux$IS)) {
      
      IS.num <- as.numeric(gsub("IS", "", is))
      
      # Number of instruments
      if (!is.null(snps[[is]][[prot]][[prot]])) {
        targets[targets$Prot == prot, paste0("NIns", IS.num)] <- length(snps[[is]][[prot]][[prot]])
      }
      
      # Pleiotropic regions
      aux.har <- har[[name]][[IS.num]]
      snps.prot <- NULL
      if (!is.null(SNP[[name]][[is]][[prot]][[prot]])) {
        snps.prot <- SNP[[name]][[is]][[prot]][[prot]]
      }
      
      regs <- character(0)
      
      if (!is.null(snps.prot) && length(snps.prot) > 0) {
        for (snp in seq_along(snps.prot)) {
          snp.chr <- aux.har$chr.exposure[aux.har$rsid == snps.prot[snp]][1]
          snp.pos <- aux.har$pos.exposure[aux.har$rsid == snps.prot[snp]][1]
          
          snp.reg <- sapply(
            ple.regions,
            function(x) as.numeric(snp.chr) == x[1] & snp.pos >= x[2] & snp.pos <= x[3]
          )
          if (any(snp.reg, na.rm = TRUE)) regs <- c(regs, names(snp.reg)[snp.reg])
        }
      }
      
      if (length(regs) > 0) {
        targets[prot, paste0("PleiotropicRegions", IS.num)] <- paste0(unique(regs), collapse = ";")
      }
      
      # Primary / secondary / robust
      if ("Wald ratio" %in% aux$Method[aux$IS == is]) {
        targets[prot, paste0("Primary", IS.num)] <- -log10(aux$PVAL[aux$IS == is & aux$Method == "Wald ratio"])
        
      } else if ("Inverse variance weighted" %in% aux$Method[aux$IS == is]) {
        a <- -log10(aux$PVAL[aux$IS == is & aux$Method == "Inverse variance weighted"])
        targets[prot, paste0("Primary", IS.num)] <- a
        
        het <- HET[[name]][[is]][[prot]][[prot]]
        if (!is.null(het)) {
          if ("pval" %in% names(het)) {
            targets[prot, paste0("Heterogeneity_check", IS.num)] <- ifelse(het$pval[1] > 0.05, TRUE, FALSE)
          } else if ("method" %in% names(het)) {
            targets[prot, paste0("Heterogeneity_check", IS.num)] <- ifelse(
              het[het$method == "Inverse variance weighted", "Q_pval"] > 0.05,
              TRUE, FALSE
            )
          }
        }
        
        ple <- PLE[[name]][[is]][[prot]][[prot]]
        if (!is.null(ple) && nrow(ple) > 0 && !is.na(ple$pval[1])) {
          targets[prot, paste0("Pleiotropy_check", IS.num)] <- ifelse(ple$pval[1] > 0.05, TRUE, FALSE)
        }
        
        loo <- LOO[[name]][[is]][[prot]][[prot]]
        if (!is.null(loo) && all(!is.nan(loo$p))) {
          targets[prot, paste0("LOO_check", IS.num)] <- ifelse(
            all(loo$p < 0.1 & (all(sign(loo$b) == 1) | all(sign(loo$b) == -1))),
            TRUE, FALSE
          )
        }
        
      } else if ("IVW" %in% aux$Method[aux$IS == is]) {
        targets[prot, paste0("Primary", IS.num)] <- -log10(aux$PVAL[aux$IS == is & aux$Method == "IVW"])
        targets[prot, paste0("Secondary", IS.num)] <- mean(
          -log10(aux$PVAL[aux$IS == is & aux$Method %in% c("MR-Egger", "Weighted median", "Simple median")])
        )
        
        het <- HET[[name]][[is]][[prot]][[prot]]
        if (!is.null(het)) {
          if ("pval" %in% names(het)) {
            targets[prot, paste0("Heterogeneity_check", IS.num)] <- ifelse(het$pval[1] > 0.05, TRUE, FALSE)
          } else if ("method" %in% names(het)) {
            targets[prot, paste0("Heterogeneity_check", IS.num)] <- ifelse(
              het[het$method == "Inverse variance weighted", "Q_pval"] > 0.05,
              TRUE, FALSE
            )
          }
        }
        
        ple <- PLE[[name]][[is]][[prot]][[prot]]
        if (!is.null(ple) && nrow(ple) > 0 && !is.na(ple$pval[1])) {
          targets[prot, paste0("Pleiotropy_check", IS.num)] <- ifelse(ple$pval[1] > 0.05, TRUE, FALSE)
        }
        
        loo <- LOO[[name]][[is]][[prot]][[prot]]
        if (!is.null(loo) && all(!is.nan(loo$p))) {
          targets[prot, paste0("LOO_check", IS.num)] <- ifelse(
            all(loo$p < 0.1 & (all(sign(loo$b) == 1) | all(sign(loo$b) == -1))),
            TRUE, FALSE
          )
        }
      }
      
      if ("RAPS" %in% aux$Method[aux$IS == is]) {
        targets[prot, paste0("RAPS", IS.num)] <- -log10(aux$PVAL[aux$IS == is & aux$Method == "RAPS"])
        # keep original positional assumption for MoE row
        targets[prot, paste0("MoE", IS.num)] <- -log10(aux$PVAL[aux$IS == is][7])
        targets[prot, paste0("Secondary", IS.num)] <- mean(
          -log10(aux$PVAL[aux$IS == is & aux$Method %in% c("MR-Egger", "Weighted median", "Simple median")])
        )
      }
      
      # Steiger
      ste <- STE[[name]][[is]][[prot]]
      if (!is.null(ste)) {
        ste <- ste[ste$exposure == prot, ]
        if (nrow(ste) > 0) {
          if (sum(ste$correct_causal_direction == FALSE, na.rm = TRUE) / length(ste$correct_causal_direction) < 0.2) {
            targets[prot, paste0("Steiger_check", IS.num)] <- TRUE
          }
        }
      }
    }
    
    # Primary concordance
    estimates <- aux$Estimate[aux$Method %in% c("Wald ratio", "Inverse variance weighted", "IVW")]
    pvals <- aux$PVAL[aux$Method %in% c("Wald ratio", "Inverse variance weighted", "IVW")]
    
    if (length(estimates) > 0 && (all(sign(estimates) == 1) | all(sign(estimates) == -1))) {
      targets[prot, "PrimaryConcordant"] <- TRUE
      targets[prot, "PrimaryAverage"] <- mean(-log10(pvals))
    } else {
      targets[prot, "Primary1"] <- 0
      targets[prot, "Primary2"] <- 0
      targets[prot, "Primary3"] <- 0
    }
    
    # Robust concordance
    if ("RAPS" %in% aux$Method) {
      idx <- which(aux$Method == "RAPS")
      if (length(idx) > 0 && all(idx + 1 <= nrow(aux))) {
        estimates <- aux$Estimate[c(idx, idx + 1)]
        pvals <- aux$PVAL[c(idx, idx + 1)]
        if (all(sign(estimates) == 1) | all(sign(estimates) == -1)) {
          targets[prot, "RobustConcordant"] <- TRUE
          targets[prot, "RobustAverage"] <- mean(-log10(pvals))
        } else {
          targets[prot, "RobustConcordant"] <- FALSE
          targets[prot, c("RAPS1", "RAPS2", "RAPS3", "MoE1", "MoE2", "MoE3")] <- 0
        }
      }
    }
    
    # Secondary concordance
    if ("IVW" %in% aux$Method) {
      idx <- which(aux$Method %in% c("Simple median", "Weighted median", "MR-Egger"))
      if (length(idx) > 0) {
        estimates <- aux$Estimate[idx]
        pvals <- aux$PVAL[idx]
        if (all(sign(estimates) == 1) | all(sign(estimates) == -1)) {
          targets[prot, "SecondaryConcordant"] <- TRUE
          targets[prot, "SecondaryAverage"] <- mean(-log10(pvals))
        } else {
          targets[prot, "SecondaryConcordant"] <- FALSE
          targets[prot, c("Secondary1", "Secondary2", "Secondary3")] <- 0
        }
      }
    }
    
    # All direction concordant
    if ("(intercept)" %in% aux$Method) {
      estimates <- aux$Estimate[aux$Method != "(intercept)"]
    } else {
      estimates <- aux$Estimate
    }
    if (length(estimates) > 0 && (all(sign(estimates) == 1) | all(sign(estimates) == -1))) {
      targets[prot, "AllDirectionConcordant"] <- TRUE
    }
    
    # Proportion concordant
    if (isTRUE(targets[prot, "PrimaryConcordant"])) {
      correct.sign <- unique(sign(aux$Estimate[aux$Method %in% c("Wald ratio", "Inverse variance weighted", "IVW")]))
      correct.sign <- correct.sign[!is.na(correct.sign)][1]
      targets[prot, "Directionality"] <- correct.sign
      
      aux2 <- aux
      if ("(intercept)" %in% aux2$Method) {
        aux2 <- aux2[aux2$Method != "(intercept)", ]
      }
      
      correct.estimates <- aux2$Estimate[aux2$PVAL < 0.05 & sign(aux2$Estimate) == correct.sign]
      targets[prot, "ProportionConcordant"] <- length(correct.estimates) / nrow(aux2)
    }
    
    # Coloc
    if (!is.null(COL[[name]][[prot]])) {
      probs <- sapply(COL[[name]][[prot]], function(x) x$summary["PP.H4.abf"])
      targets[prot, "SimpleColoc_04"] <- any(probs > 0.4)
      targets[prot, "SimpleColoc_08"] <- any(probs > 0.8)
    }
    
    if (!is.null(COL2[[name]][[prot]])) {
      targets[prot, "LDcheck_v1"] <- any(COL2[[name]][[prot]] == TRUE)
    }
    
    if (!is.null(COL3[[name]][[prot]])) {
      targets[prot, "LDcheck_v2"] <- any(COL3[[name]][[prot]] == TRUE)
    }
  }
  
  targets.list[[gsub("_.*", "", name)]] <- targets
}

saveRDS(targets.list, file = paste0(root_data_dir, proj_dir, "/06_Scoring/statisticsObject.rds"))

# Scoring ----
cat("Scoring targets...\n")

for (name in studies) {
  
  targets <- targets.list[[gsub("_.*", "", name)]]
  
  t <- -log10(0.05 / nrow(targets))
  targets$RankingScore <- 0
  
  # Consistent directionality
  con <- targets$Prot[which(targets$AllDirectionConcordant == TRUE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2
  
  # Proportion concordant
  targets$RankingScore <- targets$RankingScore + 2 * targets$ProportionConcordant
  
  # Robust significance
  for (j in c(1.3, t, 7, 10)) {
    con <- targets$Prot[which(targets$RobustConcordant == TRUE & targets$RobustAverage > j)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2
  }
  if (length(con) > 0) {
    normps <- targets[con, "RobustAverage"] / max(targets[con, "RobustAverage"])
    targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 6 * normps
  }
  
  # Primary only
  for (j in c(1.3, t, 7, 10)) {
    con <- targets$Prot[which(
      is.na(targets$RobustConcordant) &
        is.na(targets$SecondaryConcordant) &
        targets$PrimaryConcordant == TRUE &
        targets$PrimaryAverage > j
    )]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  if (length(con) > 0) {
    normps <- targets[con, "PrimaryAverage"] / max(targets[con, "PrimaryAverage"])
    targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 6 * normps
  }
  
  # Primary + Secondary
  for (j in c(1.3, t, 7, 10)) {
    con <- targets$Prot[which(
      is.na(targets$RobustConcordant) &
        !is.na(targets$SecondaryConcordant) &
        targets$PrimaryConcordant == TRUE &
        targets$PrimaryAverage > j
    )]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  if (length(con) > 0) {
    normps <- targets[con, "PrimaryAverage"] / max(targets[con, "PrimaryAverage"])
    targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2 * normps
  }
  
  for (j in c(1.3, t, 7, 10)) {
    con <- targets$Prot[which(
      is.na(targets$RobustConcordant) &
        targets$SecondaryConcordant == TRUE &
        targets$SecondaryAverage > j
    )]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  if (length(con) > 0) {
    normps <- targets[con, "SecondaryAverage"] / max(targets[con, "SecondaryAverage"])
    targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2 * normps
  }
  
  # Checks
  for (i in 1:3) {
    con <- targets$Prot[which(targets[, paste0("Heterogeneity_check", i)] == FALSE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 1
    
    con <- targets$Prot[which(targets[, paste0("Heterogeneity_check", i)] == TRUE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  
  for (i in 1:3) {
    con <- targets$Prot[which(targets[, paste0("Pleiotropy_check", i)] == FALSE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 1
    
    con <- targets$Prot[which(targets[, paste0("Pleiotropy_check", i)] == TRUE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  
  for (i in 1:3) {
    con <- targets$Prot[which(targets[, paste0("LOO_check", i)] == FALSE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 1
    
    con <- targets$Prot[which(targets[, paste0("LOO_check", i)] == TRUE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  }
  
  for (i in 1:3) {
    con <- targets$Prot[which(targets[, paste0("Steiger_check", i)] == FALSE)]
    if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 3
  }
  
  # LD check
  con <- targets$Prot[which(targets$LDcheck_v1 == FALSE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 3
  
  con <- targets$Prot[which(targets$LDcheck_v1 == TRUE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 3
  
  # Coloc
  con <- targets$Prot[which(targets$SimpleColoc_04 == FALSE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 2
  con <- targets$Prot[which(targets$SimpleColoc_04 == TRUE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2
  
  con <- targets$Prot[which(targets$SimpleColoc_08 == FALSE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] - 2
  con <- targets$Prot[which(targets$SimpleColoc_08 == TRUE)]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2
  
  # Tier
  targets$Tier <- tier[match(targets$Prot, names(tier))]
  targets$Cis <- targets$Tier %in% c(1, 2, 4)
  
  con <- targets$Prot[targets$Tier == 1]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 8
  con <- targets$Prot[targets$Tier == 2]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 6
  con <- targets$Prot[targets$Tier == 3]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 4
  con <- targets$Prot[targets$Tier == 4]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 2
  con <- targets$Prot[targets$Tier == 5]
  if (length(con) > 0) targets[con, "RankingScore"] <- targets[con, "RankingScore"] + 1
  
  # Pleiotropic region penalty
  targets$RankingScore <- targets$RankingScore -
    2 * !is.na(targets$PleiotropicRegions1) -
    2 * !is.na(targets$PleiotropicRegions2) -
    2 * !is.na(targets$PleiotropicRegions3)
  
  # Sort
  targets <- targets[order(targets$RankingScore, decreasing = TRUE), ]
  
  # Zero if no signal at all
  con <- targets$Prot[
    !apply(
      targets[, c(
        paste0("Primary", 1:3),
        paste0("Secondary", 1:3),
        paste0("RAPS", 1:3),
        paste0("MoE", 1:3)
      )],
      1,
      function(x) any(x > 1.3, na.rm = TRUE)
    )
  ]
  if (length(con) > 0) targets[con, "RankingScore"] <- 0
  
  con <- targets$Prot[targets$Directionality == 0]
  if (length(con) > 0) targets[con, "RankingScore"] <- 0
  
  targets$RankingScore[targets$RankingScore < 0] <- 0
  
  # Extra annotations
  targets <- data.frame(
    Symbol = hgnc_names[match(targets$Prot, names(hgnc_names))],
    targets
  )
  
  targets$Druggable <- targets$Symbol %in% druggable$hgnc_names
  targets$NewAptamer <- targets$Prot %in% gsub("\\.", "_", gsub("seq.", "", newlist$newAptamers))
  targets$NewSymbol <- targets$Symbol %in% newlist$newSymbols
  
  saveRDS(targets, file = paste0(out_dir, "/", name, ".rds"))
  
  # preserve updated object in memory too
  targets.list[[gsub("_.*", "", name)]] <- targets
}

# Save final combined object
saveRDS(targets.list, file = paste0(out_dir, "/targets_list_scored.rds"))

# README ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = paste0(root_data_dir, proj_dir, "/README"),
  subanalisi = 8,
  titol = "Final scoring",
  autors = "Sergio",
  lloc = "08_scoring.R",
  temps = temps,
  descripcio = "Final prioritization score integrating MR, colocalization, LD checks, sensitivity analyses and instrument tier.",
  altres = ""
)

cat("08_scoring.R finished successfully.\n")
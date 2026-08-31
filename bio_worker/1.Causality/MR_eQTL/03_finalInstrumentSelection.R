#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(readxl)
  library(jsonlite)
})

# Configuration ----
source("~/.imiducb/config.R")
source("~/.imiducb/mr_pipeline_config.R")
# Loads: root_data_dir, raw_data_dir, volume_dir, analyses_dir, reference_source_dir, tmp_dir
# Sources: readme_generator.R and register_analysis.R

mode <- "eQTL"
cfg <- load_mr_config(mode)
proj_dir <- file.path(analyses_dir, cfg$proj_dir)
proj_dir0 <- file.path(analyses_dir, cfg$proj_dir0)

message("Running instrument filtering for mode: ", mode)

tic <- Sys.time()

# TSS annotation ----

ensembl.tss <- readRDS(paste0(root_data_dir, proj_dir, "/Annotations/All_ENS_SYMB_TSS_v37.rds"))

# Pleiotropic regions ----

transhub <- readxl::read_xlsx(path = file.path(root_data_dir, proj_dir, "Annotations/Vosa_SupTables.xlsx"), sheet = 7, skip = 1)
ple.snps <- transhub$SNP[transhub$`nr. of genes affected` > 5]

# Identify pleiotropic SNPs across genes ----

ins_full <- readRDS(file.path(root_data_dir, proj_dir, "02_Instruments/IS3_full.rds"))
tab <- table(ins_full$SNP)
pleiotropics <- names(which(tab > 5))

# Instrument filtering ----

trackGenes <- vector(mode = "list", length = 3)
INS <- list()

for (i in 1:3) {
  message("Processing IS", i)

  ins <- readRDS(file.path(root_data_dir, proj_dir, paste0("02_Instruments/IS", i, "_full.rds")))
  ins$exposure <- gsub(".*a-", "", ins$exposure)

  # Remove pleiotropic variants (>5 genes)
  ins <- ins[!ins$SNP %in% pleiotropics, ]

  # Keep strongest SNP when duplicated (extreme p-values)
  ins <- ins[order(ins$pval.exposure), ]
  ins <- ins[!(ins$pval.exposure > 1e-200 & duplicated(ins$SNP)), ]
  ins <- ins[order(ins$exposure), ]

  # Annotate pleiotropic SNPs (trans hotspots)
  ins$PleSNP <- ins$SNP %in% ple.snps

  # Create TSS variable
  ins$TSS <- NA

  # Remove pleiotropic-region SNPs per gene if alternatives exist
  for (g in unique(ins$exposure)) {
    idx <- which(ins$exposure == g)

    ins[idx, "TSS"] <- mean(ensembl.tss$transcription_start_site[ensembl.tss$ensembl_gene_id == g])

    aux <- ins[idx, ]

    if (sum(aux$PleSNP) > 0 && sum(aux$PleSNP) < nrow(aux)) {
      trackGenes[[i]][[g]] <- aux$SNP[aux$PleSNP]

      ins <- ins[-idx[aux$PleSNP]]
    }
  }

  ins$TSS <- round(ins$TSS)
  ins$rsid <- ins$SNP
  INS[[i]] <- ins
}

# Save instruments ----

ins_dir <- file.path(root_data_dir, proj_dir, "02_Instruments")
dir.create(ins_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(INS, file = file.path(ins_dir, "ALL.rds"))
saveRDS(trackGenes, file = file.path(ins_dir, "genesToReanalyze.rds"))

# Build LD sets ----

ld_sets <- list()
for (i in 1:3) {
  ins <- INS[[i]]
  for (g in unique(ins$exposure)) {
    snps <- sort(ins$rsid[ins$exposure == g])
    if (length(snps) > 1) {
      key <- paste0(g, "_IS", i)
      ld_sets[[key]] <- snps
    }
  }
}

message("Total LD sets: ", length(ld_sets))

# LD precomputation ----

ld_dir <- file.path(root_data_dir, proj_dir, "LD_reference")
dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)

plink_bin <- cfg$tools$plink_bin
plink_ref <- file.path(root_data_dir, cfg$reference$plink_ref)
refPanel <- fread(paste0(plink_ref, ".bim"))

tmp_snps <- file.path(ld_dir, "tmp_snps.txt")
tmp_out <- file.path(ld_dir, "tmp_ld")

for (name in names(ld_sets)) {
  outfile <- file.path(ld_dir, paste0(name, ".rds"))
  # if (file.exists(outfile)) next

  snps <- ld_sets[[name]]
  snps <- snps[snps %in% refPanel$V2]

  message("Computing LD for ", name, " (", length(snps), " SNPs)")

  write.table(snps, file = tmp_snps, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cmd <- paste(plink_bin, "--bfile", plink_ref, "--extract", tmp_snps, "--r square", "--out", tmp_out)
  system(cmd)
  ld_file <- paste0(tmp_out, ".ld")
  if (!file.exists(ld_file)) {
    warning("LD file missing for ", name)
    next
  }
  ld <- as.matrix(fread(ld_file))
  saveRDS(list(snps = snps, ld = ld), outfile)
}

# Cleanup temporary files ----

file.remove(tmp_snps)
file.remove(paste0(tmp_out, ".ld"))
file.remove(paste0(tmp_out, ".log"))
file.remove(paste0(tmp_out, ".nosex"))

# Register analysis ----
register_analysis(
  analysis_type = "instrument_selection",
  analysis_subtype = "eQTL",
  analysis_method = "custom",
  analysis_design = "final_selection",
  parameters = list(mode = mode),
  datasets = list(),
  reference_sources = as.list(unique(c(cfg$vcf_id, cfg$reference$plink_ref_id))),
  input_files = list(),
  analysis_path = ins_dir,
  analysis_code = script_path,
  status = "done",
  computational_running_time = temps,
  description = "eQTL instruments are filtered to remove pleiotropic variants and trans-hotspot regions, then LD matrices are precomputed for each gene and instrument set."
)

# README ----

toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 3,
  titol = "Instrument filtering and LD precomputation",
  autors = "Sergio",
  lloc = "03_finalInstrumentSelection.R",
  temps = temps,
  descripcio = "eQTL instruments are filtered to remove pleiotropic variants (>5 genes) and trans-hotspot regions. LD matrices are precomputed per gene and instrument set using 1000G EUR.",
  altres = ""
)

message("Step 03 completed successfully.")

#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Configuration ----
source("~/.imiducb/config.R")
source("~/.imiducb/mr_pipeline_config.R")
# Loads: root_data_dir, raw_data_dir, volume_dir, analyses_dir, reference_source_dir, tmp_dir
# Sources: readme_generator.R and register_analysis.R

set.seed(1)

# Project settings ----
cfg <- load_mr_config("pQTL")
proj_dir <- file.path(analyses_dir, cfg$proj_dir)
proj_dir0 <- file.path(analyses_dir, cfg$proj_dir0)

tic <- Sys.time()

# Set Plink ----
if (root_data_dir == "/mnt/proteome/bioinformatics") {
  # Binary
  plink_bin <- "/mnt/nfs_home/smartinez/tools/plink/plink"
  # set_plink(path = "/mnt/nfs_home/smartinez/tools/plink/plink")
  Sys.setenv(PATH = paste("/home/sergio/bin/", Sys.getenv("PATH"), sep = ":"))

  # Reference
  plink_ref <- paste0(root_data_dir, "/Data/Molecular/Biological_Annotations/1000G/EUR")
  # sudo mkdir /tmp/ramdisk
  # sudo mount -t tmpfs -o size=1G tmpfs /tmp/ramdisk
  # plink_ref <- "/dev/shm/tmp/EUR"
} else {
  # Pending
  # Sys.setenv(PATH=paste("/home/sergio/bcftools/", Sys.getenv("PATH"), sep=":"))
}

# Pleiotropic regions ----
ple.regions <- list(
  MHC = c(6, 29645000, 33365000),
  ABO = c(9, 136131052, 136150605),
  CFH = c(1, 196621173, 196716634),
  ARHGEF3 = c(3, 56761448, 57113296),
  BCHE = c(3, 165490692, 165555211),
  VTN = c(17, 26694305, 26697328),
  CFD = c(19, 859664, 863641),
  APOE = c(19, 45409053, 45412650)
)

# Merge chunk outputs ----
INS <- vector("list", 3)

for (i in 1:8) {
  fn <- file.path(root_data_dir, proj_dir, "02_Instruments", paste0("ALL_full", i, ".rds"))

  if (!file.exists(fn)) stop("Missing chunk file: ", fn)

  a <- readRDS(fn)

  for (j in 1:3) {
    if (is.null(INS[[j]])) {
      INS[[j]] <- a[[j]]
    } else {
      INS[[j]] <- rbind(INS[[j]], a[[j]])
    }
  }
}

saveRDS(INS, file = file.path(root_data_dir, proj_dir, "02_Instruments", "ALL_full.rds"))

# Filter pleiotropic instruments ----
INS <- readRDS(file = file.path(root_data_dir, proj_dir, "02_Instruments", "ALL_full.rds"))

# SNPs associated with >5 proteins
tab <- table(INS[[3]]$SNP)
pleiotropics <- names(tab[tab > 5])

trackProteins <- vector("list", 3)

for (i in 1:3) {
  ins <- INS[[i]]

  # Remove highly pleiotropic variants
  ins <- ins[!ins$SNP %in% pleiotropics, ]

  # Remove duplicates if strongest pQTL identifiable
  ins <- ins[order(ins$pval.exposure), ]
  ins <- ins[!(ins$pval.exposure > 1e-200 & duplicated(ins$SNP)), ]
  ins <- ins[order(ins$exposure), ]

  # Annotate pleiotropic regions ----
  ins$PleRegion <- NA
  for (j in seq_along(ple.regions)) {
    chr <- ple.regions[[j]][1]
    start <- ple.regions[[j]][2]
    end <- ple.regions[[j]][3]

    idx <- ins$chr.exposure == chr &
      ins$pos.exposure >= start &
      ins$pos.exposure <= end

    ins$PleRegion[idx] <- names(ple.regions)[j]
  }

  # Remove pleiotropic-region SNPs if other instruments exist ----
  for (p in unique(ins$exposure)) {
    aux <- ins[ins$exposure == p, ]
    nPle <- sum(!is.na(aux$PleRegion))
    if (nPle > 0 & nPle < nrow(aux)) {
      trackProteins[[i]][[p]] <- aux$SNP[!is.na(aux$PleRegion)]
      ins <- ins[!(ins$exposure == p & !is.na(ins$PleRegion)), ]
    }
  }
  INS[[i]] <- ins
}

# Save final instruments ----
saveRDS(INS, file = file.path(root_data_dir, proj_dir, "02_Instruments", "ALL.rds"))
saveRDS(trackProteins, file = file.path(root_data_dir, proj_dir, "02_Instruments", "proteinsToReanalyze.rds"))

# Compute LD matrix between all instruments ----
ld_sets <- list()
for (i in 1:3) {
  ins <- INS[[i]]
  for (prot in unique(ins$exposure)) {
    snps <- sort(ins$rsid[ins$exposure == prot])
    if (length(snps) > 1) {
      key <- paste0(prot, "_IS", i)
      ld_sets[[key]] <- snps
    }
  }
}
length(ld_sets)

ld_dir <- file.path(root_data_dir, proj_dir, "LD_reference")
dir.create(ld_dir, showWarnings = FALSE)
for (name in names(ld_sets)) {
  snps <- ld_sets[[name]]
  outfile <- file.path(ld_dir, paste0(name, ".rds"))
  if (file.exists(outfile)) next
  write.table(snps, file = file.path(ld_dir, "tmp_snps.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  system(paste(
    plink_bin,
    "--bfile", plink_ref,
    "--extract", file.path(ld_dir, "tmp_snps.txt"),
    "--r square",
    "--out", file.path(ld_dir, "tmp_ld")
  ))
  ld <- as.matrix(fread(file.path(ld_dir, "tmp_ld.ld")))
  saveRDS(list(snps = snps, ld = ld), outfile)
}

# README ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 3,
  titol = "Instrument selection II",
  autors = "Sergio",
  lloc = "03_finalInstrumentSelection.R",
  temps = temps,
  descripcio = "The selected instruments are merged across chunks. Variants associated with >5 proteins are removed. For variants appearing in 2–5 proteins, only the strongest pQTL is kept when distinguishable. Instruments located in known pleiotropic regions are removed when alternative instruments exist. For each final set of instruments we compute the LD correlation matrix with plink and keep it in separate files for fast access later.",
  altres = ""
)

# Register analysis ----
register_analysis(
  analysis_type = "instrument_selection",
  analysis_subtype = "pQTL",
  analysis_method = "custom",
  analysis_design = "final_selection",
  parameters = list(),
  datasets = list(),
  reference_sources = as.list(unique(c(unname(unlist(lapply(cfg$gwas, function(x) x$reference_id))), cfg$reference$plink_ref_id))),
  input_files = list(),
  analysis_path = file.path(proj_dir, "02_Instruments"),
  analysis_code = script_path,
  status = "done",
  computational_running_time = temps,
  description = "Merge chunk outputs, remove pleiotropic variants and pleiotropic regions when alternatives exist, and compute LD matrices for the final pQTL instrument sets."
)

message("Step 03 completed successfully.")

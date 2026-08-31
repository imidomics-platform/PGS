#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))
args <- commandArgs(trailingOnly = TRUE)
chunk <- suppressWarnings(as.integer(args[1]))
if (length(chunk) != 1 || is.na(chunk) || !chunk %in% 1:8) {
  stop("chunk argument must be an integer in 1:8, e.g. Rscript 02_instrumentSelection.R 1")
}

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(stringr)
  library(rtracklayer)
  library(liftOver)
  library(ieugwasr)
  library(genetics.binaRies)
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

# Keep your current PATH logic (institution-specific)
if (root_data_dir == "/mnt/proteome/bioinformatics/") {
  Sys.setenv(PATH = paste("~/tools/", Sys.getenv("PATH"), sep = ":"))
  # Tooling / environment
  plink_bin <- "~/tools/plink/plink"
} else {
  Sys.setenv(PATH = paste("/home/sergio/METAL/build/bin/", Sys.getenv("PATH"), sep = ":"))
  Sys.setenv(PATH = paste("/home/sergio/bcftools/", Sys.getenv("PATH"), sep = ":"))
}

# Paths / IO ----
in_preproc_dir <- file.path(root_data_dir, proj_dir, "01_pQTLs_Preprocessed", "02_Filtered_Formatted")
files_rds_path <- file.path(root_data_dir, proj_dir, "01_pQTLs_Preprocessed", "files_for_clumping.rds")

out_dir <- file.path(root_data_dir, proj_dir, "02_Instruments")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "log.txt")

# temp files (chunk-scoped; safe for SLURM array)
fn_positions <- file.path(out_dir, paste0("input_temp_", chunk))
fn_bcf_out <- file.path(out_dir, paste0("output_temp_", chunk))

if (!file.exists(files_rds_path)) {
  stop("Missing input file list: ", files_rds_path)
}

# LD reference panel ----
ramdisk_prefix <- "/tmp/ramdisk/EUR"
shared_prefix <- file.path(root_data_dir, "Data/Molecular/Biological_Annotations/1000G/EUR")

if (file.exists(paste0(ramdisk_prefix, ".bim"))) {
  message("Using LD reference from /tmp/ramdisk")
  plink_ref_prefix <- ramdisk_prefix
  eur_bim_path <- paste0(ramdisk_prefix, ".bim")
} else {
  message("Using LD reference from shared storage")
  plink_ref_prefix <- shared_prefix
  eur_bim_path <- paste0(shared_prefix, ".bim")
}

# VEP runtime detection ----
vep_mode <- "none"
if (nzchar(Sys.which("apptainer"))) {
  vep_mode <- "apptainer"
  message("Using VEP via Apptainer")
} else if (nzchar(Sys.which("docker"))) {
  vep_mode <- "docker"
  message("Using VEP via Docker")
} else if (nzchar(Sys.which("vep"))) {
  vep_mode <- "local"
  message("Using local VEP installation")
} else {
  warning("No VEP runtime detected. PAV annotation will be skipped.")
}

# Load objects ----
featureTable <- readRDS(file.path(root_data_dir, "/Data/Molecular/Proteomics/Plasma/2021_IMX_P4V4_SOMAscan41_IA/Metadata/featureTable.rds"))

# 1000G reference panel bim
if (!file.exists(eur_bim_path)) {
  stop("Missing EUR.bim at ", eur_bim_path, " (expected /tmp/ramdisk/EUR.bim).")
}
eur.bim <- fread(eur_bim_path, header = FALSE)
eur.bim$ID <- paste0(eur.bim$V1, ":", eur.bim$V4)
rsyes <- grepl("^rs", eur.bim$V2)
panel.ids <- eur.bim$ID[rsyes]

# TSS mapping
ensembl.tss <- readRDS(file.path(root_data_dir, "/Data/Molecular/Biological_Annotations/refTSS/SOMAscanProts_FinalTSSMapping_lift37.rds"))
ensembl.tss <- ensembl.tss[ensembl.tss$CHR %in% c(1:22, "X", "Y"), ]

# ENSID dict for VEP filtering
dict <- readRDS(file.path(root_data_dir, "/Data/Molecular/Biological_Annotations/GeneIDs_Dictionaries/SOMAscanProtsV41_ensemblIDs.rds"))

# dbSNP assembly report (load once)
chroms_path <- file.path(root_data_dir, "/Data/Molecular/Biological_Annotations/dbsnp/GCF_000001405.25.assembly_report.txt")
dbsnp_vcf <- file.path(root_data_dir, "/Data/Molecular/Biological_Annotations/dbsnp/GCF_000001405.25.gz")

if (!file.exists(chroms_path) || !file.exists(dbsnp_vcf)) {
  stop("Missing dbSNP resources: ", chroms_path, " or ", dbsnp_vcf)
}
chroms <- fread(chroms_path)

# Chunk selection ----
files <- readRDS(files_rds_path)

s <- floor(seq(from = 1, to = length(files), length.out = 8 + 1))
if (chunk == 1) {
  chosen.files <- files[1:s[2]]
} else {
  chosen.files <- files[(s[chunk] + 1):s[chunk + 1]]
}

message("Chunk ", chunk, " processing ", length(chosen.files), " proteins.")

# Helpers ----
log_line <- function(...) {
  write.table(
    paste0(...),
    file = log_file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, append = TRUE
  )
}

choose_tier <- function(ins, exp.dat) {
  # Returns list(dat=..., tier=...)
  # Tier 1: cis, p<1e-8, NConsistent>1
  ins2 <- ins[ins$pval.exposure < 1e-8 & ins$NConsistent > 1, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 1
    return(ins2)
  }

  # Tier 2: cis, p<1e-8
  ins2 <- ins[ins$pval.exposure < 1e-8, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 2
    return(ins2)
  }

  # Tier 3: strong trans p<1e-20 (from cis-subset 'ins' per your original)
  ins2 <- ins[ins$pval.exposure < 1e-20, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 3
    return(ins2)
  }

  # Tier 4: softer cis p<1e-5
  ins2 <- ins[ins$pval.exposure < 1e-5, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 4
    return(ins2)
  }

  # Tier 5: trans p<1e-9, NConsistent>1
  ins2 <- exp.dat[exp.dat$pval.exposure < 1e-9 & exp.dat$NConsistent > 1, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 5
    return(ins2)
  }

  # Tier 6: trans p<1e-9
  ins2 <- exp.dat[exp.dat$pval.exposure < 1e-9, ]
  if (nrow(ins2) > 0) {
    ins2$Tier <- 6
    return(ins2)
  }

  return(NULL)
}

# Collect all outputs (3 LD thresholds)
INS <- vector("list", 3)

# Main loop ----
for (file in chosen.files) {
  message("Starting ", file)

  prot_Apt <- gsub("1.tbl", "", file)
  protApt <- gsub("_", ".", paste0("seq_", prot_Apt))
  protSymb <- featureTable$EntrezGeneSymbol[featureTable$AptName == protApt]

  exp_path <- file.path(in_preproc_dir, paste0(prot_Apt, ".rds"))
  if (!file.exists(exp_path)) {
    log_line(protSymb, " missing preprocessed file: ", exp_path)
    next()
  }

  exp.dat <- readRDS(exp_path)

  # Keep only SNPs in the reference panel (panel.ids are chr:pos for rs variants in EUR.bim)
  exp.dat <- exp.dat[exp.dat$SNP %in% panel.ids, ]
  if (nrow(exp.dat) == 0) {
    log_line(protSymb, " does not have potential instruments in reference panel")
    next()
  }

  # CIS region selection: +/- 500kb around TSS
  tss <- ensembl.tss$TSS[ensembl.tss$AptName == protApt]
  cis_chr <- ensembl.tss$CHR[ensembl.tss$AptName == protApt]

  if (length(tss) == 0 || length(cis_chr) == 0) {
    log_line(protSymb, " has not been annotated to TSS")
    next()
  }

  exp.dat$TSS <- tss[1]
  exp.dat$cisChromosome <- cis_chr[1]

  ins <- exp.dat[exp.dat$chr.exposure == exp.dat$cisChromosome & abs(exp.dat$pos.exposure - exp.dat$TSS) < 5e5, ]

  # Tier selection (your logic, made deterministic)
  ins2 <- choose_tier(ins, exp.dat)
  if (is.null(ins2) || nrow(ins2) == 0) {
    log_line(protSymb, " does not have neither cis nor trans instruments")
    next()
  }

  # Annotate rsID using dbSNP + bcftools ----
  # Note: in your code SNP is chr:pos (panel.ids). ins2$chr.exposure + pos.exposure used for querying dbSNP VCF.
  chrom.aux <- chroms$`RefSeq-Accn`[match(as.character(ins2$chr.exposure), chroms$`Assigned-Molecule`)]

  fwrite(
    x = cbind(chrom.aux, ins2$pos.exposure),
    file = fn_positions,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  # bcftools query
  cmd <- paste0(
    "bcftools view ", shQuote(dbsnp_vcf),
    " -H -R ", shQuote(fn_positions),
    " -V indels > ", shQuote(fn_bcf_out)
  )
  system(cmd)

  aux <- tryCatch(
    fread(fn_bcf_out, header = FALSE, data.table = FALSE),
    error = function(e) NULL
  )

  if (is.null(aux) || nrow(aux) == 0) {
    # If no dbSNP match, keep going but rsid will be NA
    ins2$rsid <- NA_character_
  } else {
    aux$Marker <- paste0(chroms$`Assigned-Molecule`[match(aux$V1, chroms$`RefSeq-Accn`)], ":", aux$V2)
    ins2$rsid <- aux$V3[match(ins2$SNP, aux$Marker)]
  }

  ins2$pval <- ins2$pval.exposure

  # Keep only those with an rsid for clumping (ieugwasr expects rsid)
  ins2 <- ins2[!is.na(ins2$rsid), ]
  if (nrow(ins2) == 0) {
    log_line(protSymb, " has no rsIDs after dbSNP annotation; skipping")
    next()
  }

  # Annotate PAVs ----
  # Only for Tier <= 2, i.e. strong cis instruments
  if (unique(ins2$Tier) > 2) {
    ins2$PAV <- FALSE
  } else {
    ins2 <- ins2[order(ins2$chr.exposure, ins2$pos.exposure), ]
    ins2$PAV <- FALSE

    # Prepare VEP input
    vep_dir <- file.path(path.expand("~"), "vep_data")
    dir.create(vep_dir, showWarnings = FALSE, recursive = TRUE)

    tfn <- paste0("input_temp_", chunk)
    tfn2 <- paste0("output_temp_", chunk)

    x <- cbind(
      ins2$chr.exposure, ins2$pos.exposure, ".",
      ins2$other_allele.exposure, ins2$effect_allele.exposure,
      ".", ".", "."
    )
    write.table(x, file = file.path(vep_dir, tfn), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

    # Run VEP via docker (kept as-is, only made paths explicit)
    if (vep_mode == "apptainer") {
      cmd <- paste(
        "apptainer exec docker://ensemblorg/ensembl-vep",
        "vep --cache --offline",
        "-i", file.path(vep_dir, tfn),
        "-o", file.path(vep_dir, tfn2),
        "--force_overwrite --tab"
      )
    } else if (vep_mode == "docker") {
      cmd <- paste(
        "docker run -v", paste0(vep_dir, ":/data:Z"),
        "ensemblorg/ensembl-vep",
        "vep --cache --offline",
        "-i /data/", tfn,
        "-o /data/", tfn2,
        "--force_overwrite --tab"
      )
    } else if (vep_mode == "local") {
      cmd <- paste(
        "vep --cache --offline",
        "-i", file.path(vep_dir, tfn),
        "-o", file.path(vep_dir, tfn2),
        "--force_overwrite --tab"
      )
    } else {
      cmd <- NULL
    }
    if (!is.null(cmd)) {
      system(cmd)
    }

    res <- tryCatch(
      fread(file.path(vep_dir, tfn2), header = FALSE, data.table = FALSE),
      error = function(e) NULL
    )

    if (!is.null(res) && nrow(res) > 0) {
      # Keep only variants mapping to the target ENSID for this aptamer
      ensid <- dict$ENSID[dict$AptName == protApt]
      if (length(ensid) > 0) {
        res <- res[res$V4 == ensid[1], , drop = FALSE]
      }

      if (nrow(res) > 0) {
        res$PAV <- grepl(
          "(frameshift_variant)|(splice_acceptor_variant)|(splice_donor_variant)|(stop_gained)|(missense_variant)|(protein_altering_variant)|(splice_region_variant)|(inframe_deletion)|(inframe_insertion)",
          res$V7
        )

        # Flag PAVs among ins2
        for (i in seq_len(nrow(ins2))) {
          id <- ins2$SNP[i] # SNP is chr:pos
          if (id %in% res$V2) {
            if (TRUE %in% res$PAV[res$V2 == id]) ins2$PAV[i] <- TRUE
          }
        }

        # PAVs in LD with other instruments (if possible)
        if (sum(ins2$rsid %in% eur.bim$V2) > 1 && any(ins2$PAV)) {
          ins_ld <- ld_matrix(
            variants = ins2$rsid,
            plink_bin = plink_bin,
            bfile = plink_ref_prefix
          )

          for (id in ins2$rsid[ins2$PAV]) {
            abs.lds <- abs(ins_ld[grep(id, rownames(ins_ld)), ]) > 0.8
            ids <- gsub("_.*", "", names(abs.lds[abs.lds]))
            ins2$PAV[ins2$SNP %in% ids] <- TRUE
          }
        }
      }
    }
  }

  # LD clumping with different thresholds ----
  r2 <- c(0.001, 0.1, 0.2)
  ins3 <- vector(mode = "list", length = length(r2))

  for (i in seq_along(r2)) {
    ins3[[i]] <- ld_clump(
      dat = ins2,
      clump_p = 1,
      clump_r2 = r2[i],
      clump_kb = 250000,
      plink_bin = plink_bin,
      bfile = plink_ref_prefix
    )

    # If removing PAVs still leaves enough instruments, remove them
    if ("PAV" %in% names(ins3[[i]]) && sum(ins3[[i]]$PAV, na.rm = TRUE) > 0) {
      if (sum(ins3[[i]]$PAV, na.rm = TRUE) < nrow(ins3[[i]]) - 5) {
        ins3[[i]] <- ins3[[i]][ins3[[i]]$PAV == FALSE, ]
      }
    }
  }

  saveRDS(ins3, file = file.path(out_dir, paste0(prot_Apt, ".rds")))

  toc_prot <- Sys.time()
  message("Time for ", prot_Apt, ": ", as.numeric(difftime(toc_prot, tic, units = "mins")))

  # Add to global pile (per r2)
  for (i in seq_along(ins3)) {
    if (is.null(INS[[i]])) {
      INS[[i]] <- ins3[[i]]
    } else {
      INS[[i]] <- rbind(INS[[i]], ins3[[i]])
    }
  }
}

# Save aggregated
saveRDS(INS, file = file.path(out_dir, paste0("ALL_full", chunk, ".rds")))

# README + JSON metadata ----
toc <- Sys.time()
temps <- as.numeric(difftime(toc, tic, units = "mins"))

readme(
  path = file.path(root_data_dir, proj_dir, "README"),
  subanalisi = 2,
  titol = "Instrument selection",
  autors = "Sergio",
  lloc = "02_instrumentSelection.R",
  temps = temps,
  descripcio = "We consider tiers of instruments based on whether there are genome-wide significant signals in cis or only trans. 3 LD thresholds are considered for clumping. PAVs are annotated using ensembl-vep command line tool and removed if possible while retaining enough instruments.",
  altres = "This step is computationally demanding; run in parallel with chunk argument (1:8)."
)

metadata <- list(
  script = "02_instrumentSelection.R",
  chunk = chunk,
  timestamp = as.character(Sys.time()),
  n_files_total = length(files),
  n_files_in_chunk = length(chosen.files),
  R_version = R.version.string,
  plink_bin = plink_bin
)

write_json(
  metadata,
  path = file.path(out_dir, paste0("run_metadata_chunk_", chunk, ".json")),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Step 02 completed successfully for chunk ", chunk, ".")

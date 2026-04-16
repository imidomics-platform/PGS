library(httr)
library(jsonlite)
library(sodium)
library(dplyr)
library(ggplot2)
library(readxl)
library(purrr)
library(data.table)
library(gridExtra)
library(grid)
library(ggpubr)
library(caret)
library(pROC)
library(fmsb)

#---------- Base URL configuration ---------------#
#api_host <- Sys.getenv("DB_API_HOST", "localhost")
api_host <- Sys.getenv("DB_API_HOST", "10.7.50.21")
api_port <- Sys.getenv("DB_API_PORT", "8001")
BASE_URL <- Sys.getenv("DB_API_URL", paste0("http://", api_host, ":", api_port))

#---------------- Authentication -----------------#
# For testing purposes, we use a fixed admin user ID in helpers.R (admin_user_id) to set the X-User-Id header for all API calls. 
# This allows us to bypass authentication and role checks during testing.
# In a real scenario, you would implement proper authentication and token management.
# admin_user_id <- Sys.getenv("DB_API_ADMIN_USER_ID", "900c9b55-b976-4954-ad60-434d4239d538") # AK
admin_user_id <- Sys.getenv("DB_API_ADMIN_USER_ID", "2122b2ed-86f7-4478-a442-c5b258b5b5fe") # AA

#----------------- Load helper and endpoint functions -----------------#
# Load helper functions and endpoint functions
db_functions_dir <- "../imidomics_platform/db_service/tests/R/" # replace with actual R path if different
source(file.path(db_functions_dir, "helpers.R"))
source(file.path(db_functions_dir, "datasets.R"))
source(file.path(db_functions_dir, "users.R"))
source(file.path(db_functions_dir, "projects.R"))
source(file.path(db_functions_dir, "users_projects.R"))
source(file.path(db_functions_dir, "reference_sources.R"))
source(file.path(db_functions_dir, "db.R"))

#----------------- Set up volume directory -----------------#
# volume directory configuration for file-based reference sources
volume_dir <- Sys.getenv("VOLUME_DIR", unset = "/media/bioinformatics/imidomics_platform/volume") # replace with actual volume path if different

# ===================================
# Accessing database for PGS analysis 
# ===================================

disease_of_interest <- "SS"
tissue <- "salivary_gland"

assay <- "singlecell"
project <- "SSAD"

analysis <- "a3"

excl.lowexp<-"no"

pseudo_bulk<-"yes"

# 1) Getting bulk transcriptomics dataset for the specified disease, assay type, tissue and project

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project, assay=assay, tissue=tissue), limit = 50)
dataset_path    <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

if (pseudo_bulk=="yes") {
  dataset_data_mol <- readRDS(paste0(dataset_path, "/", dataset_yaml[[grep("pseudo_bulk_per_sample_data",names(dataset_yaml),value=T)[1]]]))
} else {
  dataset_data_mol <- readRDS(paste0(dataset_path, "/", dataset_yaml[[grep("_data",names(dataset_yaml),value=T)[1]]]))
}
dataset_samples_mol<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_tech_mol<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

print(paste0(disease_of_interest," | ",assay2," | ",project2," | ",tissue," | ",dataset_info$dataset_name))

# 2) Getting reference files

gene.anno<-read.csv(paste0(volume_dir,"/reference/gene_reference_data/gene_annotation_hg38v44_hg19v19.tsv"),sep="\t")
transc.feat<-readRDS("/media/bioinformatics/Projects/Internal_Projects/2024_IMX_PlatformUpdate/InputData/Internal/GenExp_hg38_featureTable.rds") # project-specific feature table
#transc.feat<-readRDS(paste0(volume_dir,"/reference/gene_reference_data/hg38_GenExp_featureTable.rds"))

# =================================================
# DGEA
# =================================================

tmp.folder<-paste0(paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/analyses/",analysis))
dir.create(tmp.folder)

#---------------- Merging Datasets -------------#

if (tissue=="skin" | tissue=="salivary_gland") {
  df.mol.cases<-merge(dataset_samples_mol[(dataset_samples_mol$timepoint=="wk0") & dataset_samples_mol$tissue_type=="affected", c("sex","age","timepoint","sample_id","followup_id","disease")], dataset_tech_mol[,c("sample_id","Sequencing_batch")], by="sample_id")
  df.mol.ctrls<-merge(dataset_samples_mol[(dataset_samples_mol$timepoint=="uv") & dataset_samples_mol$tissue_type=="unaffected", c("sex","age","timepoint","sample_id","followup_id","disease")], dataset_tech_mol[,c("sample_id","Sequencing_batch")], by="sample_id")
  df.mol.samples<-rbind(df.mol.cases,df.mol.ctrls)
} else {
  df.mol.samples<-merge(dataset_samples_mol[(dataset_samples_mol$timepoint=="wk0" | dataset_samples_mol$timepoint=="uv"),c("sex","age","timepoint","sample_id","followup_id","disease")], dataset_tech_mol[,c("sample_id","Sequencing_batch")], by="sample_id")
}

mol.data<-as.data.frame(t(dataset_data_mol))
if (excl.lowexp=="yes") {genes2analyze<-rownames(subset(transc.feat,ExpressedPercentageCategory_Disc != "Non-expressed" & ExpressedPercentageCategory_Repl != "Non-expressed"))} else {genes2analyze<-grep("ENS",colnames(mol.data),value=T)}
mol.data.filt<-mol.data[colnames(mol.data) %in% genes2analyze]

d2a<-merge(df.mol.samples, mol.data.filt, by.x="sample_id", by.y=0)
d2a$Sequencing_batch<-as.factor(d2a$Sequencing_batch)
d2a$sex<-as.factor(d2a$sex)
d2a$IMID.num <- ifelse(d2a$disease == disease_of_interest, 1, 0)


#---------------- Association Analysis ----#

results <- lapply(genes2analyze, function(gene) {
  formula <- as.formula(paste0("IMID.num ~ `", gene, "` + age + sex + Sequencing_batch"))
  model <- glm(formula, data = d2a, family = "binomial")
  
  coefs <- coefficients(summary(model))
  
  if (!(gene %in% rownames(coefs))) {
    message("Missing coefficient for: ", gene)
    return(data.frame(
      gene = gene,
      Estimate = NA_real_,
      Std.Error = NA_real_,
      P.value = NA_real_
    ))
  }
  
  data.frame(
    gene = gene,
    Estimate = coefs[gene, "Estimate"],
    Std.Error = coefs[gene, "Std. Error"],
    P.value = coefs[gene, "Pr(>|z|)"]
  )
})
res <- do.call(rbind, results)

colnames(res)<-c("EnsemblID","Effect","Std.Error","Pvalue")
res$Effect<-as.numeric(res$Effect); res$Pvalue<-as.numeric(res$Pvalue); res$Std.Error<-as.numeric(res$Std.Error);  res$TechID<-NA
res<-merge(res,gene.anno[,c("EnsemblID","gene_symbol")],by="EnsemblID",all.x=T)
res<-res[,c(1,6,5,2:4)]
res<-res[order(res$Pvalue),]
res$FDR<-p.adjust(res$Pvalue, method="fdr")
if (pseudo_bulk=="yes") {assay2<-"transcriptomics"} 
saveRDS(res,paste0(tmp.folder,"/",disease_of_interest,"_BulkTranscriptomics_",assay2,"_",tissue,".rds"))
dim(res)

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

disease_of_interest <- "AD"
tissue <- "blood"

assay1 <- "genetics"
project1 <- "ssad_expanded"

assay2 <- "transcriptomics"
project2 <- "ssad"

analysis <- "a1"

excl.lowexp<-"no"

# 1) Getting pgs dataset & genetic pc information for the specified disease, assay type, tissue and project

dataset_path <- list.files(paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/analyses/",analysis), pattern = "hg", full.names = TRUE)
ths.all<-readRDS(paste0(dataset_path,"/PGSprofiles/ThresholdingSummary.rds"))
ths.opt<-ths.all[!is.na(ths.all$Nagelkerke.R2.All.FullDataset),"Threshold"]

dataset_genetics_pgs <- read.table(paste0(dataset_path,"/PGSprofiles/",sub("^LD0?(\\d)\\.(P\\d+)$","PGS_R2_0.\\1.\\2.profile",ths.opt)),header=T)

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project1, assay = assay1), limit = 50)
dataset_path    <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_genetics_pcs<-read.delim(paste0(dataset_path, "/", dataset_yaml$pc_data),header=T)
dataset_genetics_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)

print(paste0(disease_of_interest," | ",assay1," | ",project1," | ",dataset_info$dataset_name))

# 2) Getting molecular dataset information for the specified disease, assay type, tissue and project

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project2, assay=assay2, tissue=tissue), limit = 50)
dataset_path    <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_data_mol <- readRDS(paste0(dataset_path, "/", dataset_yaml[[grep("_data",names(dataset_yaml),value=T)[1]]]))
dataset_samples_mol<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_tech_mol<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

print(paste0(disease_of_interest," | ",assay2," | ",project2," | ",tissue," | ",dataset_info$dataset_name))

# 3) Getting reference files

gene.anno<-read.csv(paste0(volume_dir,"/reference/gene_reference_data/gene_annotation_hg38v44_hg19v19.tsv"),sep="\t")
transc.feat<-readRDS("/media/bioinformatics/Projects/Internal_Projects/2024_IMX_PlatformUpdate/InputData/Internal/GenExp_hg38_featureTable.rds") # project-specific feature table
#transc.feat<-readRDS(paste0(volume_dir,"/reference/gene_reference_data/hg38_GenExp_featureTable.rds"))


# =================================================
# mPGS Calculation 
# =================================================

tmp.folder<-paste0(paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/analyses/",analysis))

#---------------- Merging Datasets -------------#

df.gen.v0<-merge(dataset_genetics_pgs[,c("IID","SCORE")],dataset_genetics_samples[,c("sample_id","followup_id")], by.x="IID", by.y="sample_id")
df.gen.v1<-merge(df.gen.v0,dataset_genetics_pcs[,c("sample_id","PC1","PC2")], by.x="IID", by.y="sample_id")
df.gen<-df.gen.v1[,c(3:5,2)]
colnames(df.gen)<-c("followup_id","gen.pc1","gen.pc2","pgs")

df.mol.samples<-merge(dataset_samples_mol[dataset_samples_mol$timepoint=="WK0" & dataset_samples_mol$disease==disease_of_interest,c("sex","age","sample_id","followup_id")], dataset_tech_mol[,c("sample_id","Sequencing_batch")], by="sample_id")
df.all.samples<-merge(df.mol.samples,df.gen,by="followup_id")

mol.data<-as.data.frame(t(dataset_data_mol))
if (excl.lowexp=="yes") {genes2analyze<-rownames(subset(transc.feat,ExpressedPercentageCategory_Disc != "Non-expressed" & ExpressedPercentageCategory_Repl != "Non-expressed"))} else {genes2analyze<-grep("ENS",colnames(mol.data),value=T)}
mol.data.filt<-mol.data[colnames(mol.data) %in% genes2analyze]

d2a<-merge(df.all.samples, mol.data.filt, by.x="sample_id", by.y=0)
d2a$pgs<-as.numeric(scale(d2a$pgs))
d2a$Sequencing_batch<-as.factor(d2a$Sequencing_batch)
d2a$sex<-as.factor(d2a$sex)

#---------------- Association Analysis ----#

results <- lapply(genes2analyze, function(gene) {
  formula <- as.formula(paste(gene, "~ pgs + gen.pc1 + gen.pc2 + age + sex + Sequencing_batch"))
  model <- lm(formula, data=d2a)
  out <- c(gene,coefficients(summary(model))["pgs","Estimate"],coefficients(summary(model))["pgs","Std. Error"],coefficients(summary(model))["pgs","Pr(>|t|)"])
})

res<-as.data.frame(do.call(rbind, results))
colnames(res)<-c("EnsemblID","Effect","Std.Error","Pvalue")
res$Effect<-as.numeric(res$Effect); res$Pvalue<-as.numeric(res$Pvalue); res$Std.Error<-as.numeric(res$Std.Error);  res$TechID<-NA
res<-merge(res,gene.anno[,c("EnsemblID","gene_symbol")],by="EnsemblID",all.x=T)
res<-res[,c(1,6,5,2:4)]
res<-res[order(res$Pvalue),]
res$FDR<-p.adjust(res$Pvalue, method="fdr")
saveRDS(res,paste0(tmp.folder,"/PGS_",disease_of_interest,"_",assay2,"_",tissue,".rds"))
dim(res)

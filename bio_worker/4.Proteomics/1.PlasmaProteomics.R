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
library(microbenchmark)
library(peakRAM)

########  Configuration  ########

#### Load configuration ####
suppressWarnings(local_config <- yaml::read_yaml("config.yml"))

# directories
volume_dir <- local_config$global$volume_dir
tmp_dir <- local_config$global$tmp_dir
output_dir <- paste0(tmp_dir,"TargetID/")
if(!dir.exists(output_dir)){dir.create(output_dir)}

# database API
api_host <- local_config$database_api$host
api_port <- local_config$database_api$port
BASE_URL <- paste0("http://", api_host, ":", api_port)

#### Load helper and endpoint functions ####

# Load helper functions and endpoint functions
db_functions_dir <- local_config$database_wrapper$R_dir

# source all .R files in the db_functions_dir
DB_R_files <- list.files(db_functions_dir, pattern = "*.R", full.names=T)
temp <- lapply(DB_R_files, source)

#### Set User ID ####

db_user_name <- "aaterido"

user_ids <- db_list_users()
admin_user_id <- user_ids$user_id[user_ids$user_name == db_user_name]

source(paste0(db_functions_dir, "/helpers.R"))

########  Analysis Workflow  ########

#### Input ####

disease_of_interest <- "UC"

assay <- "proteomics"
project <- "Cross-sectional"
tissue <- "plasma"

analysis_design <- "case-control"
analysis_goal <- "TargetID"
analysis_scope<-"PlasmaProteomics"

analysis_method <- paste0("LM: ",assay)

out.path<-paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds")
analysis.path<-gsub("/home/aaterido/Escritorio/GitHub_Repositories/","/",paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/",basename(rstudioapi::getSourceEditorContext()$path)))

#### Getting proteomics dataset information for the specified input ####

dataset_info <- db_list_datasets(list(disease = disease_of_interest, dataset_project=project, assay=assay, tissue=tissue))
dataset_id <- dataset_info$dataset_id
dataset_path    <- paste0(volume_dir,dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_samples_mol<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_samples_tech<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

dataset_analysis_data <- readRDS(paste0(dataset_path, "/", dataset_yaml$analysis_data))
dataset_analysis_covs <- unlist(strsplit(dataset_yaml$analysis_covariates,", "))

if (dataset_yaml$pcs_data!="") {
  dataset_analysis_pcs <-  read.delim(paste0(dataset_path, "/",dataset_yaml$pcs_data),header=T)
} else {
  dataset_analysis_pcs <- ""
}

print(paste0(disease_of_interest," | ",assay," | ",project," | ",tissue," | ",dataset_info$dataset_name))

#### Getting reference information ####

tech.source<-dataset_info$technology_platform
prot.anno_reference_id <- db_list_reference_sources(list(source = tech.source, is_default = TRUE))$reference_id
prot.anno_reference <- db_get_reference_source(prot.anno_reference_id)
corresp<-read.csv(paste0(volume_dir,prot.anno_reference$reference_path),sep="\t")

if (tech.source=="SomaScan_Assay_v4.1") {
    corresp<-corresp[,c("AptName","Ensembl_ID","Updated_Symbol")]
}
if (tech.source=="Olink Explore HT 7.0.5") {
    corresp<-corresp[,c("OlinkID","Ensembl_ID","MappedGene")]
}
colnames(corresp)<-c("TechID","Ensembl_ID","Symbol")

#### Association Analysis ####

# Getting Proteomics Data (Filtering by wk0) & Analysis Covariates

id.merge<-colnames(dataset_samples_mol)[colnames(dataset_samples_mol) %in% colnames(dataset_samples_tech)]

df.mol.samples<-dataset_samples_mol[dataset_samples_mol$timepoint=="wk0",c("sex","age","disease",id.merge)]

id.merge2<-colnames(dataset_samples_tech)[colnames(dataset_samples_tech) %in% colnames(dataset_analysis_pcs)]
dataset_vars<-merge(dataset_samples_tech, dataset_analysis_pcs, by=id.merge2)

analysis.covs2get<-colnames(dataset_vars)[colnames(dataset_vars) %in% dataset_analysis_covs]

id.merge3<-colnames(df.mol.samples)[colnames(df.mol.samples) %in% colnames(dataset_vars)]
df.prot<-merge(df.mol.samples,dataset_vars[,c(id.merge3, analysis.covs2get)], by=id.merge3)

rownames(df.prot) <- df.prot[[1]]
df.prot.samples <- df.prot[ , -1]
df.prot<-merge(df.prot.samples,t(dataset_analysis_data),by=0)
rownames(df.prot) <- df.prot[[1]]
df.prot <- df.prot[ , -1]
rownames(df.prot) <- sub("^(([^-]+-){2}[^-]+)-[^-]+$", "\\1", rownames(df.prot))

# Merging omics layers

df.prot$sex<-as.factor(df.prot$sex)
df.prot$age<-as.numeric(df.prot$age)

df.prot$IMID <- factor(ifelse(df.prot$disease == disease_of_interest, "1", "0"))

dim(df.prot)

covs2a<-colnames(df.prot)[colnames(df.prot) %in% c(dataset_analysis_covs)]

print(table(df.prot$disease))
print(covs2a)
print(dataset_analysis_covs)

# Statistical Analysis

mol2analyze<-colnames(df.prot)[!(colnames(df.prot) %in% c("disease","IMID",covs2a))]

results <- lapply(mol2analyze, function(mol) {
  formula <- as.formula(paste0("IMID ~ ","`", mol, "` + ",paste(covs2a,collapse=" + ")))
  model <- glm(formula, data=df.prot, family="binomial")
  coef.model<-rownames(coefficients(summary(model)))
  if (mol %in% coef.model) {
    out <- c(mol,coefficients(summary(model))[mol,"Estimate"],coefficients(summary(model))[mol,"Std. Error"],coefficients(summary(model))[mol,"Pr(>|z|)"])
  } else {
    out <- c(mol,NA,NA,NA)
  }
})

res<-as.data.frame(do.call(rbind, results))
colnames(res)<-c("TechID","Effect","Std.Error","Pvalue")
res$Effect<-as.numeric(res$Effect); res$Pvalue<-as.numeric(res$Pvalue); res$Std.Error<-as.numeric(res$Std.Error);
res<-merge(res,corresp,by="TechID",all.x=T)
res<-res[order(res$Pvalue),]
res<-res[,c(5,1,6,2,4)]
colnames(res)<-c("EnsemblID","TechID","Symbol","Effect","Pvalue")

saveRDS(res,out.path)

#### Register the analysis: create a new analysis in the database  ####

# Payload fields: http://10.7.50.21:8000/docs#/analyses/create_analysis_analyses__post

new_analysis_payload <- list(
  analysis_goal = analysis_goal,
  analysis_scope = analysis_scope,
  analysis_type = "",
  analysis_subtype = paste0(assay,"_",tissue),
  analysis_method = analysis_method,
  analysis_design = analysis_design,
  parameters = list("Data.Molecular"=dataset_yaml$analysis_data,"Data.Covariates"=paste(covs2a,collapse=" + ")),
  omics_layer = assay,
  disease_name = disease_of_interest,
  datasets = list(dataset_id),
  reference_sources = as.list(c(prot.anno_reference_id)),
  input_files = list()
)

new_analysis <- db_create_analysis(new_analysis_payload)

update_analysis_payload <- list(
  analysis_path = out.path,
  analysis_code = analysis.path
)

db_update_analysis(new_analysis$analysis_id, update_analysis_payload)

print(table(df.prot$disease))
print(covs2a)
print(dataset_analysis_covs)

View(db_list_analyses()[,c("analysis_id","analysis_path","disease_name")])







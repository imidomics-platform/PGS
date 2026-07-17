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
DB_R_files <- list.files(db_functions_dir, pattern = "*.R", full.names = TRUE)
temp <- lapply(DB_R_files, source)

#### Set User ID ####

db_user_name <- "aaterido"

user_ids <- db_list_users()
admin_user_id <- user_ids$user_id[user_ids$user_name == db_user_name]

source(paste0(db_functions_dir, "/helpers.R"))


########  Analysis Workflow  ########

#### Input ####

disease_of_interest <- "PS"

assay1 <- "genetics"
project1 <- "Cross-sectional"

assay2 <- "proteomics"
project2 <- "Cross-sectional"

analysis_design <- "case-only"
analysis_goal <- "TargetID"
analysis_method <- paste0("LM: ",assay2," ~ ",assay1)
if (assay2=="transcriptomics") {analysis_scope<-"tPGS"; tissue <- "blood"}
if (assay2=="proteomics") {analysis_scope<-"pPGS"; tissue <- "plasma"}

#### Getting genetic dataset information for the specified input ####

df.a<-as.data.frame(db_list_analyses())
params<-subset(df.a, analysis_goal == "PGS" & disease_name == disease_of_interest & !is.na(parameters), select = "parameters")[[1]]
ths.opt<-fromJSON(params)$Data.Optimal$Threshold
input.dir<-subset(df.a, analysis_goal == "PGS" & disease_name == disease_of_interest & !is.na(parameters), select = "analysis_path")[[1]]

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project1, assay = assay1))
dataset_path <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

if (grepl("PGS", ths.opt)) {
  dataset_pgs_cat <- paste0(dataset_path, "/", dataset_yaml$pgs.cat_data)
  dataset_genetics_pgs<-read.delim(dataset_pgs_cat, header = TRUE)[,c("sample_id",ths.opt)]
  colnames(dataset_genetics_pgs)[2]<-"SCORE"
} else {
  dataset_genetics_pgs <- read.table(paste0(input.dir,sub("^LD0?(\\d)\\.(P\\d+)$","PGS_R2_0.\\1.\\2.profile",ths.opt)),header=T)
  colnames(dataset_genetics_pgs)[2]<-"sample_id"
}

dataset_genetics_pcs<-read.delim(paste0(dataset_path, "/", dataset_yaml$pc_data),header=T)
dataset_genetics_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)

print(paste0(disease_of_interest," | ",assay1," | ",project1," | ",dataset_info$dataset_name))

#### Getting molecular dataset information for the specified input ####

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project2, assay=assay2, tissue=tissue))
dataset_id <- dataset_info$dataset_id
dataset_path    <- paste0(volume_dir,dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_samples_mol<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_samples_tech<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

dataset_analysis_data <- readRDS(paste0(dataset_path, "/", dataset_yaml$analysis_data))
#dataset_analysis_covs <- unlist(strsplit(dataset_yaml$analysis_covariates,", "))

print(paste0(disease_of_interest," | ",assay2," | ",project2," | ",tissue," | ",dataset_info$dataset_name))

#### Getting reference information ####

gene.anno_reference_id <- db_list_reference_sources(list(reference_type = "annotation", reference_name="Gene Annotation", is_default = TRUE))$reference_id
gene.anno_reference <- db_get_reference_source(gene.anno_reference_id)
corresp<-read.csv(paste0(volume_dir,gene.anno_reference$reference_path),sep="\t")

#### Association Analysis ####

# Merging Datasets

df.gen<-merge(dataset_genetics_pgs[,c("sample_id","SCORE")],dataset_genetics_pcs[,c("sample_id","PC1","PC2")], by="sample_id")
colnames(df.gen)<-c("followup_id","pgs","gen.pc1","gen.pc2")

if (tissue=="skin" | tissue=="salivary_gland") {
  df.mol.samples<-merge(dataset_samples_mol[dataset_samples_mol$timepoint=="wk0" & dataset_samples_mol$tissue_type=="affected" & dataset_samples_mol$disease==disease_of_interest,c("sex","age","timepoint","sample_id","followup_id")], dataset_samples_tech[,c("sample_id","Sequencing_batch")], by="sample_id")
} else {
  df.mol.samples<-merge(dataset_samples_mol[dataset_samples_mol$timepoint=="wk0" & dataset_samples_mol$disease==disease_of_interest,c("sex","age","timepoint","sample_id","followup_id")], dataset_samples_tech[,c("sample_id",grep("age|sex",unlist(strsplit(dataset_yaml$analysis_covariates,", ")),invert=T,value=T))], by="sample_id")
}
df.all.samples<-merge(df.mol.samples,df.gen,by.x="sample_id",by.y="followup_id")

d2a<-merge(df.all.samples, as.data.frame(t(dataset_analysis_data)), by.x="sample_id", by.y=0)
d2a$pgs<-as.numeric(scale(d2a$pgs))
d2a[[grep("age|sex",unlist(strsplit(dataset_yaml$analysis_covariates,", ")),invert=T,value=T)]]<-as.factor(d2a[[grep("age|sex",unlist(strsplit(dataset_yaml$analysis_covariates,", ")),invert=T,value=T)]])
d2a$sex<-as.factor(d2a$sex)
dim(d2a)

# Statistical Analysis

genes2analyze<-grep("^ENSG", colnames(d2a), value=T)

results <- lapply(genes2analyze, function(gene) {
  formula <- as.formula(paste0("`", gene, "` ~ pgs + ",gsub(","," +",paste0(dataset_yaml$analysis_covariates,", gen.pc1, gen.pc2"))))
  model<- lm(formula, data = d2a)
  out <- c(gene,coefficients(summary(model))["pgs","Estimate"],coefficients(summary(model))["pgs","Std. Error"],coefficients(summary(model))["pgs","Pr(>|t|)"])
})

res<-as.data.frame(do.call(rbind, results))
colnames(res)<-c("EnsemblID","Effect","Std.Error","Pvalue")
res$Effect<-as.numeric(res$Effect); res$Pvalue<-as.numeric(res$Pvalue); res$Std.Error<-as.numeric(res$Std.Error);  res$TechID<-NA
res<-merge(res,corresp[,c("EnsemblID","gene_symbol")],by="EnsemblID",all.x=T)
res<-res[,c(1,6,5,2:4)]
res<-res[order(res$Pvalue),]
res$FDR<-p.adjust(res$Pvalue, method="fdr")

saveRDS(res,paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds"))

####  Register the analysis: create a new analysis in the database  ####

# Payload fields: http://10.7.50.21:8000/docs#/analyses/create_analysis_analyses__post

new_analysis_payload <- list(
  analysis_goal = analysis_goal,
  analysis_scope = analysis_scope,
  analysis_type = analysis_type,
  analysis_subtype = paste0(assay1,",",assay2,"_",tissue),
  analysis_method = analysis_method,
  analysis_design = analysis_design,
  parameters = list("Data.Molecular"=dataset_yaml$analysis_data,"Data.Covariates"=paste0(dataset_yaml$analysis_covariates,", gen.pc1, gen.pc2"),"Data.Filter"=paste0("PGS:",ths.opt)),
  omics_layer = paste0(assay1,",",assay2),
  disease_name = disease_of_interest,
  datasets = list(dataset_id),
  reference_sources = as.list(c(gene.anno_reference_id)),
  input_files = list()
)

new_analysis <- db_create_analysis(new_analysis_payload)

update_analysis_payload <- list(
  analysis_path = paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds"),
  analysis_code = gsub("/home/aaterido/Escritorio/GitHub_Repositories/","/",paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/",basename(rstudioapi::getSourceEditorContext()$path)))
)

db_update_analysis(new_analysis$analysis_id, update_analysis_payload)







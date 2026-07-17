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

disease_of_interest <- "RA"
tissue <- "blood"
assay1 <- "transcriptomics"
project1 <- "Cross-sectional"
analysis_design <- "case-control"
analysis_goal <- "TargetID"
analysis_subtype <- "DGEA"
analysis_scope <- "BulkBlood"

#### Getting dataset information for the specified input ####

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project1, assay=assay1, tissue=tissue))
dataset_id <- dataset_info$dataset_id
dataset_path    <- paste0(volume_dir,dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_samples_tech<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

dataset_analysis_data <- readRDS(paste0(dataset_path, "/", dataset_yaml$analysis_data))
dataset_analysis_covs <- unlist(strsplit(dataset_yaml$analysis_covariates,", "))

print(paste0(disease_of_interest," | ",assay1," | ",project1," | ",tissue," | ",dataset_info$dataset_name))

#### Getting reference information ####

gene.anno_reference_id <- db_list_reference_sources(list(reference_type = "annotation", reference_name="Gene Annotation", is_default = TRUE))$reference_id
gene.anno_reference <- db_get_reference_source(gene.anno_reference_id)
corresp<-read.csv(paste0(volume_dir,gene.anno_reference$reference_path),sep="\t")

gene.featureTable_id <- db_list_reference_sources(list(source = project1, disease=disease_of_interest, is_default = TRUE))$reference_id
gene.featureTable_reference <- db_get_reference_source(gene.featureTable_id)
featureTable<-read.csv(paste0(volume_dir,gene.featureTable_reference$reference_path),sep="\t")

#### Differential Gene Expression Analysis ####

filt.genes.criterion<-c("Non-expressed","PartiallyExpressed")
filt.genes.out<-featureTable[(featureTable$ExpressedPercentageCategory %in% filt.genes.out),"feature_id"]
RNA<-t(dataset_analysis_data[!(rownames(dataset_analysis_data) %in% filt.genes.out),])

s2a<-merge(dataset_samples[dataset_samples$timepoint=="wk0",],dataset_samples_tech,by="sample_id"); rownames(s2a)<-s2a$sample_id
df2a<-merge(s2a[,c("disease",dataset_analysis_covs)],RNA,by=0)

df2a$IMID <- factor(ifelse(df2a$disease == disease_of_interest, "1", "0"))
genes<-grep("^ENSG", names(df2a), value=T)

covs.text<-paste(dataset_analysis_covs,collapse=" + ")

results <- lapply(genes, function(gene) {
  formula <- as.formula(paste0("IMID ~ ","`", gene, "` + ",covs.text))
  model <- glm(formula, data=df2a, family="binomial")
  coef.model<-rownames(coefficients(summary(model)))
  if (gene %in% coef.model) {
    out <- c(disease_of_interest,gene,coefficients(summary(model))[gene,"Estimate"],coefficients(summary(model))[gene,"Std. Error"],coefficients(summary(model))[gene,"Pr(>|z|)"])
  } else {
    out <- c(disease_of_interest,gene,NA,NA,NA)
  }
})

res<-as.data.frame(do.call(rbind, results))[,c(2,3,5)]
colnames(res)<-c("Ensembl_ID","Effect","Pvalue")
res$Effect<-as.numeric(res$Effect); res$Pvalue<-as.numeric(res$Pvalue)
res$TechID<-NA
res<-merge(res,corresp[,c("EnsemblID","gene_symbol")],by.x="Ensembl_ID",by.y="EnsemblID", all.x=T)
res<-res[,c(1,4,5,2,3)]
colnames(res)[3]<-"Symbol"
res<-res[order(res$Pvalue),]
res<-res[!is.na(res$Pvalue),]

saveRDS(res,paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds"))

####  Estimate computational resources  ####

set.seed(1)
sample_genes <- sample(genes, min(20, length(genes)))

# Benchmark Runtime
mb <- microbenchmark(lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial"); coefficients(summary(model))[gene, "Estimate"]}),times=5)
times_sec <- mb$time / 1e9

# CPU Timing
cpu_time <- system.time({tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})})

# Memory Usage (MB aprox)
gc_before <- gc()
tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})
gc_after <- gc()
mem_used_mb <- sum(gc_after[,2] - gc_before[,2])*8 / 1024

# Peak memory
peak <- peakRAM({tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})})
peak_ram_mb <- peak$Peak_RAM_Used_MiB

# Disk usage (save results)
disk_usage_mb <- file.size(paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds")) / (1024^2)

# Derived metrics
median_batch_time <- median(times_sec)
genes_tested <- length(sample_genes)
total_genes <- length(genes)

time_per_gene <- median_batch_time / genes_tested
estimated_total_time <- time_per_gene * total_genes
throughput <- genes_tested / median_batch_time

estimated_disk_total <- disk_usage_mb * (total_genes / genes_tested)

# Final table

resource_table <- data.frame(Metric = c("Median batch time (s)", "User CPU time (s)", "Elapsed time (s)", "Memory used (MB)", "Peak memory (MB)", "Disk usage sample (MB)", "Estimated disk usage total (MB)", "Genes per batch", "Total genes", "Time per gene (s)", "Throughput (genes/sec)", "Estimated total time (min)"),
                             Value = c(median_batch_time, cpu_time["user.self"], cpu_time["elapsed"], mem_used_mb, peak_ram_mb, disk_usage_mb, estimated_disk_total, genes_tested, total_genes, time_per_gene, throughput, estimated_total_time/60))


####  Register the analysis: create a new analysis in the database  ####

# Payload fields: http://10.7.50.21:8000/docs#/analyses/create_analysis_analyses__post

new_analysis_payload <- list(
  analysis_goal = analysis_goal,
  analysis_scope = analysis_scope,
  analysis_type = analysis_type,
  analysis_subtype = paste0(assay1,"_",tissue),
  analysis_method = analysis_method,
  analysis_design = analysis_design,
  parameters = list("Data.Molecular"=dataset_yaml$analysis_data,"Data.Covariates"=dataset_yaml$analysis_covariates,"Data.Filter"=paste(filt.genes.criterion,collapse=",")),
  omics_layer = assay1,
  disease_name = disease_of_interest,
  datasets = list(dataset_id),
  reference_sources = as.list(c(gene.anno_reference_id,gene.featureTable_id)),
  input_files = list()
)

new_analysis <- db_create_analysis(new_analysis_payload)

update_analysis_payload <- list(
  analysis_path = paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds"),
  analysis_code = gsub("/home/aaterido/Escritorio/GitHub_Repositories/","/",paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/",basename(rstudioapi::getSourceEditorContext()$path)))
)

db_update_analysis(new_analysis$analysis_id, update_analysis_payload)

db_list_analyses()

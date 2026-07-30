
########  Libraries  ########

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
library(limma)

########  Configuration  ########

#### Load configuration ####
suppressWarnings(local_config <- yaml::read_yaml("config.yml"))

# directories
volume_dir <- local_config$global$volume_dir
tmp_dir <- local_config$global$tmp_dir

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

# Datasets to be used according to AI-powered algorithm

# AD  - d106 - tissue_type, timepoint, technical_replicate - LimmaVoomDreamWeights - data/RNA.rds
# AD  - d107 - tissue_type, timepoint                      - LimmaVoomDreamWeights - data/RNA.rds
# AD  - d108 - tissue_type                                 - LimmaVoom             - data/RNA.rds                                  
# SS  - d109 - sex, age                                    - LimmaVoom             - data/RNA.rds
# SS  - d110 - sex, age, tissue_region                     - LimmaVoom             - data/RNA.rds
# CD  - d115 - sex, age                                    - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# CD  - d117 - sex, age                                    - LimmaVoom             - data/RNA.rds
# CD  - d118 - sex, age                                    - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# CD  - d124 -                                             - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# UC  - d116 - sex, age                                    - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# UC  - d125 -                                             - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# RA  - d119 - sex, age                                    - LimmaLinearModel      - preprocessed/RNA_NormNoAdj.rds
# PSA - d123 -                                             - LimmaVoom             - data/RNA.rds
# PS  - d122 -                                             - LimmaVoom             - data/RNA.rds 
# SLE - d114 - tissue_region, technical_replicate          - LimmaVoomDreamWeights - preprocessed/RNA_NormNoAdj.rds












disease_of_interest <- "RA"
tissue <- "blood"
assay1 <- "transcriptomics"
project1 <- "SSAD"
analysis_design <- "case-control"
analysis_goal <- "TargetID"
analysis_type <- "DGEA"
analysis_scope <- "BulkBlood"
analysis_method <- "LimmaLinearModel"                             # LimmaLinearModel vs LogisticRegression
analysis_covs <-c("sex","age","Sequencing_batch")                 # 6 IMIDs: c("sex","age","plateId"); SSAD: c("sex", "age", "Sequencing_batch")












analysis_subtype <- paste0(analysis_method,"_",paste(analysis_covs,collapse="+"))

output_dir <- paste0(tmp_dir,"/",analysis_goal,"/")
if(!dir.exists(output_dir)){dir.create(output_dir)}
analysis_path <- paste0(output_dir,analysis_scope,"_",disease_of_interest,"_",analysis_method,".rds")
analysis_code <- gsub("/home/aaterido/Escritorio/GitHub_Repositories/","/",paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/",basename(rstudioapi::getSourceEditorContext()$path)))

#### Get dataset information for the specified input ####

dataset_info <- db_list_datasets(list(disease = disease_of_interest, dataset_project=project1, assay=assay1, tissue=tissue))
dataset_id <- dataset_info$dataset_id
dataset_path    <- paste0(volume_dir,dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)

dataset_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_samples_tech<-read.delim(paste0(dataset_path, "/samples/samples_technology.tsv"),header=T)

dataset_analysis_data <- readRDS(paste0(dataset_path, "/", dataset_yaml$analysis_data))

print(paste0(disease_of_interest," | ",assay1," | ",project1," | ",tissue," | ",dataset_info$dataset_name))

#### Get reference information ####

gene.anno_reference_id <- db_list_reference_sources(list(reference_type = "annotation", reference_name="Gene Annotation", is_default = TRUE))$reference_id
gene.anno_reference <- db_get_reference_source(gene.anno_reference_id)
corresp<-read.csv(paste0(volume_dir,gene.anno_reference$reference_path),sep="\t")

gene.featureTable_id <- db_list_reference_sources(list(source = project1, disease=disease_of_interest, is_default = TRUE))$reference_id
gene.featureTable_reference <- db_get_reference_source(gene.featureTable_id)
featureTable<-read.csv(paste0(volume_dir,gene.featureTable_reference$reference_path),sep="\t")

#### Filter Data for Analysis ####

filt.genes.criterion<-c("Non-expressed","PartiallyExpressed")
filt.genes.out<-featureTable[(featureTable$ExpressedPercentageCategory %in% filt.genes.criterion),"feature_id"]
RNA<-t(dataset_analysis_data[!(rownames(dataset_analysis_data) %in% filt.genes.out),])

s2a<-merge(dataset_samples[dataset_samples$timepoint=="wk0",],dataset_samples_tech,by="sample_id"); rownames(s2a)<-s2a$sample_id
df2a<-merge(s2a[,c("disease",analysis_covs)],RNA,by=0)

genes<-grep("^ENSG", names(df2a), value=T)

covs.text<-paste(analysis_covs,collapse=" + ")

#### Register the analysis: create a new analysis in the database ####

# Payload fields: http://10.7.50.21:8000/docs#/analyses/create_analysis_analyses__post

new_analysis_payload <- list(
  analysis_goal = analysis_goal,
  analysis_scope = analysis_scope,
  analysis_type = analysis_type,
  analysis_subtype = analysis_subtype,
  analysis_method = analysis_method,
  analysis_design = analysis_design,
  omics_layer = assay1,
  disease_name = disease_of_interest,
  parameters = list("DataMolecular"=dataset_yaml$analysis_data,
                    "DataFilter"=paste(filt.genes.criterion,collapse=","),
                    "AnalysisCovariates"=paste(analysis_covs,collapse=",")),
  datasets = list(dataset_id),
  reference_sources = as.list(c(gene.anno_reference_id,gene.featureTable_id)),
  input_files = list()
)
new_analysis <- db_create_analysis(new_analysis_payload)

#### Execute the analysis ####

if (analysis_method=="LogisticRegression") {
  df2a$IMID <- factor(ifelse(df2a$disease == disease_of_interest, "1", "0"))
  
  df2a$disease <- factor(df2a$disease)
  df2a$sex <- factor(df2a$sex)
  
  if ("plateId" %in% colnames(df2a)) {
    df2a$plateId <- factor(df2a$plateId)
  } else if ("Sequencing_batch" %in% colnames(df2a)) {
    df2a$Sequencing_batch <- factor(df2a$Sequencing_batch)
  } else {
    stop("No batch variable found (plateId or Sequencing_batch)")
  }
  
  
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
}
if (analysis_method=="LimmaLinearModel") {
  
  meta <- df2a[, c("disease", analysis_covs)]
  meta$disease <- factor(meta$disease)
  meta$sex <- factor(meta$sex)
  
  if ("plateId" %in% colnames(meta)) {
    meta$plateId <- factor(meta$plateId)
  } else if ("Sequencing_batch" %in% colnames(meta)) {
    meta$Sequencing_batch <- factor(meta$Sequencing_batch)
  } else {
    stop("No batch variable found (plateId or Sequencing_batch)")
  }
  
  meta$disease <- relevel(meta$disease, ref = "CTRL")
  
  expr <- as.matrix(df2a[, !(colnames(df2a) %in% c("Row.names", "disease", analysis_covs))])
  expr <- t(expr)
  
  design <- model.matrix(as.formula(paste0("~ disease + ",covs.text)), data = meta)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  coef<-paste0("disease",disease_of_interest)
  res <- topTable(fit, coef = coef, number = Inf)
  
  res$gene <- rownames(res)
  res <- res[, c("gene", colnames(res)[!colnames(res) %in% "gene"])]
  res<-merge(res,corresp[,c("EnsemblID","gene_symbol")],by.x="gene",by.y="EnsemblID", all.x=T)
  res$TechID<-NA
  res<-res[,c("gene","TechID","gene_symbol","logFC","P.Value")]
  colnames(res)<-c("Ensembl_ID","TechID","Symbol","Effect","Pvalue")
  res<-res[order(res$Pvalue),]
  res<-res[!is.na(res$Pvalue),]
  
}
saveRDS(res,analysis_path)

#### Update Analysis Table ####

update_analysis_payload<-list(
  analysis_path = analysis_path,
  analysis_code = analysis_code
)
db_update_analysis(new_analysis$analysis_id, update_analysis_payload)

#### (OPTIONAL) Estimate computational resources  ####

# set.seed(1)
# sample_genes <- sample(genes, min(20, length(genes)))
# 
# # Benchmark Runtime
# mb <- microbenchmark(lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial"); coefficients(summary(model))[gene, "Estimate"]}),times=5)
# times_sec <- mb$time / 1e9
# 
# # CPU Timing
# cpu_time <- system.time({tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})})
# 
# # Memory Usage (MB aprox)
# gc_before <- gc()
# tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})
# gc_after <- gc()
# mem_used_mb <- sum(gc_after[,2] - gc_before[,2])*8 / 1024
# 
# # Peak memory
# peak <- peakRAM({tmp <- lapply(sample_genes, function(gene) {formula<-as.formula(paste0("IMID ~ ", gene, " + ", covs.text)); model<-glm(formula, data=df2a, family="binomial")})})
# peak_ram_mb <- peak$Peak_RAM_Used_MiB
# 
# # Disk usage (save results)
# disk_usage_mb <- file.size(paste0(output_dir,analysis_scope,"_",disease_of_interest,".rds")) / (1024^2)
# 
# # Derived metrics
# median_batch_time <- median(times_sec)
# genes_tested <- length(sample_genes)
# total_genes <- length(genes)
# 
# time_per_gene <- median_batch_time / genes_tested
# estimated_total_time <- time_per_gene * total_genes
# throughput <- genes_tested / median_batch_time
# 
# estimated_disk_total <- disk_usage_mb * (total_genes / genes_tested)
# 
# # Final table
# 
# resource_table <- data.frame(Metric = c("Median batch time (s)", "User CPU time (s)", "Elapsed time (s)", "Memory used (MB)", "Peak memory (MB)", "Disk usage sample (MB)", "Estimated disk usage total (MB)", "Genes per batch", "Total genes", "Time per gene (s)", "Throughput (genes/sec)", "Estimated total time (min)"),
#                              Value = c(median_batch_time, cpu_time["user.self"], cpu_time["elapsed"], mem_used_mb, peak_ram_mb, disk_usage_mb, estimated_disk_total, genes_tested, total_genes, time_per_gene, throughput, estimated_total_time/60))
# 

#### Inspection Analysis Results  ####

print(disease_of_interest)
dim(res)
table(df2a$disease)
head(res)

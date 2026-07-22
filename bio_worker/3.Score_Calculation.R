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
library(tidyr)

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

disease_of_interest <- "8 IMIDs"
tissue <- ""
assay1 <- ""
project1 <- ""
analysis_design <- ""
analysis_goal <- "TargetID"
analysis_subtype <- ""
analysis_scope <- "Scoring Individual Dimensions"

########  Analysis Workflow  ########

diseases<-grep("oJIA",unique(db_list_datasets()$disease),invert=T, value=T)

gene.anno_reference_id <- db_list_reference_sources(list(reference_type = "annotation", reference_name="Gene Annotation", is_default = TRUE))$reference_id
gene.anno_reference <- db_get_reference_source(gene.anno_reference_id)
corresp<-read.csv(paste0(volume_dir,gene.anno_reference$reference_path),sep="\t")

corresp_long <- corresp %>% separate_rows(aliases, sep = "\\|") %>% mutate(match_symbol = aliases) %>% bind_rows(corresp %>% mutate(match_symbol = gene_symbol)) %>% distinct(match_symbol, .keep_all = TRUE)

####  1.Causality  ####

dimensions<-c("eMR","pMR")

for (d in dimensions) {
  Scores<-list()
  for (imid in diseases) {
    if (d=="eMR") {
      x<-readRDS(paste0(output_dir,d,"_",imid,".rds"))[,c("ENSGene","RankingScore","Directionality")]
    } else {
      x<-readRDS(paste0(output_dir,d,"_",imid,".rds"))[,c("Symbol","RankingScore","Directionality")]
      x<-as.data.frame(x %>% left_join(corresp_long, by = c("Symbol" = "match_symbol")))[,c("EnsemblID","RankingScore","Directionality")]
    }
    x$Score<-(x$RankingScore - min(x$RankingScore)) / (max(x$RankingScore) - min(x$RankingScore))
    x2save<-x[,c(1,4,3)]; rownames(x2save)<-NULL; colnames(x2save)<-c("EnsemblID","Score","Direction")
    Scores[[imid]]<-x2save
  }
  saveRDS(Scores,paste0(output_dir,"Scores_",d,".rds"))
}

####  2.Genetic Regulation  ####

dimensions<-c("tPGS","pPGS")

for (d in dimensions) {
  Scores<-list()
  for (imid in diseases) {
      x<-readRDS(paste0(output_dir,d,"_",imid,".rds"))
      x<-x[!is.na(x$Pvalue),]
      x$Score<- (-log10(x$Pvalue) - min(-log10(x$Pvalue))) / (max(-log10(x$Pvalue)) - min(-log10(x$Pvalue)))
      x$Direction<-sign(x$Effect)
      x2save<-x[,c(1,2,6,7)]; rownames(x2save)<-NULL; colnames(x2save)<-c("EnsemblID","TechID","Score","Direction")
      Scores[[imid]]<-x2save
  }
  saveRDS(Scores,paste0(output_dir,"Scores_",d,".rds"))
}

####  3.Transcriptomics  ####

dimensions<-c("BulkBlood")

for (d in dimensions) {
  Scores<-list()
  for (imid in diseases) {
    x<-readRDS(paste0(output_dir,d,"_",imid,".rds"))
    x<-x[!is.na(x$Pvalue),]
    x$Score<- (-log10(x$Pvalue) - min(-log10(x$Pvalue))) / (max(-log10(x$Pvalue)) - min(-log10(x$Pvalue)))
    x$Direction<-sign(x$Effect)
    x2save<-x[,c(1,6,7)]; rownames(x2save)<-NULL; colnames(x2save)<-c("EnsemblID","Score","Direction")
    Scores[[imid]]<-x2save
  }
  saveRDS(Scores,paste0(output_dir,"Scores_",d,".rds"))
}

####  4.Proteomics  ####

dimensions<-c("PlasmaProteomics")

for (d in dimensions) {
  Scores<-list()
  for (imid in diseases) {
    x<-readRDS(paste0(output_dir,d,"_",imid,".rds"))
    x<-x[!is.na(x$Pvalue),]
    x$Score<- (-log10(x$Pvalue) - min(-log10(x$Pvalue))) / (max(-log10(x$Pvalue)) - min(-log10(x$Pvalue)))
    x$Direction<-sign(x$Effect)
    x2save<-x[,c(1,2,6,7)]; rownames(x2save)<-NULL; colnames(x2save)<-c("EnsemblID","TechID","Score","Direction")
    Scores[[imid]]<-x2save
  }
  saveRDS(Scores,paste0(output_dir,"Scores_",d,".rds"))
}




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
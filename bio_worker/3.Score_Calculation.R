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

# =========================================
# Accessing results for scoring calculation
# =========================================

disease_of_interest <- "AD"
tissue <- "blood"

assay1 <- "genetics"
project1 <- "ssad_expanded"

assay2 <- "transcriptomics"
project2 <- "ssad"

analysis <- "a1"

# ==================================
# mPGS Scoring 
# ==================================

tmp.folder<-paste0(paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/analyses/",analysis))

res<-readRDS(paste0(tmp.folder,"/PGS_",disease_of_interest,"_",assay2,"_",tissue,".rds"))
res$Score <- (-log10(res$Pvalue) - min(-log10(res$Pvalue))) / (max(-log10(res$Pvalue)) - min(-log10(res$Pvalue)))
res$Direction<-sign(res$Effect)
out<-res[,c("EnsemblID","TechID","gene_symbol","Effect","Pvalue","Score","Direction")]
saveRDS(out,paste0(tmp.folder,"/PGS_",disease_of_interest,"_",assay2,"_",tissue,"_Scores.rds"))


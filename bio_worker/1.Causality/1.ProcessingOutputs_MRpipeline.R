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

####  eMR  ####

system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/AD_BuduAggrey2023.rds ",output_dir,"eMR_AD.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/SS_Sakaue2021.rds ",output_dir,"eMR_SS.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/SLE_Bentham2015.rds ",output_dir,"eMR_SLE.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/PS_Dand2025.rds ",output_dir,"eMR_PS.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/PSA_Soomro2022.rds ",output_dir,"eMR_PSA.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/RA_Ishigaki2022.rds ",output_dir,"eMR_RA.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/UC_deLange2017.rds ",output_dir,"eMR_UC.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_eQTL/06_Scoring/CD_deLange2017.rds ",output_dir,"eMR_CD.rds"))

####  pMR  ####

system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/AD_BuduAggrey2023.rds ",output_dir,"pMR_AD.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/SS_Sakaue2021.rds ",output_dir,"pMR_SS.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/SLE_Bentham2015.rds ",output_dir,"pMR_SLE.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/PS_Dand2025.rds ",output_dir,"pMR_PS.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/PSA_Soomro2022.rds ",output_dir,"pMR_PSA.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/RA_Ishigaki2022.rds ",output_dir,"pMR_RA.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/UC_deLange2017.rds ",output_dir,"pMR_UC.rds"))
system(paste0("cp /media/bioinformatics/Projects/Internal_Projects/2026_IMX_TargetIdentification_vMR3/01_MR_pQTL/06_Scoring/CD_deLange2017.rds ",output_dir,"pMR_CD.rds"))





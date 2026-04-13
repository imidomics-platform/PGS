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

# ===========================================
# Accessing database for PGS analysis example 
# ===========================================

disease_of_interest <- "AD"
assay <- "genetics"
project <- "ssad_expanded"

# 1) Getting reference information for the specified disease and reference type
# ==============================================================================

# List reference sources if you want to see all available references in the system and you can filter them.
refs <- list_reference_sources()

# List reference sources with filters
GWAS_ref_summary <- list_reference_sources(list(reference_type = "reference GWAS summary statistics", disease = disease_of_interest), limit = 50)
GWAS_PGS_pvals <- list_reference_sources(list(reference_type = "PGS variant p-values", disease = disease_of_interest), limit = 50)
GWAS_PGS_weights <- list_reference_sources(list(reference_type = "PGS variant weights", disease = disease_of_interest), limit = 50)
GWAS_ref_info <- list_reference_sources(list(reference_type = "reference GWAS and associated-PGS information", disease = disease_of_interest), limit = 50)
# Get PGS thresholds reference source (not disease-specific)
PGS_thresholds <- list_reference_sources(list(reference_name = "Significance Thresholds for PGS calculation"), limit = 50)

# If there are multiple reference sources for the same disease and reference type, you can check for the is_default flag to identify the default reference source for that disease and reference type
# or filter based on reference_name manually to keep only one reference source.
GWAS_ref_summary_default <- list_reference_sources(list(reference_type = "reference GWAS summary statistics", disease = disease_of_interest, is_default = TRUE), limit = 50)
GWAS_PGS_pvals_default <- list_reference_sources(list(reference_type = "PGS variant p-values", disease = disease_of_interest, is_default = TRUE), limit = 50)
GWAS_PGS_weights_default <- list_reference_sources(list(reference_type = "PGS variant weights", disease = disease_of_interest, is_default = TRUE), limit = 50)
GWAS_ref_info_default <- list_reference_sources(list(reference_type = "reference GWAS and associated-PGS information", disease =  disease_of_interest, is_default = TRUE), limit = 50)

# Get the file paths for the reference sources (assuming the first reference source is the one we want to use if there are multiple)
GWAS_ref_summary_default_path <- paste0(volume_dir, "/", GWAS_ref_summary_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_PGS_pvals_default_path <- paste0(volume_dir, "/", GWAS_PGS_pvals_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_PGS_weights_default_path <- paste0(volume_dir, "/", GWAS_PGS_weights_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_ref_info_default_path <- paste0(volume_dir, "/", GWAS_ref_info_default$reference_path[1]) # assuming the first reference source is the one we want to use
PGS_thresholds_path <- paste0(volume_dir, "/", PGS_thresholds$reference_path[1]) # assuming the first reference source is the one we want to use

GWAS_ref_info <- read.delim(GWAS_ref_info_default_path, header = TRUE, stringsAsFactors = FALSE)
GWAS_ref_id <- GWAS_ref_info$GWAS_id
GWAS_reference_genome <- GWAS_ref_info$REF.ANNO
GWAS_pgs.cat.id <- GWAS_ref_info$PGS.CAT.ID

# 2) Getting dataset information for the specified disease, assay type and project
# Assuming the first dataset is the one we want to use
# ================================================================================

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project, assay = assay), limit = 50)
dataset_path    <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")

# read dataset YAML file to get the path to the processed PGS file for the specified disease
dataset_yaml <- yaml::read_yaml(dataset_yaml_path)
dataset_gwas_hg19_path <- paste0(dataset_path, "/", dataset_yaml$gwas_hg19_data)
dataset_gwas_hg38_path <- paste0(dataset_path, "/", dataset_yaml$gwas_hg38_data)
dataset_pc_table_path <- paste0(dataset_path, "/", dataset_yaml$pc_data)
dataset_pgs_cat <- paste0(dataset_path, "/", dataset_yaml$pgs.cat_data)
dataset_pgs_dict <- paste0(dataset_path, "/", dataset_yaml$pgs.cat_dict)
dataset_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_pcs<-read.delim(paste0(dataset_path, "/", dataset_yaml$pc_data),header=T)

print(paste0(disease_of_interest," | ",assay," | ",project," | ",dataset_info$dataset_name, " | ", GWAS_ref_id))

# =================================================
# PGS Calculation 
# =================================================

tmp.folder<-paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/")

#---------------- Functions -----------------#

read_and_rename <- function(file) {
  df <- read.table(file, header=T, stringsAsFactors=F)[,c(2,6)]
  colnames(df) <- c("sample_id",gsub("PGS_R2_0.","LD0",gsub(".profile","",tail(strsplit(file, "/")[[1]], 1))))
  return(df)
}

#---------------- Harmonization -------------#

filt.data<-readRDS(GWAS_ref_summary_default_path)
if(GWAS_reference_genome=="hg19") {dataset_gwas<-dataset_gwas_hg19_path} else {dataset_gwas<-dataset_gwas_hg38_path}
bim.imx<-read.table(paste0(dataset_gwas,".bim")); colnames(bim.imx)<-c("CHR","SNP","CM","POS","A1","A2")

dfm<-merge(bim.imx,filt.data[,c("SNP","REF.ALLELE","OTH.ALLELE","STDERR","EFFECT","PVALUE")], by="SNP")
dfm<-dfm[(nchar(dfm$REF.ALLELE)<2 & nchar(dfm$OTH.ALLELE)<2),]
dfm<-dfm[(toupper(dfm$REF.ALLELE)==dfm$A1 | toupper(dfm$REF.ALLELE)==dfm$A2),]

#---------------- Thresholding Approach ----#

tmp.gwas.ref.folder<-paste0(tmp.folder,project,"_",disease_of_interest,"_",GWAS_ref_info$GWAS_id,"_",GWAS_reference_genome)
dir.create(tmp.gwas.ref.folder)
write.table(data.frame("SNPS.REF"=nrow(filt.data), "SNPS.IMX"=nrow(bim.imx), "SNPS.COMMON.QC"=nrow(dfm), "SNPS.COMMON.QC.PERC"=100*(nrow(dfm)/nrow(bim.imx))), paste0(tmp.gwas.ref.folder,"/CommonSNPs.txt"), quote=F, row.names=F, sep="\t")

ths.ld<-c(0.2,0.5,0.9)
ths.kb<-250

for (th in ths.ld) {
  system(paste0("plink_v1.9 --bfile ",dataset_gwas, " --clump ",GWAS_PGS_pvals_default_path, " --clump-p1 1 --clump-p2 1 --clump-kb ",ths.kb," --clump-r2 ",th, " --out ",tmp.gwas.ref.folder,"/Clumped_R2_",th))
  system(paste0("awk '{print $3}' ",tmp.gwas.ref.folder,"/Clumped_R2_",th,".clumped > ",tmp.gwas.ref.folder,"/Clumped_SNPs_R2_",th,".txt"))
  system(paste0("plink_v1.9 --bfile ",dataset_gwas," --score ",GWAS_PGS_weights_default_path," --q-score-range ",PGS_thresholds_path," ",GWAS_PGS_pvals_default_path," --extract ",tmp.gwas.ref.folder,"/Clumped_SNPs_R2_",th,".txt --out ",tmp.gwas.ref.folder,"/PGS_R2_",th))
  system(paste0("rm ",tmp.gwas.ref.folder,"/*.log | rm ",tmp.gwas.ref.folder,"/*.clumped | rm ",tmp.gwas.ref.folder,"/*.nosex | rm ",tmp.gwas.ref.folder,"/*.nopred | rm ",tmp.gwas.ref.folder,"/Clumped_*"))
}

df.list<-map(list.files(tmp.gwas.ref.folder, pattern="*.profile", full.names=T), read_and_rename)

pgs.ext<-read.delim(dataset_pgs_cat, header = TRUE)
pgs.dict<-read.delim(dataset_pgs_dict, header = TRUE)
pgs.data <- merge(purrr:::reduce(df.list, full_join, by="sample_id"),pgs.ext, by="sample_id")

#---------------- K-Fold Cross-Validation ----#

k<-5

# Getting optimal combination of thresholds

colnames(dataset_pcs)[2:21]<-paste0("Gen.",colnames(dataset_pcs)[2:21])
data0<-merge(dataset_samples, dataset_pcs, by="sample_id")
data<-merge(data0,pgs.data, by="sample_id")

set.seed(123)
folds<-createFolds(data$disease, k=k, list=T)
cols2test<-grep("LD",colnames(data),value=T)[!(grep("LD",colnames(data),value=T) %in% c("LD02","LD05","LD09"))]

results<-data.frame(Threshold=cols2test, Nagelkerke.R2.All=NA, Nagelkerke.R2.Gen=NA, Nagelkerke.R2.Epi=NA)
for (col in cols2test) {
  R2.GEN<-list(); R2.EPI<-list(); R2.ALL<-list();
  for (i in 1:k) {
    train.data <- data[unlist(folds[-i]), ]
    test.data <- data[unlist(folds[i]), ]
    train.data$IMID.Num<-ifelse(train.data$disease==disease_of_interest, 1, 0)
    test.data$IMID.Num<-ifelse(test.data$disease==disease_of_interest, 1, 0)
    train.data$Th2a<-train.data[[col]]
    test.data$Th2a<-test.data[[col]]
    
    model.all<-glm(IMID.Num~Th2a+Gen.PC1+Gen.PC2+sex+age, data=train.data, family=binomial)
    model.epi<-glm(IMID.Num~sex+age, data=train.data, family=binomial)
    model.gen<-glm(IMID.Num~Th2a+Gen.PC1+Gen.PC2, data=train.data, family=binomial)
    
    R2.GEN[[i]]<-NagelkerkeR2(model.gen)$R2
    R2.EPI[[i]]<-NagelkerkeR2(model.epi)$R2
    R2.ALL[[i]]<-NagelkerkeR2(model.all)$R2
  }
  results[results$Threshold == col, "Nagelkerke.R2.All"]<-mean(unlist(R2.ALL))
  results[results$Threshold == col, "Nagelkerke.R2.Gen"]<-mean(unlist(R2.GEN))
  results[results$Threshold == col, "Nagelkerke.R2.Epi"]<-mean(unlist(R2.EPI))
}

count<-nrow(results)
pgs.cat.id<-unlist(strsplit(GWAS_pgs.cat.id,"_"))[unlist(strsplit(GWAS_pgs.cat.id,"_")) %in% colnames(pgs.ext)]
empty_rows <- data.frame(matrix(NA, nrow=length(pgs.cat.id), ncol = ncol(results)))
colnames(empty_rows) <- colnames(results)
results <- rbind(results, empty_rows)
for (j in pgs.cat.id) {
  count=count+1
  data$PGS<-data[[j]]
  data$IMID.Num<-ifelse(data$disease==disease_of_interest, 1, 0)
  model.all<-glm(IMID.Num~PGS+Gen.PC1+Gen.PC2+sex+age, data=data, family=binomial)
  model.epi<-glm(IMID.Num~sex+age, data=data, family=binomial)
  model.gen<-glm(IMID.Num~PGS+Gen.PC1+Gen.PC2, data=data, family=binomial)
  results$Threshold[count]<-j
  results$Nagelkerke.R2.All[count]<-NagelkerkeR2(model.all)$R2
  results$Nagelkerke.R2.Gen[count]<-NagelkerkeR2(model.gen)$R2
  results$Nagelkerke.R2.Epi[count]<-NagelkerkeR2(model.epi)$R2
}

# Applying optimal thresholds to full dataset

opt.ths<-results[which.max(results$Nagelkerke.R2.All),"Threshold"]

data$Optimal<-data[[opt.ths]]
data$IMID.Num<-ifelse(data$disease==disease_of_interest, 1, 0)
model.final.all<-glm(IMID.Num~Optimal+Gen.PC1+Gen.PC2+sex+age, data=data, family=binomial)
model.final.gen<-glm(IMID.Num~Optimal+Gen.PC1+Gen.PC2, data=data, family=binomial)
model.final.epi<-glm(IMID.Num~age+sex, data=data, family=binomial)

results$N.SNPs<-NA
results$Nagelkerke.R2.All.FullDataset<-NA
results$Nagelkerke.R2.Gen.FullDataset<-NA
results$Nagelkerke.R2.Epi.FullDataset<-NA
results$Nagelkerke.Pval.All.FullDataset<-NA
results$Nagelkerke.Pval.Gen.FullDataset<-NA
results[results$Threshold==opt.ths,"Nagelkerke.R2.All.FullDataset"]<-NagelkerkeR2(model.final.all)$R2
results[results$Threshold==opt.ths,"Nagelkerke.R2.Gen.FullDataset"]<-NagelkerkeR2(model.final.gen)$R2
results[results$Threshold==opt.ths,"Nagelkerke.R2.Epi.FullDataset"]<-NagelkerkeR2(model.final.epi)$R2
results[results$Threshold==opt.ths,"Nagelkerke.Pval.All.FullDataset"]<-coefficients(summary(model.final.all))["Optimal","Pr(>|z|)"]
results[results$Threshold==opt.ths,"Nagelkerke.Pval.Gen.FullDataset"]<-coefficients(summary(model.final.gen))["Optimal","Pr(>|z|)"]

for (i in 1:nrow(results)) {
  th<-results$Threshold[i]
  if (startsWith(th,"LD")) {
    fname<-paste0(tmp.gwas.ref.folder,'/PGS_R2_', substr(gsub("LD","",unlist(strsplit(th,"\\."))[1]), 1, 1),
                  ".",
                  substr(gsub("LD","",unlist(strsplit(th,"\\."))[1]), 2, nchar(gsub("LD","",unlist(strsplit(th,"\\."))[1]))),
                  ".",
                  unlist(strsplit(th,"\\."))[2],".profile")
    results$N.SNPs[i]<-round(mean(fread(fname)$CNT)/2)
  } else {
    x<-pgs.dict[pgs.dict$score==th,"variants_used"]
    if (length(x) == 0) x<-NA
    results$N.SNPs[i]<-x
  }
}

# Summary optimal threshold

p.opt<-ggboxplot(data, x="IMID.Num", y="Optimal", color="IMID.Num", add="jitter", fill="IMID.Num", add.params=list(size=1,jitter=0.33)) +
  scale_color_manual(values=rep("black",2)) +
  scale_fill_manual(values=c("0"="#0F7D0F", "1"="#AF0C0C")) +
  scale_x_discrete(labels = c("0" = "CTRL", "1" = disease_of_interest)) +
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=11, face="bold")) +
  labs(x="",
       y=paste0("Optimal ",disease_of_interest,"-PGS"),
       title=paste0(opt.ths," (n=",results[results$Threshold==opt.ths,"N.SNPs"]," SNPs)\n",
                    "Nagelkerke.R2=",round(results[results$Threshold==opt.ths,"Nagelkerke.R2.All.FullDataset"],2),
                    ";  Pval=",format(results[results$Threshold==opt.ths,"Nagelkerke.Pval.All.FullDataset"], scientific=T, digits=3)))

dir.create(paste0(tmp.gwas.ref.folder,"/PGSprofiles"))
system(paste0("mv ",tmp.gwas.ref.folder,"/* ",tmp.gwas.ref.folder,"/PGSprofiles/"))

saveRDS(results, paste0(tmp.gwas.ref.folder,"/PGSprofiles/ThresholdingSummary.rds"))

jpeg(paste0(tmp.gwas.ref.folder,"/PGSprofiles/Optimal.jpeg"), res=300, height=1500, width=1800, type="cairo")
plot(p.opt)
dev.off()

#system(paste0("tar -czvf ",tmp.gwas.ref.folder,"/PGS_",disease_of_interest,"_Profiles.tar.gz ",tmp.gwas.ref.folder,"/PGSprofiles"))
#system(paste0("rm -r ",tmp.gwas.ref.folder,"/PGSprofiles/"))

# =======================================
# Input Files for Downstream Analyses
# =======================================

# tmp/prject_disease_ref_ref.anno/PGSprofiles/ThresholdingSummary.rds
# tmp/prject_disease_ref_ref.anno/PGSprofiles/Optimal.jpeg
# tmp/prject_disease_ref_ref.anno/PGSprofiles/PGS_*.profile optimal indicated in ThresholdingSummary.rds

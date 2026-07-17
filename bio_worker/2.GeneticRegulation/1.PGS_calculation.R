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
output_dir <- paste0(tmp_dir,"PGS/")
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

disease_of_interest <- "UC"
assay1 <- "genetics"
tissue <- "blood"
project1 <- "Cross-sectional"

analysis_goal<-"PGS"
analysis_scope<-""
analysis_type<-""
analysis_subtype<-paste0(assay1,"_",tissue)
analysis_method<-"Thresholding Approach and K-Fold Cross-Validation"
analysis_design<-""

dataset_covariates<-"Gen.PC1, Gen.PC2, sex, age"

#### Getting reference information for the specified input ####

GWAS_ref_summary <- list_reference_sources(list(reference_type = "reference GWAS summary statistics", disease = disease_of_interest))
GWAS_PGS_pvals <- list_reference_sources(list(reference_type = "PGS variant p-values", disease = disease_of_interest))
GWAS_PGS_weights <- list_reference_sources(list(reference_type = "PGS variant weights", disease = disease_of_interest))
GWAS_ref_info <- list_reference_sources(list(reference_type = "reference GWAS and associated-PGS information", disease = disease_of_interest))
PGS_thresholds <- list_reference_sources(list(reference_name = "Significance Thresholds for PGS calculation"))

GWAS_ref_summary_default <- list_reference_sources(list(reference_type = "reference GWAS summary statistics", disease = disease_of_interest, is_default = TRUE))
GWAS_PGS_pvals_default <- list_reference_sources(list(reference_type = "PGS variant p-values", disease = disease_of_interest, is_default = TRUE))
GWAS_PGS_weights_default <- list_reference_sources(list(reference_type = "PGS variant weights", disease = disease_of_interest, is_default = TRUE))
GWAS_ref_info_default <- list_reference_sources(list(reference_type = "reference GWAS and associated-PGS information", disease =  disease_of_interest, is_default = TRUE))

GWAS_ref_summary_default_path <- paste0(volume_dir, "/", GWAS_ref_summary_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_PGS_pvals_default_path <- paste0(volume_dir, "/", GWAS_PGS_pvals_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_PGS_weights_default_path <- paste0(volume_dir, "/", GWAS_PGS_weights_default$reference_path[1]) # assuming the first reference source is the one we want to use
GWAS_ref_info_default_path <- paste0(volume_dir, "/", GWAS_ref_info_default$reference_path[1]) # assuming the first reference source is the one we want to use
PGS_thresholds_path <- paste0(volume_dir, "/", PGS_thresholds$reference_path[1]) # assuming the first reference source is the one we want to use

GWAS_ref_info <- read.delim(GWAS_ref_info_default_path, header = TRUE, stringsAsFactors = FALSE)
GWAS_ref_id <- GWAS_ref_info$GWAS_id
GWAS_reference_genome <- GWAS_ref_info$REF.ANNO
GWAS_pgs.cat.id <- GWAS_ref_info$PGS.CAT.ID

#### Getting dataset information for the specified input ####

dataset_info <- list_datasets(list(disease = disease_of_interest, dataset_project=project1, assay = assay1))
dataset_path    <- paste0(volume_dir, "/", dataset_info$dataset_path[1]) 
dataset_version <- dataset_info$dataset_version[1]
dataset_yaml_path <- paste0(dataset_path, "/", dataset_version, ".yaml")

dataset_yaml <- yaml::read_yaml(dataset_yaml_path)
dataset_gwas_hg19_path <- paste0(dataset_path, "/", dataset_yaml$gwas_hg19_data)
dataset_gwas_hg38_path <- paste0(dataset_path, "/", dataset_yaml$gwas_hg38_data)
dataset_pc_table_path <- paste0(dataset_path, "/", dataset_yaml$pc_data)
dataset_pgs_cat <- paste0(dataset_path, "/", dataset_yaml$pgs.cat_data)
dataset_pgs_dict <- paste0(dataset_path, "/", dataset_yaml$pgs.cat_dict)
dataset_samples<-read.delim(paste0(dataset_path, "/samples/samples.tsv"),header=T)
dataset_pcs<-read.delim(paste0(dataset_path, "/", dataset_yaml$pc_data),header=T)

print(paste0(disease_of_interest," | ",assay1," | ",project1," | ",dataset_info$dataset_name, " | ", GWAS_ref_id))

####  PGS Calculation  ####

tmp.gwas.ref.folder<-paste0(output_dir,disease_of_interest,"_",GWAS_ref_info$GWAS_id,"_",GWAS_reference_genome)
dir.create(tmp.gwas.ref.folder)

# Functions

read_and_rename <- function(file) {
  df <- read.table(file, header=T, stringsAsFactors=F)[,c(2,6)]
  colnames(df) <- c("sample_id",gsub("PGS_R2_0.","LD0",gsub(".profile","",tail(strsplit(file, "/")[[1]], 1))))
  return(df)
}

# 1. Harmonization

filt.data<-readRDS(GWAS_ref_summary_default_path)
if(GWAS_reference_genome=="hg19") {dataset_gwas<-dataset_gwas_hg19_path} else {dataset_gwas<-dataset_gwas_hg38_path}
bim.imx<-read.table(paste0(dataset_gwas,".bim")); colnames(bim.imx)<-c("CHR","SNP","CM","POS","A1","A2")

dfm<-merge(bim.imx,filt.data[,c("SNP","REF.ALLELE","OTH.ALLELE","STDERR","EFFECT","PVALUE")], by="SNP")
dfm<-dfm[(nchar(dfm$REF.ALLELE)<2 & nchar(dfm$OTH.ALLELE)<2),]
dfm<-dfm[(toupper(dfm$REF.ALLELE)==dfm$A1 | toupper(dfm$REF.ALLELE)==dfm$A2),]

# 2. Thresholding Approach

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

# 3. K-Fold Cross-Validation to get optimal combination of thresholds

k<-5

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

# 4. Applying optimal thresholds to full dataset

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

data.optimal<-results[!is.na(results$Nagelkerke.R2.All.FullDataset),]

####  Register the analysis: create a new analysis in the database  ####

# Payload fields: http://10.7.50.21:8000/docs#/analyses/create_analysis_analyses__post

new_analysis_payload <- list(
  analysis_goal = analysis_goal,
  analysis_scope = analysis_scope,
  analysis_type = analysis_type,
  analysis_subtype = paste0(assay1,"_",tissue),
  analysis_method = analysis_method,
  analysis_design = analysis_design,
  parameters = list("Data.Molecular"=dataset_gwas,"Data.Covariates"=dataset_covariates,"Data.Optimal"=data.optimal),
  omics_layer = assay1,
  disease_name = disease_of_interest,
  datasets = list(dataset_id),
  reference_sources = as.list(c(GWAS_ref_summary_default$reference_id,GWAS_PGS_pvals_default$reference_id,GWAS_PGS_weights_default$reference_id,GWAS_ref_info_default$reference_id,PGS_thresholds$reference_id)),
  input_files = list()
)

new_analysis <- db_create_analysis(new_analysis_payload)

update_analysis_payload <- list(
  analysis_path = paste0(tmp.gwas.ref.folder,"/PGSprofiles/"),
  analysis_code = gsub("/home/aaterido/Escritorio/GitHub_Repositories/","/",paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/",basename(rstudioapi::getSourceEditorContext()$path)))
)

db_update_analysis(new_analysis$analysis_id, update_analysis_payload)

db_list_analyses()


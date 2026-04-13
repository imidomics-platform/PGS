library(ggpubr)
library(dplyr)
library(fmsb)
library(limma)
library(edgeR)
library(Matrix)
library(caret)
library(grid)
library(gridExtra)
library(TwoSampleMR)
library(viridis)
library(ggstatsplot)
library(magick)
library(ggplotify)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggrepel)
library(ggplot2)
library(fgsea)
library(bc3net)
library(khroma)
library(igraph)
library(tidygraph)
library(ggraph)

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

#----------------- Functions -----------------#

evaluate_hypergeometric <- function(scores.df, eval_geneset, fillZeros=F, percentile=F, npoints=100, by = 10){
  require(ggplot2)
  require(tidyr)
  
  # First column of scores.df must be Symbol, following columns are scores to evaluate
  if(names(scores.df)[1]!="Symbol") stop("First column of scores.df must be named `Symbol`")
  
  # Percentile dataframe
  tab <- data.frame(Pos=1:npoints)
  if(percentile == T){
    tab$Percentile <- tab$Pos/nrow(tab)
  } else{
    tab$Top <- seq(from=10, by=by, length.out=npoints)
  }
  
  # Compute enrichment
  # Compute for every column other than the first
  for(i in 1:nrow(tab)){
    for(a in 2:ncol(scores.df)){
      df <- scores.df[,c(1,a)]
      df <- df[order(df[,2],decreasing = T, na.last = T),]
      if(fillZeros==T) df[is.na(df[,2]), 2] <- 0
      df <- df[!is.na(df[,2]),]
      
      # Random reorder if multiple 0s
      if(sum(df[,2]==0)>1){
        firstZero <- which(df[,2]==0)[1]
        df[firstZero:nrow(df),] <- df[sample(x = firstZero:nrow(df), size = length(firstZero:nrow(df)), replace = F),]
      }
      
      # Candidate set
      if(percentile == T){
        candidate <- df$Symbol[1:(tab$Percentile[i]*nrow(df))]
      } else{
        candidate <- df$Symbol[1:tab$Top[i]]
      }
      
      # Reference set
      geneset <- eval_geneset[eval_geneset %in% df$Symbol]
      
      # Compute metric 
      res <- enrichment(genes = candidate, reference = df$Symbol, genesets = list(PositiveControls=geneset))
      tab[[colnames(scores.df)[a]]][i] <- -log10(res$pval)
    }
  }
  
  # Dataframe for plotting
  df2 <- gather(data = tab[1:(nrow(tab)/2),], key = category, value = Hypergeometric_logP, names(tab)[-c(1,2)])
  
  return(df2)
}

tol.palette <- list()
tol.palette[[1]]=c("#4477AA")
tol.palette[[2]]=c("#4477AA", "#CC6677")
tol.palette[[3]]=c("#4477AA", "#DDCC77", "#CC6677")
tol.palette[[4]]=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol.palette[[5]]=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol.palette[[6]]=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol.palette[[7]]=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol.palette[[8]]=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol.palette[[9]]=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol.palette[[10]]=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol.palette[[11]]=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol.palette[[12]]=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
tol.palette[[14]]=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol.palette[[15]]=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol.palette[[18]]=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
tol.palette[[21]]= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

tolsPalette <- function(x){
  
  if(x <= 12){
    palette <- tol.palette[[x]]
  } else if(x > 12 & x <= 23){
    iridescent <- colour("iridescent")
    palette <- iridescent(x)
  } else if(x > 23){
    smoothRainbow <- colour("smooth rainbow")
    palette <- smoothRainbow(x)
  }
  
  return(palette)
}

count_na <- function(x) {sum(is.na(x))}

params<-list("npoints"=100,"by"=7)

# =========================================
# Accessing reference data
# =========================================

ensembl.tss<-readRDS(paste0(volume_dir, "/reference/gene_reference_data/All_ENS_SYMB_TSS_v37.rds"))
newpos<-readRDS(paste0(volume_dir, "/reference/therapeutic_data_repository/gset_knownDrugs.rds")); newpos$IMID<-unique(unlist(newpos))
newneg<-readRDS(paste0(volume_dir, "/reference/therapeutic_data_repository/curated_KDT_v1.rds")); newneg$IMID<-unique(unlist(newneg))
nei<-readRDS(paste0(volume_dir, "/reference/gene_reference_data/Omnipath_Omnipath_Interactors_filteredAug22.rds"))
newnei <- lapply(names(newpos), function(x) {unique(unlist(nei[newpos[[x]]]))});  names(newnei)<-names(newpos)

# =========================================
# Accessing results for scoring calculation
# =========================================

disease_of_interest <- "SS"
tissue <- "blood"

assay1 <- "genetics"
project1 <- "ssad_expanded"

assay2 <- "transcriptomics"
project2 <- "ssad"

analysis <- "a2"

# ==================================
# mPGS Scoring 
# ==================================

tmp.folder<-paste0(paste0("/media/bioinformatics/imidomics_platform/volume_migration/tmp/analyses/",analysis))
scores.df<-readRDS(paste0(tmp.folder,"/PGS_",disease_of_interest,"_",assay2,"_",tissue,"_Scores.rds"))[,c("gene_symbol","Score")]
colnames(scores.df)<-c("Symbol","Score")

tab.dis<-evaluate_hypergeometric(scores.df = scores.df, eval_geneset = unique(c(newpos[[disease_of_interest]], newnei[[disease_of_interest]])), percentile=F, fillZeros=T, npoints = params$npoints, by = params$by)
tab.all<-evaluate_hypergeometric(scores.df = scores.df, eval_geneset = unique(c(newpos[["IMID"]], newnei[["IMID"]])), percentile = F, fillZeros=T, npoints = params$npoints, by = params$by)
tab.dis$KDTs<-"Disease-Specific Targets"
tab.all$KDTs<-"All IMID Targets"
tab<-rbind(tab.dis,tab.all)

ncolors<-nlevels(factor(tab$category))
pal.colors<-tolsPalette(ncolors)

p.dim<-ggplot(tab, aes(x=Top, y=Hypergeometric_logP, linetype=KDTs, color=category)) + 
  geom_line(aes(size=0.1), lwd=1, position=position_dodge(width=5)) + scale_color_manual(values=pal.colors) + theme_classic() +
  ylim(0, max(tab$Hypergeometric_logP)+1) +
  geom_line(aes(y = -log10(0.05)), color = "red", linetype = "dashed") +
  guides(colour="none", size="none") + labs(x="\nTop Ranking Positions\n", y="\nEnrichment in Known Drug Targets (-log10(Pval))\n") +
  theme(legend.title = element_blank(), legend.position=c(0.75,0.18))

hits.in.top<-unique(c(newpos[[disease_of_interest]], newnei[[disease_of_interest]]))[unique(c(newpos[[disease_of_interest]], newnei[[disease_of_interest]])) %in% scores.df[1:300,c("Symbol")]]
out.v0<-list()
count=0
for(j in 1:length(hits.in.top)) {
  query<-hits.in.top[j]
  for (target in newpos[[disease_of_interest]]) {
    if (query %in% nei[newpos[[disease_of_interest]]][[target]]) {count=count+1; out.v0[[count]]<-c(query,target)}
  }
}
out<-as.data.frame(do.call(rbind,out.v0))
colnames(out)<-c("Interactor","Known.Target")

g <- graph_from_data_frame(out, directed = TRUE)
tg <- as_tbl_graph(g)
V(tg)$color<-"#FCDB70"
V(tg)$color[V(tg)$name %in% unique(out$Known.Target)]<-"#B31515" 

p.net<-ggraph(tg, layout = "fr") + labs(title="") + theme(plot.margin = margin(20, 2, 2, 20)) +
  geom_edge_link() + geom_node_point(aes(color = color), size=6) +
  geom_node_text(aes(label=name), repel=T, size = 3) +
  scale_color_identity() +
  theme_void()

ggsave(filename = paste0(tmp.folder, "/PGS_", disease_of_interest, "_", assay2, "_", tissue, "_Evaluation.jpeg"),
  plot = grid.arrange(p.dim, p.net, ncol = 2, widths = c(1.4, 1.3), top = textGrob(paste0(disease_of_interest, "\n"), gp = gpar(fontsize=40, fontface = "bold"), hjust = 0.5)),
  width = 39, height = 20, units = "cm", dpi = 300)



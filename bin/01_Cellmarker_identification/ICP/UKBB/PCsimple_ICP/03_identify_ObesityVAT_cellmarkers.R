library(data.table)
library(xgboost)
library(stringr)
library(GeneralisedCovarianceMeasure)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/GCM.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/Github/get_cellmarkers.R")

scRNA <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/ObesityVAT/ObesityVAT_scRNA_exp.csv") # dim(scRNA): 222077  21487
scRNA <- scRNA[,-1] # the 1st column is index
meta <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/06_obese_adipose/h5ad/adiposeVAT_training_meta.csv")
# Select genes that overlap with UKBB datasets
# scRNA <- scRNA[,-1]
UKBB_genes_dict <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/06_ukbb/ukbb_tissue_genes/ukbb_genes_adiposeVAT.csv")
scRNA_genes <- colnames(scRNA)
UKBB_genes <- as.vector(unlist(UKBB_genes_dict$gene))
common_genes <- intersect(scRNA_genes,UKBB_genes)
scRNA <- subset(scRNA,select=common_genes)


celltypes_all <- as.vector(unlist(meta$CellType))
celltypes <- unique(celltypes_all)
bmi <- as.vector(unlist(meta$bmi))
disease_state <- rep(0,length(bmi))
case_index <- which(bmi>30)
ctrl_index <- which(bmi<=30)
disease_state[case_index] <- 1
# envir <- disease_state
# envir <- as.vector(unlist(meta$disease_state))
# case_index <- which(envir=="T2D")
# ctrl_index <- which(envir=="Control")
# envir_binary <- rep(0,length(envir))
# envir_binary[case_index] <- 1
envir_binary <- disease_state
# use xgboost to build a prediction model and extract variable importance for T2D
data_fir <- scRNA[case_index,]
data_fir <- as.matrix(data_fir)
#data_fir <- as.matrix(scRNA_ctrl[,(1:(feat_len-1))])
class(data_fir) <- "numeric"
label_fir <- celltypes_all[case_index]
#label_fir_binary <- rep(0,length(label_fir))

# use xgboost to build a prediction model and extract variable importance for controls
data_sec <- scRNA[ctrl_index,]
data_sec <- as.matrix(data_sec)
#data_fir <- as.matrix(scRNA_ctrl[,(1:(feat_len-1))])
class(data_sec) <- "numeric"
label_sec <- celltypes_all[ctrl_index]

label_binary <- rep(0,length(celltypes_all))
sample_num <- length(label_binary)

res_list <- get_cellmarkers(scRNA,celltypes_all, envir_binary,data_fir,label_fir,data_sec,label_sec,celltypes,sample_num)
out_filepath <- "/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/UKBB/PCsimple_ICP/"
save(res_list,file=paste0(out_filepath,"ObesityVAT_cellmarker_res.RData"))
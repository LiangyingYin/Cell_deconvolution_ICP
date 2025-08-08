source("/mnt/home/yinly/projects/Cell_deconvolution/bin/Github/get_cellmarkers.R")

library(data.table)
library(xgboost)
# library(stringr)
library(GeneralisedCovarianceMeasure)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(pcalg)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/GCM.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/PCSelect_Parallel_Update.R")

# Read single cell dataset
scRNA_h5ad <- readH5AD('/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/09_pbmc_data/01_preprocess_h5ad/pbmc_6810kABC.h5ad',  reader = c("R"))

counts(scRNA_h5ad) <- assay(scRNA_h5ad, "X")
scRNA_h5ad <- scuttle::logNormCounts(scRNA_h5ad)
scRNA <- as.matrix(assay(scRNA_h5ad, "logcounts"))
celltypes_all_origi <- as.vector(unlist(scRNA_h5ad$cell_type))
valid_index <- which(celltypes_all_origi!="Megakaryocytes")
scRNA <- scRNA[,valid_index]
celltypes_all <- celltypes_all_origi[valid_index]
scRNA <- t(scRNA)

# Select genes that overlap with bulk datasets
bulk_genes_dict <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/04_real_application/Newman/Newman_exp.txt")
scRNA_genes <- colnames(scRNA)
bulk_genes <- colnames(bulk_genes_dict)[-1]
common_genes <- intersect(scRNA_genes,bulk_genes)
scRNA <- subset(scRNA,select=common_genes)

celltypes <- unique(celltypes_all)

# Get Batch info
batch_origi <- as.vector(unlist(scRNA_h5ad$IID))
batch <- batch_origi[valid_index] # remove those for Megakaryocytes
batch_8k_index <- which(batch=="pbmc8k")
batch_10k_index <- which(batch=="pbmc10k")
batch2_index <- sort(c(batch_8k_index,batch_10k_index))
batch_all_index <- seq(1,length(batch))
batch1_index <- setdiff(batch_all_index,batch2_index)
envir_binary <- rep(0,length(batch))
envir_binary[batch2_index] <- 1

# use xgboost to build a prediction model and extract important cell markers for each cell type
# build environment-specific prediction model using xgboost and extract the cellmarkers for each celltype
data_fir <- scRNA[batch1_index,]
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- celltypes_all[batch1_index]

data_sec <- scRNA[batch2_index,]
data_sec <- as.matrix(data_sec)
class(data_sec) <- "numeric"
label_sec <- celltypes_all[batch2_index]

sample_num <- length(celltypes_all)
res_list <- get_cellmarkers(data_fir,label_fir,data_sec,label_sec,celltypes,sample_num)
library(data.table)
library(xgboost)
library(stringr)
library(GeneralisedCovarianceMeasure)
library(pcalg)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/GCM.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/PCSelect_Parallel_Update.R")

# Merge the two importance matrices
merge_importance <- function(rank1, rank2) {
    # Get union of all features
    all_features <- union(rank1$Feature, rank2$Feature)
    
    # Create result dataframe with all features
    result <- data.frame(Feature = all_features,
                        Gain = numeric(length(all_features)),
                        Cover = numeric(length(all_features)),
                        Frequency = numeric(length(all_features)))
    
    # Fill in values
    for (i in 1:nrow(result)) {
        feature <- result$Feature[i]
        idx1 <- which(rank1$Feature == feature)
        idx2 <- which(rank2$Feature == feature)
        
        # If feature exists in both matrices, take average
        if (length(idx1) > 0 && length(idx2) > 0) {
            result$Gain[i] <- mean(c(rank1$Gain[idx1], rank2$Gain[idx2]))
            result$Cover[i] <- mean(c(rank1$Cover[idx1], rank2$Cover[idx2]))
            result$Frequency[i] <- mean(c(rank1$Frequency[idx1], rank2$Frequency[idx2]))
        }
        # If feature only exists in first matrix
        else if (length(idx1) > 0) {
            result$Gain[i] <- rank1$Gain[idx1]
            result$Cover[i] <- rank1$Cover[idx1]
            result$Frequency[i] <- rank1$Frequency[idx1]
        }
        # If feature only exists in second matrix
        else if (length(idx2) > 0) {
            result$Gain[i] <- rank2$Gain[idx2]
            result$Cover[i] <- rank2$Cover[idx2]
            result$Frequency[i] <- rank2$Frequency[idx2]
        }
    }
    
    return(result)
}

scRNA <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/PBMC/PBMC_scRNA_exp_5810k_5.csv") # dim(scRNA): 222077  21487
scRNA <- scRNA[,-1] # the 1st column is index
meta_origi <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/PBMC/PBMC_scRNA_meta_5810k_5.csv")
# celltypes_all <- as.vector(unlist(meta$cell_type))
celltypes_all_origi <- as.vector(unlist(meta_origi$celltype))
CD4_index <- which(celltypes_all_origi=="CD4Tcells")
CD8_index <- which(celltypes_all_origi=="CD8Tcells")
T_index <- c(CD4_index,CD8_index)
# NK_index <- which(celltypes_all_origi=="NKcells")
# T_index <- c(CD4_index,CD8_index,NK_index)
valid_index <- T_index
# valid_index <- which(celltypes_all_origi!="Unknown")
scRNA <- scRNA[valid_index,]
celltypes_all <- celltypes_all_origi[valid_index]
meta <- meta_origi[valid_index,]
# Select genes that overlap with UKBB datasets
# UKBB_genes_dict <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/06_ukbb/ukbb_tissue_genes/ukbb_genes_cortex.csv")
# scRNA_genes <- colnames(scRNA)
# UKBB_genes <- as.vector(unlist(UKBB_genes_dict$gene))
# common_genes <- intersect(scRNA_genes,UKBB_genes)
# scRNA <- subset(scRNA,select=common_genes)
# Revised by Liangying Yin, may 7, 2025
# Make sure genes are consistent between scRNA-seq and bulk
# Select genes that overlap with bulk datasets
bulk_genes_dict <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/04_real_application/Monaco/monaco_pbmc_T.txt")
scRNA_genes <- colnames(scRNA)
bulk_genes <- colnames(bulk_genes_dict)[-1]
common_genes <- intersect(scRNA_genes,bulk_genes)
scRNA <- subset(scRNA,select=common_genes)
#' @author  Liangying Yin, May 6, 2025
#' @description  Conduct hierachy-aware marker gene identification based on estimated/known cell type hierachy
# CD4_index <- which(celltypes_all=="CD4Tcells")
# CD8_index <- which(celltypes_all=="CD8Tcells")
# T_index <- c(CD4_index,CD8_index)
# celltypes_all[T_index] <- "Tcells"
celltypes <- unique(celltypes_all)
# celltypes
# envir <- meta$IID
# unique(envir)
# sixk_index <- which(envir=="pbmc6k")
# A_index <- which(envir=="pbmcA")
# length(A_index)
# B_index <- which(envir=="pbmcB")
# length(B_index)
# C_index <- which(envir=="pbmcC")
# length(C_index)
# envir_one_index <- c(sixk_index,B_index)
# envir_one_index <- sort(envir_one_index)
# envir_two_index <- c(A_index,C_index)
# envir_two_index <- sort(envir_two_index)
# envir_one_index <-which(envir=="pbmc8k")
# envir_two_index <- which(envir=="pbmc10k")

# Revised by Liangying, April 28,2025
# batch_one_index <- which(meta$batch=="batch1")
# batch_two_index <- which(meta$batch=="batch2")
# batch_three_index <- which(meta$batch=="batch3")
# envir_one_index <- batch_one_index
# envir_two_index <- c(batch_two_index,batch_three_index)
# envir_two_index <- sort(envir_two_index)

# envir_binary <- rep(0,length(celltypes_all))
# envir_binary[envir_two_index] <- 1

data_fir <- scRNA # Use the whole data to identify predictors then use  PCSimple to remove redundant markers
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- celltypes_all

# celltypes <- unique(label_fir)
for(i in 1:length(celltypes)){
    #i <- 4
    celltype <- celltypes[i]
    cat("This is the identification of cell marker for",celltype,"\n")
    
    # train xgboost for the T2D  dataset
    label_fir_binary <- rep(0,length(label_fir))
    celltype_index <- which(label_fir==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the dataset. \n")
    label_fir_binary[celltype_index] <- 1
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified cellmarkers in the dataset.\n")

    merged_rank <- dtrain_fir_rank
    # merged_rank <- dtrain_fir_rank[order(dtrain_fir_rank$Gain,decreasing = TRUE),]
    

    # Y_binary <- label_binary
    assoc_feats <- as.vector(unlist(merged_rank$Feature))
    X_assoc <- scRNA[, ..assoc_feats]  # For data.table
    
    dfnew = data.frame(outcome=label_fir_binary, X_assoc) # binarized celltype labels
    library(coop)
    expr_corMat = pcor(dfnew[,-1])
    #now add back the correlation between the outcome and each feature
    #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)}  )
    precompute_corMat = cbind(cor_with_outcome, expr_corMat)
    precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat)
    n <- nrow (dfnew)
    V <- colnames(dfnew) # labels aka node names
    ## estimate local causal network structure around the response variable
    pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
                                        dm=dfnew[,-1],
                                        method = c("parallel"),
                                        mem.efficient = FALSE,
                                        num_workers = 10,
                                        # alpha=0.001,
                                        alpha=0.05,
                                        corMat = precompute_corMat,
                                        max_ord=3,
                                        corMethod = "standard",
                                        verbose = TRUE, directed = TRUE)
    # The following codes are used to summarize causally relevant genes
    assoc <- pcSimple.fit$G
    genes <- names(pcSimple.fit$G)
    zMin <- pcSimple.fit$zMin
    res <- cbind(genes,assoc,zMin)
    res_causal <- res[which(assoc==TRUE),]
    if(!is.null(dim(res_causal))){
        causal_genes <- res_causal[,1]
    }else{
        causal_genes <- res_causal[1]
    }    
    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/PBMC/Bulk_matched/5810K_5_PCSimple_Hierarchy_CD48/Monaco/pbmc5810k_5_no_unknown_major_cellmarkers_PCSimple_",celltype,".RData")
    save(celltype,causal_genes,causal_genes,res,file=outputfile)
    #save(celltype,dtrain_fir_rank,dtrain_sec_rank,dtrain_thi_rank,dtrain_four_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")

}


library(data.table)
library(xgboost)
library(stringr)
library(GeneralisedCovarianceMeasure)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/GCM.R")

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
data_fir <- scRNA
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- celltypes_all

sample_num <- length(label_fir)

for(i in 1:length(celltypes)){
    celltype <- celltypes[i]
    cat("This is the identification of cell marker for",celltype,"\n")
    
    # train xgboost 
    label_fir_binary <- rep(0,length(label_fir))
    celltype_index <- which(label_fir==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_fir_binary[celltype_index] <- 1
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified cellmarkers in the 1st dataset.\n")
    merged_rank <- dtrain_fir_rank # Here we only train the xgboost for celltype annotation

    # employ IGCM to identify the invariant cellmarkers
    # celltype_all_index <- which(label_binary==celltype)
    # label_binary[celltype_all_index] <- 1

    # Y_binary <- label_binary
    Y_binary <- label_fir_binary
    assoc_feats <- as.vector(unlist(merged_rank$Feature))
    X_assoc <- subset(scRNA,select=assoc_feats)
    initial_assoc_len <- ncol(X_assoc)
    alpha <-0.05
    res_mat <- matrix(1,nrow=ncol(X_assoc),ncol=3)
    colnames(res_mat)<- c("variables","stat","p")
 
    resid_mat <- NULL
    for(j in initial_assoc_len:1){
        X<-envir_binary
        present_index_set <- 1:j
        Z <- subset(X_assoc,select=present_index_set)
        gcm_obj <- gcm_test(X=X,Y=Y_binary,Z=Z,regr.par=list(max_depth=4,max_delta_step=6), regr.method = "xgboost")
        #gcm_obj <- gcm_test(X=X,Y=Y_binary,Z=Z,regr.par=list(max_depth=3,max_delta_step=10), regr.method = "xgboost")
        gcm_obj_stat <- gcm_obj$test.statistic
        gcm_obj_p <- gcm_obj$p.value
        gcm_R <- gcm_obj$R
        resid_mat <- cbind(resid_mat,gcm_R)
        cat("The statistic and pvalue for independence test are respectively",gcm_obj_stat,"and",gcm_obj_p,"\n")
        res_mat[(initial_assoc_len-j+1),1] <- j
        res_mat[(initial_assoc_len-j+1),2] <- gcm_obj_stat
        res_mat[(initial_assoc_len-j+1),3] <- gcm_obj_p
    }
    IGCM_mat <- matrix(0,nrow=nrow(res_mat)-1,ncol=5)
    colnames(IGCM_mat) <- c("base_set","reduced_set","stat_change","stat_change_Z","stat_change_P")
    causal_set_index <- NULL
    change_alpha <- 0.05
    for(k in 1:(nrow(res_mat)-1)){
        base_gcm_stat <- res_mat[k,2]
        next_gcm_stat <- res_mat[(k+1),2]
        R1 <- resid_mat[,k]
        R2 <- resid_mat[,(k+1)]
        gcm_stat_stat_cov <- stat_change_detection(R1,R2,sample_num) # for Tns that follows standard normal distribution, the correlation also equals the covariance
        gcm_stat_change <- abs(next_gcm_stat) - abs(base_gcm_stat)
        gcm_stat_change_var <- 2-4/pi*(gcm_stat_stat_cov*asin(gcm_stat_stat_cov)+sqrt(1-gcm_stat_stat_cov^2))
        if(abs(gcm_stat_change)>0){
        gcm_stat_change_Z <- gcm_stat_change/sqrt(gcm_stat_change_var) # we directly use the calculated change!
        gcm_stat_change_P <- pnorm(gcm_stat_change_Z,lower.tail = FALSE) # Here, we use one tail test
        IGCM_mat[k,1] <- nrow(res_mat)-k+1
        IGCM_mat[k,2] <- nrow(res_mat)-k
        IGCM_mat[k,3] <- gcm_stat_change
        IGCM_mat[k,4] <- gcm_stat_change_Z
        IGCM_mat[k,5] <- gcm_stat_change_P
        cat("The stats change between the first",nrow(res_mat)-k+1,"and",nrow(res_mat)-k,"variables is",gcm_stat_change,"with a significance of",gcm_stat_change_P,"\n")
        if(gcm_stat_change_P<=change_alpha){
           cat("The first",nrow(res_mat)-k+1,"variables is the potential causal variable set.\n")
           causal_index <- nrow(res_mat)-k+1
           causal_set_index <- c(causal_set_index,causal_index)
        }
        }else{
            cat("The stats change between the first",nrow(res_mat)-k+1,"and",nrow(res_mat)-k,"variables is",gcm_stat_change,"\n")
        }
    }       

    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/PBMC/6810KABC_ICP/6810KABC_cellmarkers_",celltype,".RData")
    save(celltype,dtrain_fir_rank,merged_rank,IGCM_mat,causal_set_index,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")
    
}


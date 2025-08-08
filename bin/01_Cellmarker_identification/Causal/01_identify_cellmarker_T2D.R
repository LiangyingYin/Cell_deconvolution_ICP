library(data.table)
library(xgboost)
library(stringr)
library(GeneralisedCovarianceMeasure)
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


scRNA <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/T2D/T2D_scRNA_exp.csv") # dim(scRNA): 222077  21487
scRNA <- scRNA[,-1]
meta <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/02_diabetes_T2D_pancreas/h5ad/diabetes_training_meta.csv")
celltypes_all <- as.vector(unlist(meta$cell_type))
celltypes <- unique(celltypes_all)
envir <- as.vector(unlist(meta$disease_state))
case_index <- which(envir=="T2D")
ctrl_index <- which(envir=="Control")
envir_binary <- rep(0,length(envir))
envir_binary[case_index] <- 1
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
# celltypes <- unique(label_fir)

for(i in 1:length(celltypes)){
    #i <- 4
    celltype <- celltypes[i]
    cat("This is the identification of cell marker for",celltype,"\n")
    
    # # train xgboost for the T2D  dataset
    # label_fir_binary <- rep(0,length(label_fir))
    # celltype_index <- which(label_fir==celltype)
    # cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    # label_fir_binary[celltype_index] <- 1
    # dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    # dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    # dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    # cat("There are",nrow(dtrain_fir_rank),"identified cellmarkers in the 1st dataset.\n")
    # #cellmarker_sec_list[[i]] <- dtrain_sec_rank
    
    # # train xgboost for the control dataset
    # label_sec_binary <- rep(0,length(label_sec))
    # celltype_index <- which(label_sec==celltype)
    # label_sec_binary[celltype_index] <- 1
    # cat("There are",length(celltype_index),celltype,"cells in the 2nd dataset. \n")
    # dtrain_sec <- xgb.DMatrix(data = data_sec, label = label_sec_binary)
    # dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    # dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
    # cat("There are",nrow(dtrain_sec_rank),"identified cellmarkers in the 1st dataset.\n")
    # #cellmarker_thi_list[[i]] <- dtrain_thi_rank
    # merged_rank <- merge_importance(dtrain_fir_rank, dtrain_sec_rank)
    # merged_rank <- merged_rank[order(merged_rank$Gain,decreasing = TRUE),]
    

    # # employ IGCM to identify the invariant cellmarkers
    # celltype_all_index <- which(label_binary==celltype)
    # label_binary[celltype_all_index] <- 1

    # Y_binary <- label_binary
    # assoc_feats <- as.vector(unlist(merged_rank$Feature))
    # X_assoc <- scRNA[, ..assoc_feats]  # For data.table
    
    # outcome_celltype <- outcome[celltype_index]
    # dfnew = data.frame(outcome=label_binary, X_assoc)
    dfnew = data.frame(outcome=label_binary, scRNA)
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
    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/Full_gene_set/T2D/Diabetes_cellmarkers_PCSimple0.05_",celltype,".RData")
    save(celltype,causal_genes,res,file=outputfile)
    #save(celltype,dtrain_fir_rank,dtrain_sec_rank,dtrain_thi_rank,dtrain_four_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")
    
}


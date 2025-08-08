library(data.table)
library(xgboost)
library(stringr)
library(SingleCellExperiment)
library(scuttle)
library(pcalg)
library(GeneralisedCovarianceMeasure)
library(coop)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/GCM.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/PCSelect_Parallel_Update.R")
#source("/mnt/home/yinly/projects/Cell_deconvolution/bin/Github/01_cellmarker_identification_method.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/Github/get_cellmarkers.R")

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






load("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/res/14_addon_simulation/2celltypes/toydata_batch.RData")
normal <- NULL
normal_celltype <- NULL
normal_batch <- NULL
patient <- NULL
patient_celltype <- NULL
patient_batch <- NULL
# We only take the first 4 patients and controls for analysis, as there are only 6 samples
for(i in 1:4){
    subject_normal <- as.matrix(logcounts(normal_list[[i]]))
    subject_normal_meta <- colData(normal_list[[i]])
    subject_normal_genes <- rownames(normal_list[[i]])
    subject_normal_final <- t(subject_normal)
    subject_normal_celltype <- as.vector(unlist(subject_normal_meta$cell_type))
    normal_celltype <- c(normal_celltype,subject_normal_celltype)
    subject_normal_batch <- as.vector(unlist(subject_normal_meta$batch))
    normal_batch <- c(normal_batch,subject_normal_batch)
    normal <- rbind(normal,subject_normal_final)

    subject_patient <- as.matrix(logcounts(patient_list[[i]]))
    subject_patient_genes <- rownames(patient_list[[i]])
    subject_patient_meta <- colData(patient_list[[i]])
    subject_patient_celltype <- as.vector(unlist(subject_patient_meta$cell_type))
    subject_patient_batch <- as.vector(unlist(subject_patient_meta$batch))
    patient_celltype <- c(patient_celltype,subject_patient_celltype)
    patient_batch <- c(patient_batch,subject_patient_batch)
    subject_patient_final <- t(subject_patient)
    patient <- rbind(patient,subject_patient_final)

}

colnames(normal) <- subject_normal_genes
colnames(patient) <- subject_patient_genes
scRNA <- rbind(normal,patient)
colnames(scRNA) <- colnames(normal)
celltype <- c(normal_celltype,patient_celltype)
batch <- c(normal_batch,patient_batch) # unique(batch): CTRL,STIM
# scRNA_all <- cbind(scRNA,as.vector(unlist(celltype)),as.vector(unlist(batch)))
scRNA_all <- cbind(scRNA,celltype,batch)
colnames(scRNA_all) <- c(colnames(normal),"celltype","batch")
scRNA_all <- as.data.frame(scRNA_all)
# scRNA_all$status <- c(rep("CTRL",nrow(normal)),rep("STIM",nrow(patient)))
# Use batch as environment
batch1_index <- which(scRNA_all$batch=="10x Chromium (v2)")
batch1 <- scRNA_all[batch1_index,] 
batch2_index <- which(scRNA_all$batch=="10x Chromium (v3)")
batch2 <- scRNA_all[batch2_index,] 
# get binary environment variable coded as 0,1
envir_binary <- rep(0,nrow(scRNA_all))
envir_binary[batch2_index] <- 1
feat_len <- ncol(scRNA_all) - 2
cat("There are",feat_len,"genes in the simulated scRNA-seq dataset.\n")
# Get environment specific data
# batch 1
data_fir <- subset(batch1,select=(1:feat_len))
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- as.vector(unlist(scRNA_all$celltype))[batch1_index]
# batch 2
data_sec <- subset(batch2,select=(1:feat_len))
data_sec <- as.matrix(data_sec)
class(data_sec) <- "numeric"
label_sec <- as.vector(unlist(scRNA_all$celltype))[batch2_index]
celltypes_all <- scRNA_all$celltype

celltypes <- unique(scRNA_all$celltype)
sample_num <- nrow(scRNA_all)

for(i in 1:length(celltypes)){
    celltype <- celltypes[i]
    cat("This is the identification of cell marker for",celltype,"\n")
    
    # train xgboost for the 1st dataset
    label_fir_binary <- rep(0,length(label_fir))
    celltype_index <- which(label_fir==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_fir_binary[celltype_index] <- 1
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified cellmarkers in the 1st dataset.\n")
    # merged_rank <- dtrain_fir_rank # Here we only train the xgboost for celltype annotation

    # Y_binary <- label_fir_binary
    # assoc_feats <- as.vector(unlist(dtrain_fir_rank$Feature))
    # X_assoc <- subset(data_fir,select=assoc_feats)
    # # employ PCSimple to identify the causal cellmarkers
    # dfnew = data.frame(outcome=label_fir_binary, X_assoc) # binarized celltype labels

    # expr_corMat = pcor(dfnew[,-1])
    # #now add back the correlation between the outcome and each feature
    # #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    # cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)}  )
    # precompute_corMat = cbind(cor_with_outcome, expr_corMat)
    # precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat)
    # n <- nrow (dfnew)
    # V <- colnames(dfnew) # labels aka node names
    # ## estimate local causal network structure around the response variable
    # #library(parallel)
    # pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
    #                                     dm=dfnew[,-1],
    #                                     method = c("parallel"),
    #                                     mem.efficient = FALSE,
    #                                     num_workers = 10,
    #                                     # alpha=0.001,
    #                                     alpha=0.05,
    #                                     corMat = precompute_corMat,
    #                                     max_ord=3,
    #                                     corMethod = "standard",
    #                                     verbose = TRUE, directed = TRUE)
    # # The following codes are used to summarize causally relevant genes
    # assoc <- pcSimple.fit$G
    # genes <- names(pcSimple.fit$G)
    # genes <- gsub("HLA.","HLA-",genes) # There are special cases where the "-" in gene symbol was replaced by ".
    # zMin <- pcSimple.fit$zMin
    # res <- cbind(genes,assoc,zMin)
    # res_causal <- res[which(assoc==TRUE),]
    # if(!is.null(dim(res_causal))){
    #     causal_genes <- res_causal[,1]
    # }else{
    #     causal_genes <- res_causal[1]
    # }
    
    # Train xgboost for the 2nd dataset
    # train xgboost for the 1st dataset
    label_sec_binary <- rep(0,length(label_sec))
    celltype_index <- which(label_sec==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_sec_binary[celltype_index] <- 1
    dtrain_sec <- xgb.DMatrix(data = data_sec, label = label_sec_binary)
    dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
    cat("There are",nrow(dtrain_sec_rank),"identified cellmarkers in the 2nd dataset.\n")
    # merged_rank <- dtrain_fir_rank # Here we only train the xgboost for celltype annotation

    # Y_binary_sec <- label_sec_binary
    # assoc_feats_sec <- as.vector(unlist(dtrain_sec_rank$Feature))
    # X_assoc_sec <- subset(data_sec,select=assoc_feats_sec)
    # # employ PCSimple to identify the causal cellmarkers
    # dfnew = data.frame(outcome=label_sec_binary, X_assoc_sec) # binarized celltype labels
    # # library(coop)
    # expr_corMat = pcor(dfnew[,-1])
    # #now add back the correlation between the outcome and each feature
    # #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    # cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)}  )
    # precompute_corMat = cbind(cor_with_outcome, expr_corMat)
    # precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat)
    # n <- nrow (dfnew)
    # V <- colnames(dfnew) # labels aka node names
    # ## estimate local causal network structure around the response variable
    # #library(parallel)
    # pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
    #                                     dm=dfnew[,-1],
    #                                     method = c("parallel"),
    #                                     mem.efficient = FALSE,
    #                                     num_workers = 10,
    #                                     # alpha=0.001,
    #                                     alpha=0.05,
    #                                     corMat = precompute_corMat,
    #                                     max_ord=3,
    #                                     corMethod = "standard",
    #                                     verbose = TRUE, directed = TRUE)
    # # The following codes are used to summarize causally relevant genes
    # assoc <- pcSimple.fit$G
    # genes <- names(pcSimple.fit$G)
    # genes <- gsub("HLA.","HLA-",genes) # There are special cases where the "-" in gene symbol was replaced by ".
    
    # zMin <- pcSimple.fit$zMin
    # res_sec <- cbind(genes,assoc,zMin)
    # res_causal_sec <- res_sec[which(assoc==TRUE),]
    # if(!is.null(dim(res_causal_sec))){
    #     causal_genes_sec <- res_causal_sec[,1]
    # }else{
    #     causal_genes_sec <- res_causal_sec[1]
    # }
    
    # res_causal_merged <- merge_and_rank_genes(res_causal,res_causal_sec)
    # if(!is.null(dim(res_causal_merged))){
    #     causal_genes_final <- res_causal_merged[,1]
    # }else{
    #     causal_genes_final <- res_causal_merged[1]
    # }
    res_causal_merged <- merge_importance(dtrain_fir_rank,dtrain_sec_rank)
    res_causal_merged <- res_causal_merged[order(res_causal_merged$Gain,decreasing=TRUE),]
    if(!is.null(dim(res_causal_merged))){
        causal_genes_final <- res_causal_merged[,1]
    }else{
        causal_genes_final <- res_causal_merged[1]
    }
    
    # Employing PC
    # Y_binary <- label_binary
    Y_binary <- rep(0,length(celltypes_all))
    celltype_index <- which(celltypes_all==celltype)
    Y_binary[celltype_index] <- 1 # celltype label: one vs remaining all

    X_assoc <- subset(scRNA,select=causal_genes_final)
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
    if(length(causal_set_index)>=1){
        causal_gene_set_index <- causal_set_index[1]
        causal_gene_final <- res_causal_merged$Feature[1:causal_gene_set_index]
    } else{
        causal_gene_final <-NULL
    }
    
    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/Simulation/2Celltypes/2Celltypes_cellmarkers_",celltype,".RData")
    save(celltype,dtrain_fir_rank,dtrain_sec_rank,res_causal_merged,IGCM_mat,causal_genes_final,causal_set_index,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")
    
}


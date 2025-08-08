library(data.table)
library(xgboost)
library(stringr)
library(GeneralisedCovarianceMeasure)
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


scRNA <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/Depression/Depression_scRNA_exp.csv") # dim(scRNA): 222077  21487
scRNA <- scRNA[,-1] # the 1st column is index
meta <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/02_depression_cortex//h5ad/depression_training_meta.csv")
celltypes_all <- as.vector(unlist(meta$CellType_major))
celltypes <- unique(celltypes_all)
patient_type <- as.vector(unlist(meta$patient_type))
disease_state <- rep(0,length(patient_type))
case_index <- which(patient_type=="Suicide")
ctrl_index <- which(patient_type=="Control")
disease_state[case_index] <- 1
# envir <- disease_state
# envir <- as.vector(unlist(meta$disease_state))
# case_index <- which(envir=="T2D")
# ctrl_index <- which(envir=="Control")
# envir_binary <- rep(0,length(envir))
# envir_binary[case_index] <- 1
# sex <- rep(0,length(bmi))
batch <- meta$batch
batch <- as.numeric(gsub("B","",batch))
envir_one_index <- which(batch<=3)
envir_two_index <- which(batch>3) # here we encode female as 1
envir_binary <- rep(0,length(batch))
envir_binary[envir_two_index] <- 1
outcome <- disease_state

for (celltype in celltypes){
    celltype_index <- which(celltypes_all==celltype)
    scRNA_celltype <- scRNA[celltype_index]
    envir_celltype <- envir_binary[celltype_index]
    envir_celltype_one_index <- which(envir_celltype==0)
    envir_celltype_two_index <- which(envir_celltype==1)
    outcome_celltype <- outcome[celltype_index]

    data_fir <- scRNA_celltype[envir_celltype_one_index,]
    data_fir <- as.matrix(data_fir)
    class(data_fir) <- "numeric"
    label_fir <- outcome_celltype[envir_celltype_one_index]

    data_sec <- scRNA_celltype[envir_celltype_two_index,]
    data_sec <- as.matrix(data_sec)
    #data_fir <- as.matrix(scRNA_ctrl[,(1:(feat_len-1))])
    class(data_sec) <- "numeric"
    label_sec <- outcome_celltype[envir_celltype_two_index]
    
    label_binary <- outcome_celltype
    sample_num <- length(celltype_index)

    label_fir_binary <- label_fir
    cat("There are",sum(label_fir_binary),"depression patients in the 1st dataset. \n")
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified trait-associated markers in the 1st dataset.\n")
    #cellmarker_sec_list[[i]] <- dtrain_sec_rank

    # train xgboost for the control dataset
    label_sec_binary <- label_sec
    cat("There are",sum(label_sec_binary),"depression patients in the 2nd dataset. \n")
    dtrain_sec <- xgb.DMatrix(data = data_sec, label = label_sec_binary)
    dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
    cat("There are",nrow(dtrain_sec_rank),"identified trait-associated markers in the 2nd dataset.\n")
    #cellmarker_thi_list[[i]] <- dtrain_thi_rank
    merged_rank <- merge_importance(dtrain_fir_rank, dtrain_sec_rank)
    merged_rank <- merged_rank[order(merged_rank$Gain,decreasing = TRUE),]

    Y_binary <- label_binary
    assoc_feats <- as.vector(unlist(merged_rank$Feature))
    X_assoc <- scRNA_celltype[, ..assoc_feats]  # For data.table
    initial_assoc_len <- ncol(X_assoc)
    alpha <-0.05
    res_mat <- matrix(1,nrow=ncol(X_assoc),ncol=3)
    colnames(res_mat)<- c("variables","stat","p")
    # Y_fir <- Y_binary[ExpInd_fir]
    # Y_sec <- Y_binary[ExpInd_sec]
    resid_mat <- NULL
    
    for(j in initial_assoc_len:1){
        X<-envir_celltype # here the environment variable should be cell-type-specific
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
    celltype <- gsub("/","",celltype)
    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/Depression/Depression_causal_",celltype,".RData")
    save(celltype,dtrain_fir_rank,IGCM_mat,causal_set_index,file=outputfile)
    #save(celltype,dtrain_fir_rank,dtrain_sec_rank,dtrain_thi_rank,dtrain_four_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")

}


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


scRNA <- fread("/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/ObesitySAT/ObesitySAT_scRNA_exp.csv") # dim(scRNA): 222077  21487
scRNA <- scRNA[,-1] # the 1st column is index
meta <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/06_obese_adipose/h5ad/adiposeSAT_training_meta.csv")
# Select genes that overlap with UKBB datasets

UKBB_genes_dict <- fread("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/06_ukbb/ukbb_tissue_genes/ukbb_genes_adiposeSAT.csv")
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
sex <- rep(0,length(bmi))
envir_one_index <- which(meta$sex__ontology_label=="male")
envir_two_index <- which(meta$sex__ontology_label=="female") # here we encode female as 1
sex[envir_two_index] <- 1
envir_binary <- sex
outcome <- disease_state

for (celltype in celltypes){
    celltype_index <- which(celltypes_all==celltype)
    scRNA_celltype <- scRNA[celltype_index]
    envir_celltype <- envir_binary[celltype_index]
    outcome_celltype <- outcome[celltype_index]

    label_binary <- outcome_celltype
    sample_num <- length(celltype_index)

    data_fir <- scRNA_celltype
    data_fir <- as.matrix(data_fir)
    class(data_fir) <- "numeric"
    label_fir_binary <- label_binary

    cat("There are",sum(label_binary),"obesity patients in the dataset. \n")
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified trait-associated markers in the dataset.\n")
    merged_rank <- dtrain_fir_rank
    # merged_rank <- merged_rank[order(merged_rank$Gain,decreasing = TRUE),]

    Y_binary <- label_binary
    assoc_feats <- as.vector(unlist(merged_rank$Feature))
    X_assoc <- scRNA_celltype[, ..assoc_feats]  # For data.table
    
    outcome_celltype <- outcome[celltype_index]
    dfnew = data.frame(outcome=outcome_celltype, X_assoc)
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
                                        alpha=0.001,
                                        # alpha=0.05,
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
    outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/UKBB_gene_set/ObesitySAT/ObesitySAT_PCSimple_causal_",celltype,".RData")
    save(celltype,causal_genes,causal_genes,res,file=outputfile)
    #save(celltype,dtrain_fir_rank,dtrain_sec_rank,dtrain_thi_rank,dtrain_four_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")

}


library(data.table)
library(dplyr)
library(pcalg)
library(coop)
library(WGCNA)
source("/home/yinly/HTE/03_Simulation/src/PCSelect_Parallel_Update.R")

# extract phenotype variable
pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)

# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
PCs_names <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
colnames(PCs)[2:11] <- PCs_names

disorders <- c("Obesity","CAD","T2DM","HTN","Stroke","Heart_Failure")
file_dir <- "/mnt/data/share/yinly/Cell_deconvolution/CTS/all_minmax/"
files <-list.files(path=file_dir,pattern="ukbb_VAT_all_")

target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/Adipose_related_disorders_2024_06_24.csv"
target <- fread(target_dir)

causal_mat <- NULL
tissue <-"VAT"
res_prefix <- "/home/yinly/Cell_deconvolution/03_Causal_Inference_CTS/res/res_causal_genes/"
for (disorder in disorders){
    cat("The causal inference for",disorder,"begins!\n")
    IID <- target$IID
    outcome <- target[,..disorder] # colnames of variables related to depression: Depression,Psychosis(1142),Bipolar(1680),Anxiety(7495)
    outcome_all <- cbind(IID,outcome)
    colnames(outcome_all) <- c("IID","outcome")
    if(disorder=="Obesity"){
        valid_index <- which(!is.na(outcome)) # remove subjects without available bmi
        outcome_all <-  outcome_all[valid_index,]
    }
    disorder_mat <- NULL
    for (file in files){
        expression_file_dir <- paste0(file_dir,file)
        celltype <- gsub("/mnt/data/share/yinly/Cell_deconvolution/CTS/all_minmax/ukbb_VAT_all_","",expression_file_dir)
        celltype <- gsub(".csv","",celltype)
        expr = fread(expression_file_dir )
        
        data <- merge(expr,outcome_all,by="IID",all.Y=TRUE, sort=FALSE) # merge expression with outcome variable
        no_genes <- (ncol(data)-2) # no_genes 

        data = inner_join(data,PCs,by ="IID") # merge with PCs
        X_assoc <- subset(data,select=(2:(no_genes+1)))  

        #********************************
        # Residualized X
        #*********************************
        resid.mat = matrix(nrow = nrow(data), ncol= no_genes )
        attach(data)
        for (gene_index in 1:no_genes) {
            # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
            gene_expr <- as.vector(unlist(subset(X_assoc,select=gene_index)))
            lm.obj = lm(gene_expr ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
            #lm.obj = lm(X_assoc[,i] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
            resid.mat[,gene_index] = lm.obj$resid
        }
        colnames(resid.mat) <- colnames(data)[2:(no_genes+1)]

        #********************************
        # Residualized Y
        #*********************************
        # Revised by yinly on Jan 7,2020,the inner_join function changed the column name of the target_outcome_name
        # # Reset the value of
        # outcome_index = "X30690.0.0"
        lm.obj_y = lm(data$outcome ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        resid_y = lm.obj_y$resid
        detach(data)

        dfnew = data.frame(outcome = resid_y, resid.mat)
        # Revised by yinly: Sep 13, 2023
        # Description: replace pcor in coop with another program cor in WGCNA
        expr_corMat = cor(dfnew[,-1],use = "all.obs", method="pearson",nThreads = 10)
        cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {cor(dfnew[,"outcome"], col)}  )
        # expr_corMat = pcor(dfnew[,-1],inplace=TRUE)
        # cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col,inplace=TRUE)}  )
        precompute_corMat = cbind(cor_with_outcome, expr_corMat)
        precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat )
        cat("Pearson correlation matrix calculation is completed!\n")
        library(parallel)
        #********************************************************
        # Apply PC-simple algorithm around the response variable on NFBC
        #**************************************************************
        pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
                                            dm=dfnew[,-1],
                                            method = c("parallel"),
                                            mem.efficient = FALSE,
                                            num_workers = 10,
                                            alpha=0.05,
                                            #alpha=0.05,
                                            corMat = precompute_corMat,
                                            max_ord=3,
                                            corMethod = "standard",
                                            verbose = TRUE, directed = TRUE)

        outfile_prefix <- paste0("/home/yinly/Cell_deconvolution/03_Causal_Inference_CTS/res/res_",disorder,"/")
        if (!dir.exists(outfile_prefix)) {
            # If the folder doesn't exist, create it
            dir.create(outfile_prefix)
            print(paste("Folder created at:", outfile_prefix))
        } else {
            print(paste("Folder already exists at:", outfile_prefix))
        }
        filepath <- paste0(outfile_prefix,"UKBB_",tissue,"_",disorder,"_",celltype,".Rdata") # nolint
        save(pcSimple.fit, file = filepath)
        cat("The identification of causally relevant genes for",celltype,"is completed!\n")
        
        # The following codes are used to summarize causally relevant genes
        assoc <- pcSimple.fit$G
        genes <- names(pcSimple.fit$G)
        zMin <- pcSimple.fit$zMin
        res <- cbind(genes,assoc,zMin)
        res_filename <- gsub(".Rdata",".csv",filename)
        res_causal <- res[which(assoc==TRUE),]
        cell_causal <- cbind(rep(celltype,nrow(res_causal)),res_causal)
        cat("The dimension of cell_causal is",dim(cell_causal),"\n")
        disorder_mat <- rbind(disorder_mat,cell_causal)
        cat("The dimension of disorder_mat is",dim(disorder_mat),"\n")
        cat("Extraction of causal genes is completed for current celltype! \n")
    }
    colnames(disorder_mat) <- c("celltype",colnames(res_causal))
    disorder_mat_final <- cbind(rep(disorder,nrow(disorder_mat)),disorder_mat)
    cat("The dimension of disorder_mat_final is",dim(disorder_mat_final),"\n")
    causal_mat <- rbind(causal_mat,disorder_mat_final)
}
cat("The dimension of causal_mat is",dim(causal_mat),"\n")
colnames(causal_mat) <- c("disorder","celltype",colnames(res_causal))
fwrite(as.data.frame(causal_mat),paste0(res_prefix,tissue,"_causal_genes.csv"),sep="\t")


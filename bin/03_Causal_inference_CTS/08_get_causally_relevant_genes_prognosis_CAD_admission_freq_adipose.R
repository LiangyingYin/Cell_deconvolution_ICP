library(data.table)
library(dplyr)
#library(RNOmni)
#library(coop)
library(pcalg)
library(WGCNA)
source("/home/yinly/HTE/03_Simulation/src/PCSelect_Parallel_Update.R")

tissue <- "VAT"
disorder <- "Obesity"
file_dir <- "/mnt/data/share/yinly/Cell_deconvolution/CTS/CTS_B-mode/"
files <-list.files(path=file_dir,pattern="ukbb_VAT_")

prognosis_filepath <- "/home/yinly/Cell_deconvolution/03_Analyses/data/CAD_admission_freq.csv"
prognosis <- fread(prognosis_filepath) # 3 cloumns eid, CAD_hos, CAD_mortality
# prognosis_selected <- prognosis[,-2] # extract CAD_mortality 
# colnames(prognosis_selected) <- c("IID","CAD_mortality")
colnames(prognosis) <- c("IID","admission_freq")
prog <- prognosis


pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)
IID <- pheno$FID
# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
colnames(PCs)[2:11] <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")

target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/Adipose_related_disorders_2024_06_24.csv"
target <- fread(target_dir)
IID <- target$IID
outcome <- target[,..disorder] # colnames of variables related to depression: Depression,Psychosis(1142),Bipolar(1680),Anxiety(7495)
outcome_all <- cbind(IID,outcome)
colnames(outcome_all) <- c("IID","outcome")
if(disorder=="Obesity"){
    valid_index <- which(!is.na(outcome)) # remove subjects without available bmi
    outcome_all <-  outcome_all[valid_index,]
}

# Extract CAD-related mortality for further analysis
Obesity_outcome <- merge(outcome_all,prog,by="IID",sort=FALSE)

selected_index <- which(Obesity_outcome$outcome==1) # analysis CAD-related mortality in Obesity patients
outcome_all <- Obesity_outcome[selected_index,-2]
colnames(outcome_all) <- c("IID","outcome")
 

for(j in 1:length(files)){
    # var_thres <- 0.01
    expression_file_dir <- paste0(file_dir,files[j])
    # expression_file_dir = "/home/yinly/Cell_deconvolution/02_Cell_deconvolution/res/all/ukbb_adiposeVAT_all_ASPC.csv" 
    # celltype <- gsub("/home/yinly/Cell_deconvolution/02_Cell_deconvolution/res/all/ukbb_VAT_","",expression_file_dir)
    celltype <- gsub("ukbb_VAT_","",files[j])
    celltype <- gsub(".csv.gz","",celltype)
    expr = fread(expression_file_dir )

       #************************************
    # screen away genes with very low variance (this will lead to NA entries after scaling in huge
    #***********************************
    # var_list <- apply(expr[,2:ncol(expr)], 2, var)
    # gene_expr <- data.frame(expr[,2:ncol(expr)])
    # gene_expr_filtered <-  gene_expr[,var_list>var_thres]
    # expr <- cbind(expr[,1],gene_expr_filtered )
    data <- merge(expr,outcome_all,by="IID",all.Y=TRUE, sort=FALSE)

    # pheno_dir <- "/home/yinly/Cell_deconvolution/03_Analyses/data/UKBB_Extracted_Traits_30April.txt"
    # pheno <- fread(pheno_dir)

    no_genes <- (ncol(data)-2) # no_genes 



    data = inner_join(data,PCs,by ="IID")
    X_assoc <- subset(data,select=(2:(no_genes+1)))  

    #********************************
    # Residualized X
    #*********************************
    resid.mat = matrix(nrow = nrow(data), ncol= no_genes )
    attach(data)
    for (i in 1:no_genes) {
        # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
        gene_expr <- as.vector(unlist(subset(X_assoc,select=i)))
        lm.obj = lm(gene_expr ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        #lm.obj = lm(X_assoc[,i] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        resid.mat[,i] = lm.obj$resid
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
    # for test only
    # dfnew = data.frame(outcome = data$outcome, X_assoc)
    #dfnew_rankNorm = dfnew
    #gc()
    #dfnew_rankNorm = dfnew
    # Revised by yinly: Sep 13, 2023
    # Description: replace pcor in coop with another program cor in WGCNAf
    #weights <- matrix(1,nrow=(nrow(dfnew)-1),ncol=(ncol(dfnew)-1))
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

    outfile_prefix <- "/home/yinly/Cell_deconvolution/03_Causal_Inference_CTS/res/res_CAD_admission_freq_prognosis_VAT/"
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

}



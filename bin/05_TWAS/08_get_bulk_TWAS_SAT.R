library(data.table)
library(dplyr)

tissue <- "SAT"
# extract phenotype variable
pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)
# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
PCs_names <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
colnames(PCs)[2:11] <- PCs_names

# target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/Pancreas_related_disorders_2024_01_23.csv"
target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/CCSR-new-onset-COVID-19-seq-history-new-hospitalization-covariates-whole-UKBB-eid-Dec-23.csv"
target <- fread(target_dir)
colnames(target)[1] <- "IID"

# expr_dir <- paste0("/mnt/data/share/yinly/UKBB/UKBB_",tissue,"_predicted_expression.txt")
expr_dir <- paste0("/mnt/data/share/yinly/UKBB/UKBB_Adipose_Subcutaneous_predicted_expression.txt")
expr <- fread(expr_dir)

# This is for frontal related disorders
# disorders <- c("Anxiety","Bipolar","Depression","Psychosis")
disorders <- c("CAD","T2DM","HTN","Stroke","Heart_Failure")
for(disorder in disorders){
    #colnames(pheno)[1] <-"IID"
    # IID <- pheno$FID
    # Outcome <- pheno[,..disorder] 
    IID <- target$IID
    Outcome <- target[,..disorder] # colnames of variables related to depression: T1DM,Prediabetes
    Outcome <- ifelse(Outcome=="no_hx_no_new", 0, 1)
    outcome_all <- cbind(IID,Outcome)
    colnames(outcome_all) <- c("IID","outcome")

    data <- merge(expr,outcome_all,by="IID",all.Y=TRUE, sort=FALSE)
    no_genes <- (ncol(data)-2) # no_genes 

    data = inner_join(data,PCs,by ="IID")
    X_assoc <- subset(data,select=(2:(no_genes+1))) 

    PCs_active <- subset(data,select=PCs_names)
    Y <- data$outcome
    #********************************
    # Perform univariate test
    #*********************************
    res_mat <- matrix(nrow = no_genes, ncol= 3 ) # 3 columns respectively for gene name, coefficients and pvalue
    # attach(data)
    for (i in 1:no_genes) {
        # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
        X <- as.vector(unlist(subset(X_assoc,select=i)))
        input <- cbind(Y,X,PCs_active)
        lm_obj <- lm( Y ~ ., data=input )
        # lm_obj <- glm( data$outcome ~  gene_expr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        res_mat[i,1] <- colnames(X_assoc)[i] # get gene symbol for the target gene
        res_mat[i,2] <- summary(lm_obj)$coeff[2,1]
        res_mat[i,3] <- summary(lm_obj)$coeff[2,4]
    }
    # detach(data)
    #colnames(res_mat) <- c("Genes","Estimates","Pvalues")
    Pvalue <- res_mat[,3]
    Pvalue_Bonferroni <- p.adjust(Pvalue,method ="bonferroni",n = length(Pvalue))
    Pvalue_BH <- p.adjust(Pvalue,method ="BH",n = length(Pvalue))
    res_new <- cbind(res_mat,Pvalue_Bonferroni,Pvalue_BH)
    colnames(res_new) <- c("Genes","Estimates","Pvalues","Pvalue_Bonferroni","Pvalue_BH")
    tissue <- gsub("Brain_","",tissue)
    outfile_prefix <- paste0("/home/yinly/Cell_deconvolution/05_TWAS/res/res_",tissue,"_bulk/")

    if (!dir.exists(outfile_prefix)) {
        # If the folder doesn't exist, create it
        dir.create(outfile_prefix)
        print(paste("Folder created at:", outfile_prefix))
    } else {
        print(paste("Folder already exists at:", outfile_prefix))
    }
        
    fwrite(as.data.frame(res_new),paste0(outfile_prefix,disorder,"_TWAS.csv"),sep="\t")
    cat("The univariate test for",disorder, "is completed!\n")
}
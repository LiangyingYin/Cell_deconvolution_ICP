library(data.table)
library(dplyr)

tissue <- "Pancreas"
# extract phenotype variable
# pheno_dir <- "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_30April.txt" 
pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)
# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
PCs_names <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
colnames(PCs)[2:11] <- PCs_names

target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/Pancreas_related_disorders_2024_01_23.csv"
# target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/CCSR-new-onset-COVID-19-seq-history-new-hospitalization-covariates-whole-UKBB-eid-Dec-23.csv"
target <- fread(target_dir)
colnames(target)[1] <- "IID"

expr_dir <- "/mnt/data/share/yinly/UKBB/UKBB_Pancreas_predicted_expression.txt"
expr <- fread(expr_dir) # including FID and IID columes, they both indicating the subject IDs
expr <- expr[,-1]

ukbb_cell_prop <- fread("/mnt/data/share/yinly/Cell_deconvolution/Cell_composition/ukbb_diabetes_pancreas_all.tsv.gz")
cellnames <- colnames(ukbb_cell_prop)[2:ncol(ukbb_cell_prop)]
cellnames_new <- gsub(" ","_",cellnames)
colnames(ukbb_cell_prop) <- c("IID",cellnames_new)
selected_celltypes <- cellnames_new[-1] # always eliminate the first celltype when perform and adjustment

disorders <- c("Prediabetes","T1DM","T2DM","T2DM_admission","T2DM_insulin")
# disorders <- c("CAD","T2DM","HTN","Stroke","Heart_Failure")
for(disorder in disorders){
    #********************************************************
    # use target file for VAT,SAT and Pancreas
    #********************************************************
    IID <- target$IID
    Outcome <- target[,..disorder] # colnames of variables related to depression: T1DM,Prediabetes

    # Outcome <- ifelse(Outcome=="no_hx_no_new", 0, 1)   # only necessary for VAT/SAT related diseases

    outcome_all <- cbind(IID,Outcome)
    colnames(outcome_all) <- c("IID","outcome")

    data <- merge(expr,outcome_all,by="IID",all.Y=TRUE, sort=FALSE)
    no_genes <- (ncol(data)-2) # no_genes 

    data = inner_join(data,PCs,by ="IID")
    X_assoc <- subset(data,select=(2:(no_genes+1))) 

    # merge expression profiles, PCs, and cell proportion dataframes together
    data <- inner_join(data,ukbb_cell_prop,by ="IID")

    PCs_active <- subset(data,select=PCs_names)
    Y <- data$outcome
    Cell_prop <- subset(data,select=selected_celltypes) # select the preserved cell proportions
    #********************************
    # Perform univariate test
    #*********************************
    res_mat <- matrix(nrow = no_genes, ncol= 11 ) # 11 columns 
    # attach(data)
    for (i in 1:no_genes) {
        # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
        X <- as.vector(unlist(subset(X_assoc,select=i)))
        input <- cbind(Y,X,PCs_active)
        lm_obj <- lm( Y ~ ., data=input )
        # lm_obj <- glm( data$outcome ~  gene_expr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        res_mat[i,1] <- colnames(X_assoc)[i] # get gene symbol for the target gene
        res_mat[i,2] <- summary(lm_obj)$coeff[2,1] # coefficient
        res_mat[i,3] <- summary(lm_obj)$coeff[2,4] # p-value
        res_mat[i,4] <- summary(lm_obj)$coeff[2,2] # standard error
        error_var <- var(lm_obj$resid) # error variance
        res_mat[i,5] <- error_var
        input_cellprop <- cbind(Y,X,PCs_active,Cell_prop)
        lm_obj_cellprop <- lm( Y ~ ., data=input_cellprop )
        res_mat[i,6] <- summary(lm_obj_cellprop)$coeff[2,1] # coefficient
        res_mat[i,7] <- summary(lm_obj_cellprop)$coeff[2,4] # p-value
        res_mat[i,8] <- summary(lm_obj_cellprop)$coeff[2,2] # standard error
        error_var_cellprop <- var(lm_obj_cellprop$resid) # error variance
        res_mat[i,9] <- error_var_cellprop
        coeff_diff <- summary(lm_obj)$coeff[2,1] - summary(lm_obj_cellprop)$coeff[2,1] # coefficient change small model -full model
        res_mat[i,10] <- coeff_diff
        coeff_diff_var <- summary(lm_obj_cellprop)$coeff[2,2]^2 - summary(lm_obj)$coeff[2,2]^2*error_var_cellprop/error_var # formula 15 in Statistical Methods for Comparing Regression Coefficients between Models
        coeff_diff_se <- sqrt(coeff_diff_var)
        coeff_diff_Z <- coeff_diff/coeff_diff_se
        coeff_diff_P <- 2*pnorm(-abs(coeff_diff_Z))
        res_mat[i,11] <- coeff_diff_P # Pvalue for coefficient change test

    }
    # detach(data)
    Pvalue <- res_mat[,3]
    Pvalue_Bonferroni <- p.adjust(Pvalue,method ="bonferroni",n = length(Pvalue))
    Pvalue_BH <- p.adjust(Pvalue,method ="BH",n = length(Pvalue))

    Pvalue_cellprop <- res_mat[,7]
    Pvalue_cellprop_Bonferroni <- p.adjust(Pvalue_cellprop,method ="bonferroni",n = length(Pvalue_cellprop))
    Pvalue_cellprop_BH <- p.adjust(Pvalue_cellprop,method ="BH",n = length(Pvalue_cellprop))

    res_new <- cbind(res_mat,Pvalue_Bonferroni,Pvalue_BH,Pvalue_cellprop_Bonferroni,Pvalue_cellprop_BH)
    colnames(res_new) <- c("Genes","Estimates","Pvalues","SD","Error_variance","Estimates_cellprop","Pvalues_cellprop","SD_cellprop","Error_variance_cellprop",
                           "Coeff_diff","Coeff_diff_P","Pvalue_Bonferroni","Pvalue_BH","Pvalue_Bonferroni_cellprop","Pvalue_BH_cellprop")

    outfile_prefix <- paste0("/home/yinly/Cell_deconvolution/05_TWAS/res/res_",tissue,"_bulk_TWAS/")  
    
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
library(data.table)
library(dplyr)

# extract phenotype variable
pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)

# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
PCs_names <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
colnames(PCs)[2:11] <- PCs_names

file_dir <- "/mnt/data/share/yinly/Cell_deconvolution/CTS/CTS_B-mode/"
files <-list.files(path=file_dir,pattern="ukbb_pancreas_all_")

# disorders <- c("Depression","Anxiety","Bipolar","Psychosis")
# disorders <- c("T1DM","T2DM","Prediabetes","T2DM_insulin","T2DM_admission")
disorders <- "T2DM"
target_dir <- "/mnt/data/share/yinly/UKBB/Phenotypes/Pancreas_related_disorders_2024_01_23.csv"
target <- fread(target_dir)

# res_prefix <-"/home/yinly/Cell_deconvolution/05_TWAS/res_subgroup/res_CTS_TWAS/"
# TWAS_mat <- NULL
# sig_alpha <- 0.05 # the alpha used to define significantly associated genes
tissue <- "Pancreas"

for (disorder in disorders){
    cat("The TWAS analysis for",disorder,"begins!\n")
    IID <- target$IID
    outcome <- target[,..disorder] # colnames of variables related to depression: Depression,Psychosis(1142),Bipolar(1680),Anxiety(7495)
    na_index <- which(is.na(outcome))
    target[na_index] <- 0

    outcome_all <- cbind(IID,outcome)
    colnames(outcome_all) <- c("IID","outcome")

    for(j in 1:length(files)){
        # get expression file path
        expression_file_dir <- paste0(file_dir,files[j])
        celltype <- gsub("/mnt/data/share/yinly/Cell_deconvolution/CTS/CTS_B-mode/ukbb_pancreas_","",expression_file_dir)
        celltype <- gsub(".csv","",celltype)
        expr = fread(expression_file_dir )

        data <- merge(expr,outcome_all,by="IID", sort=FALSE)
        no_genes <- (ncol(data)-2) # no_genes 

        data = inner_join(data,PCs,by ="IID")
        X_assoc <- subset(data,select=(2:(no_genes+1))) 

        PCs_active <- subset(data,select=PCs_names)
        Y <- data$outcome
        cat("There are",length(which(Y==1)),disorder,"patients for",celltype,"\n")
        #********************************
        # Perform univariate test
        #*********************************
        res_mat <- matrix(nrow = no_genes, ncol= 3 ) # 3 columns respectively for gene name, coefficients and pvalue
        # attach(data)
        for (i in 1:no_genes) {
            # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
            X <- as.vector(unlist(subset(X_assoc,select=i)))
            input <- cbind(Y,X,PCs_active)
            lm_obj <- lm( Y ~ ., data=input ) # use linear regression model, which is much faster than logistic regression model
            #attach(input)
            # lm_obj <- glm( data$outcome ~  gene_expr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
            res_mat[i,1] <- colnames(X_assoc)[i] # get gene symbol for the target gene
            res_mat[i,2] <- summary(lm_obj)$coeff[2,1]
            res_mat[i,3] <- summary(lm_obj)$coeff[2,4]
        }
        # detach(data)
        Pvalues_BH <- p.adjust(as.numeric(res_mat[,3]),method ="BH",n = length(res_mat[,3]))
        Pvalue_Bonferroni <- p.adjust(as.numeric(res_mat[,3]),method ="bonferroni",n = length(res_mat[,3]))
        res_mat <- cbind(res_mat,Pvalues_BH,Pvalue_Bonferroni)
        colnames(res_mat) <- c("Genes","Estimates","Pvalues","Pvalues_BH","Pvalue_Bonferroni")
        cat("The univariate test for",celltype, "is completed!\n")
        outfile_prefix <- paste0("/home/yinly/Cell_deconvolution/04_TWAS_CTS/res/res_",tissue,"_",disorder,"/")
        if (!dir.exists(outfile_prefix)) {
            # If the folder doesn't exist, create it
            dir.create(outfile_prefix)
            print(paste("Folder created at:", outfile_prefix))
        } else {
            print(paste("Folder already exists at:", outfile_prefix))
        }
        fwrite(as.data.frame(res_mat),paste0(outfile_prefix,celltype,"_TWAS.csv"),sep="\t")
    }
    cat("The TWAS analysis for",disorder,"finishes!\n\n")
}
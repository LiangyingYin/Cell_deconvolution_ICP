library(data.table)
library(biglm)
library(dplyr)
library(stringr)


expression_file_dir = "/mnt/data/yinly/UKBB/Result/UKBB_Pancreas_predicted_expression.txt" 
expr <- fread(expression_file_dir)
expr <- expr[,-1]

pheno_dir <- "/mnt/data/yinly/UKBB/UKBB_Extracted_Traits_30April.txt" 
pheno <- fread(pheno_dir)
covarite_filepath <- "/home/yinly/HTE/Data/CCSR-new-onset-COVID-19-seq-history-new-hospitalization-covariates-whole-UKBB-eid-Dec-23.csv"
covariate <- fread(covarite_filepath)
IID <- covariate$eid
T2DM <- covariate$T2DM
unique(T2DM)
outcome <- rep(1,nrow(covariate))
ctrl_index <- which(T2DM=="no_hx_no_new")
outcome[ctrl_index] <- 0
outcome_all <- cbind(IID,outcome)
colnames(outcome_all) <- c("IID","outcome")
class(outcome_all) <-"numeric"

data <- merge(expr,outcome_all,by="IID",sort=FALSE)
no_genes <- (ncol(data)-2) # no_genes 
#feat_len <- ncol(data)
#X_assoc <- subset(data,select=(1:(feat_len-1)))

# adjust for population stratifications
PC_start_ind = which( colnames(pheno) == "22009-0.1")
PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
colnames(PCs)[1] <- "IID"
colnames(PCs)[2:11] <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")

data = inner_join(data,PCs,by ="IID")
X_assoc <- subset(data,select=(2:(no_genes+1)))  
univariate_mat = matrix(nrow = no_genes, ncol= 3 )
attach(data)
for (i in 1:no_genes) {
    # Commented by:yinly
    # Perform univariate test by linear regression to select significant genes
    gene_expr <- as.vector(unlist(subset(X_assoc,select=i)))
    lm.obj = biglm(data$outcome ~ gene_expr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,data=data)
    #resid.mat[,i] = lm.obj$resid
    univariate_mat[i,1] <- colnames(X_assoc)[i]
    univariate_mat[i,c(2,3)] <- summary(lm.obj)$mat[2,c(1,5)]
    # univariate_mat[i,2] = summary(lm.obj)$coeff[2,1]
    # univariate_mat[i,3] = summary(lm.obj)$coeff[2,4]
}
detach(data)
colnames(univariate_mat) = c("Genes","Estimates","Pvalues")
Qvalues <- p.adjust(univariate_mat[,3],method="fdr",n=nrow(univariate_mat))
univariate_mat_final <- cbind(univariate_mat,Qvalues)
#out_filepath_prefix <- "/home/yinly/Cell_deconvolution/04_Association_analyses/res/res_assoc_Depression"
out_filepath_prefix <- "/home/yinly/Cell_deconvolution/04_Association_analyses/res/res_assoc_T2DM/"
out_filepath <- paste0(out_filepath_prefix,"UKBB_Pancreas_T2DM.csv") # nolint
fwrite(as.data.frame(univariate_mat_final),file=out_filepath,sep="\t")
cat("The identification of susceptibility genes is completed!\n")

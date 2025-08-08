library(data.table)
library(dplyr)

tissue <- "SAT"

file_prefix <- paste0("/home/yinly/Cell_deconvolution/05_TWAS/res/res_",tissue,"_bulk_TWAS/")
files <- list.files(file_prefix)
sum_mat <- matrix(0,nrow=length(files),ncol=6)
colnames(sum_mat) <- c("Disorder","No.of genes with sig change","No.of sig genes(original)","No.of sig genes(original) with sig change","No.of sig genes(adjusted)","No.common sig genes")
BH_alpha <- 0.1 # alpha define the cutoff for significantly associated genes
for(i in 1:length(files)){
    file_name <- files[i]
    res <- fread(paste0(file_prefix,file_name))
    disorder <- gsub("_TWAS.csv","",file_name)
    res_sig <- res[res$Coeff_diff_P<=0.05,] # only summarize genes with significant coefficient changes
    sum_mat[i,1] <- disorder
    sum_mat[i,2] <- nrow(res_sig)
    res_origi_sig <- res[res$Pvalue_BH<=BH_alpha,]
    sum_mat[i,3] <- nrow(res_origi_sig)
    res_origi_sig_changed <- res_origi_sig[res_origi_sig$Coeff_diff_P<=0.05,]
    sum_mat[i,4] <- nrow(res_origi_sig_changed)
    res_adjust_sig <- res[res$Pvalue_BH_cellprop<=BH_alpha,]
    sum_mat[i,5] <- nrow(res_adjust_sig)
    common_genes <- intersect(unlist(res_origi_sig$Genes),unlist(res_adjust_sig$Genes))
    common_genes_final <- setdiff(common_genes,"FID")
    sum_mat[i,6] <- length(common_genes_final)
}
outfile_prefix <- "/home/yinly/Cell_deconvolution/05_TWAS/res/"
outfile <- paste0(tissue,"_coeff_change_summary_BH_alpha_0.1.csv")
fwrite(as.data.frame(sum_mat),paste0(outfile_prefix,outfile),sep="\t")

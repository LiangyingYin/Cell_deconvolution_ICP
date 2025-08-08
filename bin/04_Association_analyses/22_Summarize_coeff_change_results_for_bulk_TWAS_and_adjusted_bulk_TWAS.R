library(data.table)

tissues <- c("Frontal_Cortex","Pancreas","SAT","VAT")
file_prefix <- "/home/yinly/Cell_deconvolution/05_TWAS/res/"

Overall_sig_mat <- NULL
for(tissue in tissues){
    tissue_sig_mat <- NULL
    filepath <- paste0(file_prefix,"res_",tissue,"_bulk_TWAS/")
    cat("This is the correction for",tissue,"\n")
    files <- list.files(filepath)
    for(file in files){
        disorder <- gsub("_TWAS.csv","",file)
        cat("This is the correction for",disorder,"\n")
        twas <- fread(paste0(filepath,file))
        Coeff_diff_P <- twas$Coeff_diff_P
        sig_index <- which(Coeff_diff_P<=0.05)
        twas_sig <-twas[sig_index,] 
        disorder_mat <- cbind(rep(disorder,nrow(twas_sig)),twas_sig)
        tissue_sig_mat <- rbind(tissue_sig_mat,disorder_mat)
        cat("This correction for",disorder,"is completed!\n")
    }
    tissue_mat_new <- cbind(rep(tissue,nrow(tissue_sig_mat)),tissue_sig_mat)
    Overall_sig_mat <- rbind(Overall_sig_mat,tissue_mat_new)
}
colnames(Overall_sig_mat) <- c("Tissue","Disorder",colnames(twas))
fwrite(as.data.frame(Overall_sig_mat),"/home/yinly/Cell_deconvolution/05_TWAS/res/res_coef_change_bulk_TWAS/Coef_change_detection_results.csv",sep="\t")
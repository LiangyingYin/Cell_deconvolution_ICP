library(data.table)
library(qpcR)
# library(pcalg)

disorders <- c("Anxiety","Bipolar","Depression","Psychosis")

tissue <-"Frontal_Cortex"
TWAS_mat <- NULL
sig_alpha <- 0.05 # the alpha used to define significantly associated genes
res_prefix <-paste0("/home/yinly/Cell_deconvolution/05_TWAS/res/res_CTS_TWAS/")

for(j in 1:length(disorders)){
    disorder <- disorders[j]
    filepath <- paste0("/mnt/data/share/yinly/Cell_deconvolution/05_TWAS/res/res_",disorder,"/")
    
    # if (!dir.exists(res_prefix)) {
    #     # If the folder doesn't exist, create it
    #     dir.create(res_prefix)
    #     print(paste("Folder created at:", res_prefix))
    # } else {
    #     print(paste("Folder already exists at:", res_prefix))
    # }
    Genes <- NULL
    files <- list.files(filepath)
    res_summary <- matrix(0,nrow=length(files),ncol=2)
    for(i in 1:length(files)){
        filename <- files[i]
        celltype <- gsub("_TWAS.csv","",filename)
        cat("This is the extraction of CTS associated gene for",celltype,"\n")
        res_TWAS <- fread(paste0(filepath,filename))
        res_TWAS_sig <- res_TWAS[which(Pvalue_BH<=sig_alpha),]
        res_summary[i,1] <- celltype
        res_summary[i,2] <- nrow(res_TWAS_sig)
        Genes <- c(Genes,as.vector(unlist(res_TWAS_sig[,1])))
        cat("Extraction of causal genes is completed for current celltype! \n")
    }
    colnames(res_summary) <- c("Celltype","No.")
    unique_genes <- unique(Genes)
    TWAS_mat <- qpcR:::cbind.na(TWAS_mat, Genes)
    cat("There are",(length(Genes)-length(unique_genes)),"common genes \n")
    cat("There are",length(unique_genes),"unique genes \n")
    filename <- paste0(tissue,"_",disorder,".csv")
    fwrite(as.data.frame(res_summary),paste0(res_prefix,filename),sep="\t")
    cat("Extraction of causal genes is completed for",disorder,"\n\n")
}
print(dim(TWAS_mat))
colnames(TWAS_mat) <- c("NULL",disorders)
out_path <- paste0("/home/yinly/Cell_deconvolution/05_TWAS/res/res_CTS_TWAS_genes_details/",tissue,"_TWAS_genes_details_Mar28_Linear.csv")
fwrite(as.data.frame(TWAS_mat),out_path,sep="\t")

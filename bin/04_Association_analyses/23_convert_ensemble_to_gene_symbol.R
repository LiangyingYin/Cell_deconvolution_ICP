library(data.table)
library(stringr)

coef_twas <- fread("/home/yinly/Cell_deconvolution/05_TWAS/res/res_coef_change_bulk_TWAS/Coef_change_detection_results.csv")
dict <- fread("/mnt/data/share/yinly/Ensemble_to_Genesymbol/dict.csv")
genes <- coef_twas$Genes
genes_split <- str_split_fixed(genes, "\\.", 2)
colnames(genes_split) <- c("gene_id","version")
genes_all <- merge(genes_split,dict,by="gene_id",all.x=TRUE)
coef_twas_final <- cbind(coef_twas,genes_all)
fwrite(as.data.frame(coef_twas_final),"/home/yinly/Cell_deconvolution/05_TWAS/res/res_coef_change_bulk_TWAS/Coef_change_detection_results_gene_symbol.csv")

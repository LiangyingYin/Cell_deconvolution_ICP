library(data.table)

filepath <- "/mnt/data/yinly/Cell_deconvolution/CTS/lung/"

CTS_all <- NULL
for(i in 1:48){
    filename <- paste0("ukbb_lung_common_Bcells_",i,".csv.gz")
    CTS <- fread(paste0(filepath,filename))
    cat("This is the",i,"file for CTS profile.\n")
    CTS <- t(CTS)
    CTS_split <- CTS[-1,]
    CTS_all <- rbind(CTS_all,CTS_split)
}
#colnames(CTS_all)<- CTS_split[1,] # previously, we wrongly named the genes, so we need to change them back.
colnames(CTS_all)<- CTS[1,]
colnames(CTS_all)[1] <- "IID"
outfilepath <- "/mnt/data/yinly/Cell_deconvolution/CTS/all/"
fwrite(CTS_all,paste0(outfilepath,"ukbb_lung_common_Bcells.csv"),sep="\t")
cat("The processing of data files is completed!\n")
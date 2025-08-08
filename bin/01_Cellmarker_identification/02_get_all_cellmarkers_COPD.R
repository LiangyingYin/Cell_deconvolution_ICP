library(data.table)

filepath <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_COPD/"
files <- list.files(filepath)
celltypes <- NULL
cellmarkers_common <- vector(mode = "list", length = length(files))
cellmarkers_all <- vector(mode = "list", length = length(files))
for(i in 1:length(files)){
    filename <- files[i]
    celltype <- gsub("COPD_cellmarkers_","",filename)
    #celltype <- gsub("COVID_Lung_cellmarkers_","",filename)
    celltype <- gsub(".RData","",celltype)
    celltypes <- c(celltypes,celltype)
    filename_full <- paste0(filepath,filename)
    load(filename_full)
    cellmakers_fir <- as.vector(unlist(dtrain_fir_rank[,1]))
    cellmakers_sec <- as.vector(unlist(dtrain_sec_rank[,1]))
    #cellmakers_thi <- as.vector(unlist(dtrain_thi_rank[,1]))
    #cellmakers_four <- as.vector(unlist(dtrain_four_rank[,1]))
    celltype_common <- intersect(cellmakers_fir,cellmakers_sec)
    #celltype_common <- intersect(celltype_common,cellmakers_thi)
    #celltype_common <- intersect(celltype_common,cellmakers_four)
    cellmarkers_common[[i]] <- celltype_common
    celltype_all <- union(cellmakers_fir,cellmakers_sec)
    #celltype_all <- union(celltype_all,cellmakers_thi)
    #celltype_all <- union(celltype_all,cellmakers_four)
    cellmarkers_all[[i]] <- celltype_all
}
outputfile <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_allmarkers/COPD_Lung.RData"
save(celltypes,cellmarkers_common,cellmarkers_all,file=outputfile)
cat("The identification of cellmarker for all celltypes is completed!")


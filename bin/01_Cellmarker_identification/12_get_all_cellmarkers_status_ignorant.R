library(data.table)

filepath <- "/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_Breast_normal_SPred_final_status_ignorant/"
files <- list.files(filepath)
celltypes_Dep <- NULL
cellmarkers_all <- vector(mode = "list", length = length(files))
for(i in 1:length(files)){
    filename <- files[i]
    celltype <- gsub("BRCA_cellmarkers_","",filename)
    celltype <- gsub(".RData","",celltype)
    celltypes_Dep <- c(celltypes_Dep,celltype)
    filename_full <- paste0(filepath,filename)
    load(filename_full)
    cellmarkers <- as.vector(unlist(dtrain_rank[,1]))
    cellmarkers_all[[i]] <- cellmarkers
}
outputfile <- "/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_allmarkers/Breast_normal_SPred_final_all_status_ignorant.RData"
save(celltypes_Dep,cellmarkers_all,file=outputfile)

library(data.table)

filepath <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_dep_coarse_5Bs_status_ignorant/"
files <- list.files(filepath)
celltypes_Dep <- NULL
cellmarkers_all <- vector(mode = "list", length = length(files))
for(i in 1:length(files)){
    filename <- files[i]
    celltype <- gsub("Dep_cellmarkers_","",filename)
    celltype <- gsub(".RData","",celltype)
    celltypes_Dep <- c(celltypes_Dep,celltype)
    filename_full <- paste0(filepath,filename)
    load(filename_full)
    cellmarkers <- as.vector(unlist(dtrain_rank[,1]))
    cellmarkers_all[[i]] <- cellmarkers
}
outputfile <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_allmarkers/Dep_all_coarse_celltypes_5Bs_status_ignorant.RData"
save(celltypes_Dep,cellmarkers_all,file=outputfile)

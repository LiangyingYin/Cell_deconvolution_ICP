library(data.table)

filepath <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_dep_coarse_5Bs/"
files <- list.files(filepath)
celltypes_Dep <- NULL
cellmarkers_common <- vector(mode = "list", length = length(files))
cellmarkers_all <- vector(mode = "list", length = length(files))
for(i in 1:length(files)){
    filename <- files[i]
    celltype <- gsub("Dep_cellmarkers_","",filename)
    celltype <- gsub(".RData","",celltype)
    celltypes_Dep <- c(celltypes_Dep,celltype)
    filename_full <- paste0(filepath,filename)
    load(filename_full)
    cellmakers_fir <- as.vector(unlist(dtrain_fir_rank[,1]))
    cellmakers_sec <- as.vector(unlist(dtrain_sec_rank[,1]))
    
    celltype_common <- intersect(cellmakers_fir,cellmakers_sec)
    cellmarkers_common[[i]] <- celltype_common
    celltype_all <- union(cellmakers_fir,cellmakers_sec)
    cellmarkers_all[[i]] <- celltype_all
}
outputfile <- "/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_allmarkers/Dep_all_coarse_celltypes_5Bs.RData"
save(celltypes_Dep,cellmarkers_common,cellmarkers_all,file=outputfile)
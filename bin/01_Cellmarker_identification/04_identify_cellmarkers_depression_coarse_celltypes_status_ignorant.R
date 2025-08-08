library(data.table)
library(xgboost)
library(stringr)

celltypes_path <- '/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/Depression/Coarse_Celltypes/'
celltypes_files <- list.files(celltypes_path)
dep <- as.matrix(fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/Depression.csv",header=TRUE))
#dep <- t(dep)
cellnames <- fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/Depression/GSE144136_barcodes.csv",header=TRUE)
cellnames_all <- as.vector(unlist(cellnames[,2]))
cellnames_all_split <- str_split_fixed(cellnames_all,"\\.",2)
cell_info <- str_split_fixed(cellnames_all_split[,2],"_",4)
cellinfo <- cbind(cellnames_all,cellnames_all_split[,1],cell_info)
colnames(cellinfo) <- c("cellnames_all","celltype","patient_ID","status","batch","UMI")
selected_index <- which(cellinfo[,5]!="B6") # leave B6 for validation
# remove cells from B6 
dep <- dep[selected_index,]
# remove cells from B6
cellnames_all <- cellnames_all[selected_index]

genes <- fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/Depression/GSE144136_genes.csv",header=TRUE)
genes_names <- as.vector(unlist(genes[,2]))
celltypes <- NULL
colnames(dep) <- genes_names
cell_label <- rep(0,nrow(dep))
#dep_all <- cbind(dep,cell_label)
#len <- ncol(dep_all)
#rownames(dep_all) <- cellnames_all
len <- ncol(dep)
rownames(dep) <- cellnames_all
dep <- as.matrix(dep)
class(dep) <- 'numeric'

#for(i in 1:length(celltypes_files)){
    i <- 5
    file_name <- celltypes_files[i]
    subtype_cells <- fread(paste0(celltypes_path,file_name),header=TRUE)
    subtype_cellnames <- as.vector(unlist(subtype_cells[,2]))
    cell_label[rownames(dep)%in%subtype_cellnames] <- 1
    # ctrl_cell_label <- cell_label[ctrl_index]
    # case_cell_label <- cell_label[case_index]
    celltype <- gsub("GSE144136_","",file_name)
    celltype <- gsub("_filtered_cells.csv","",celltype)
    celltypes <- c(celltypes,celltype)
    cat("This is the identification of cell marker for",celltype,"\n")
    # cat("There are",length(ctrl_cell_label),celltype,"cells in the 1st dataset. \n")
    
    dep <- data.matrix(dep, rownames.force = NA)
    dtrain <- xgb.DMatrix(data = dep, label = cell_label)
    dtrain_model <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_rank <- xgb.importance(model = dtrain_model)
    #cellmarker_fir_list[[i]] <- dtrain_rank
    outputfile <- paste0("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_dep_coarse_5Bs_status_ignorant/Dep_cellmarkers_",celltype,".RData")
    save(celltype,dtrain_rank,file=outputfile)
    
#}
# outputfile <- paste0("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/Depression_cellmarkers.RData")
# save(Celltypes,dtrain_fir_rank,dtrain_sec_rank,dtrain_thi_rank,file=outputfile)
# outputfile <- paste0("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_dep_coarse_5Bs_status_ignorant/Dep_cellmarkers_",celltype,".RData")
# save(celltype,dtrain_rank,file=outputfile)
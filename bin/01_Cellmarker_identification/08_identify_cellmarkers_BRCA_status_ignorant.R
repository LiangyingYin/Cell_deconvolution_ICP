library(data.table)
library(xgboost)
library(stringr)

normal <- fread("/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/data/BRCA_normal.csv")
cell_info <- fread("/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/data/BRCA_normal_cellinfo.csv")
unique(cell_info$broad_cell_type)
celltypes_normal <- as.vector(unlist(cell_info$broad_cell_type))

common_genes <- fread("/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/data/BRCA_common_genes.csv",header=TRUE)
common_genes <- as.vector(unlist(common_genes))
# combined cell types, the left ones are from normal breast mammary tissue, the right ones are from the cancer
# (lumhr, lumsec, basal)-> Epithelial <-(Normal Epithelial,Cancer Epithelial)
# bcells -> B-cells <-B-cells, tcells -> T-cells <- T-cells, myloid -> Myeloid <- Myeloid
# fibroblasts -> fibroblasts <-CAFs, (lymphatic,vascular) -> Endothelial <- Endothelial
# pericytes ->PVL <-PVL
celltypes_normal[celltypes_normal=="lumhr"] <- "Epithelial"
celltypes_normal[celltypes_normal=="lumsec"] <- "Epithelial"
celltypes_normal[celltypes_normal=="basal"] <- "Epithelial"
celltypes_normal[celltypes_normal=="bcells"] <- "B-cells"
celltypes_normal[celltypes_normal=="tcells"] <- "T-cells"
celltypes_normal[celltypes_normal=="myeloid"] <- "Myeloid"
celltypes_normal[celltypes_normal=="lymphatic"] <- "Endothelial"
celltypes_normal[celltypes_normal=="vascular"] <- "Endothelial"
celltypes_normal[celltypes_normal=="pericytes"] <- "PVL"
unique(celltypes_normal)

# normal_genes <- colnames(normal)
# common_genes <- intersect(genes_names,normal_genes)
# length(common_genes)
# fwrite(as.data.frame(common_genes),"/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/BRCA_common_genes.csv",sep='\t')
# case_index <- which(cell_info$disease=="breast cancer") # "P113" "P114" "P115" "P116" "P117" "P119" "P120" "P122" "P124"
# cell_info_case <- cell_info[case_index,]
# unique(cell_info_case$donor_id)

# ctrl_index <- which(cell_info$disease=="normal")
# cell_info_ctrl <- cell_info[ctrl_index,]
# unique(cell_info_ctrl$donor_id)

#
external_ids_normal <- as.vector(unlist(unique(cell_info$donor_id)[101:126]))
normal_selected_index <- which(!(cell_info$donor_id%in%external_ids_normal))
normal_selected <- normal[normal_selected_index,]
celltypes_normal <- celltypes_normal[normal_selected_index]
cell_info_selected <- cell_info[normal_selected_index,]

# case_index <- which(cell_info_selected$disease=="breast cancer")
# case <- normal_selected[case_index,]
# label_fir <- celltypes_normal[case_index]
# case <- as.matrix(case)
# class(case) <- "numeric"
# unique(label_fir)


# ctrl_index <- which(cell_info_selected$disease=="normal")
# ctrl <- normal_selected[ctrl_index,]
# label_sec <- celltypes_normal[ctrl_index]
# #normal_selected <- subset(normal_selected,select=common_genes)
# # normal_selected <- as.matrix(normal_selected)
# # class(normal_selected) <- "numeric"
# ctrl <- as.matrix(ctrl)
# class(ctrl) <- "numeric"
# unique(label_sec)


normal_selected <- as.matrix(normal_selected)
class(normal_selected) <- "numeric"
label <- celltypes_normal
celltypes <- unique(celltypes_normal)

#for(i in 1:length(celltypes_normal)){
    i <- 7
    celltype <- celltypes[i]
    #celltype <- celltypes_normal[i]
    label_binary <- rep(0,length(label))
    cat("This is the identification of cell marker for",celltype,"\n")
    celltype_index <- which(label==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_binary[celltype_index] <- 1
    dtrain <- xgb.DMatrix(data = normal_selected, label = label_binary)
    dtrain_model <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_rank <- xgb.importance(model = dtrain_model)
    cat("There are",nrow(dtrain_rank),"identified cellmakers in the 1st dataset.\n")
    #cellmarker_fir_list[[i]] <- dtrain_fir_rank
    # train xgboost for the 2nd dataset
    # label_sec_binary <- rep(0,length(label_sec))
    # celltype_index <- which(label_sec==celltype)
    # cat("There are",length(celltype_index),celltype,"cells in the 2nd dataset. \n")
    # label_sec_binary[celltype_index] <- 1
    # dtrain_sec <- xgb.DMatrix(data = ctrl, label = label_sec_binary)
    # dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    # dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
    # cat("There are",nrow(dtrain_sec_rank),"identified cellmakers in the 1st dataset.\n")
     
    outputfile <- paste0("/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_Breast_normal_status_ignorant/BRCA_cellmarkers_",celltype,".RData")
    save(celltype,dtrain_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")
   
#}

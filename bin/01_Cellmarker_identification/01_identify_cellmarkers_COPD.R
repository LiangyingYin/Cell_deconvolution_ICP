library(data.table)
library(xgboost)
library(stringr)

scRNA <- fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/COPD_Lung.csv")
cell_info <- fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/COPD_Lung_cellinfo.csv")

# The following is only for test
colnames(cell_info)
outcome <- cell_info$Disease_Identity
celltype <- cell_info$CellType_Category
unique(celltype)
# celltype_minor <- cell_info$Subclass_Cell_Identity
# unique(celltype_minor)

donor_ids <-cell_info$Subject_Identity
scRNA_all <- cbind(scRNA,celltype,donor_ids,outcome)

# preserve one ctrl and one case for external validation
external_ids <- c("001C","002C","003C","222C","034C","065C","081C","084C","052CO","056CO","137CO") # The first 7 are ctrl ids, the remaining are case ids

# external_index <- which(scRNA_all$donor_id%in%external_ids) # length(selected_index) 14530
# external <- scRNA_all[external_index,] # dim 26645 45950
# unique(external[,"celltype"])

selected_index <- which(!scRNA_all$donor_id%in%external_ids) # length(selected_index) 14530
scRNA_selected <- scRNA_all[selected_index,] # 139114  45950
rm(scRNA_all)

ctrl_index <- which(scRNA_selected$outcome=="Control")
ctrl <- scRNA_selected[ctrl_index,]
# unique(ctrl[,"celltype"])
# unique(ctrl[,"donor_ids"]) # overall 28 ctrls, 22 remained for training

case_index <- which(scRNA_selected$outcome=="COPD")
case <- scRNA_selected[case_index,]
# unique(case[,"celltype"])
# unique(case[,"donor_ids"]) # overall 18 cases ,15 remained for training

feat_len <- ncol(scRNA_selected)- 3
rm(scRNA_selected)

data_fir <- subset(ctrl,select=(1:feat_len))
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- as.vector(unlist(ctrl$celltype))

data_sec <- subset(case,select=(1:feat_len))
data_sec <- as.matrix(data_sec)
class(data_sec) <- "numeric"
label_sec <- as.vector(unlist(case$celltype))

celltypes <- unique(label_fir)
celltypes

for(i in 1:length(celltypes)){
#for(i in 9:9){
    #i <- 1
    celltype <- celltypes[i]
    label_fir_binary <- rep(0,length(label_fir))
    cat("This is the identification of cell marker for",celltype,"\n")
    celltype_index <- which(label_fir==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_fir_binary[celltype_index] <- 1
    dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
    dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
    cat("There are",nrow(dtrain_fir_rank),"identified cellmakers in the 1st dataset.\n")
    #cellmarker_fir_list[[i]] <- dtrain_fir_rank
    # train xgboost for the 2nd dataset
    label_sec_binary <- rep(0,length(label_sec))
    celltype_index <- which(label_sec==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 2nd dataset. \n")
    label_sec_binary[celltype_index] <- 1
    dtrain_sec <- xgb.DMatrix(data = data_sec, label = label_sec_binary)
    dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
    cat("There are",nrow(dtrain_sec_rank),"identified cellmakers in the 2nd dataset.\n")
    #cellmarker_sec_list[[i]] <- dtrain_sec_rank
    
    outputfile <- paste0("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/res/res_COPD/COPD_cellmarkers_",celltype,".RData")
    save(celltype,dtrain_fir_rank,dtrain_sec_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")
}

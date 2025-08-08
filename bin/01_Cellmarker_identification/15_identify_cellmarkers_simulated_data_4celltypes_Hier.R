library(data.table)
library(xgboost)
library(stringr)
library(SingleCellExperiment)
library(scuttle)

load("/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/res/14_addon_simulation/4celltypes/toydata_batch_20genes_16.RData")
normal <- NULL
normal_celltype <- NULL
normal_batch <- NULL
patient <- NULL
patient_celltype <- NULL
patient_batch <- NULL
# We only take the first 4 patients and controls for analysis
for(i in 1:7){
    subject_normal <- as.matrix(logcounts(normal_list[[i]]))
    subject_normal_meta <- colData(normal_list[[i]])
    subject_normal_genes <- rownames(normal_list[[i]])
    subject_normal_final <- t(subject_normal)
    subject_normal_celltype <- as.vector(unlist(subject_normal_meta$cell_type))
    normal_celltype <- c(normal_celltype,subject_normal_celltype)
    subject_normal_batch <- as.vector(unlist(subject_normal_meta$batch))
    normal_batch <- c(normal_batch,subject_normal_batch)
    normal <- rbind(normal,subject_normal_final)

    subject_patient <- as.matrix(logcounts(patient_list[[i]]))
    subject_patient_genes <- rownames(patient_list[[i]])
    subject_patient_meta <- colData(patient_list[[i]])
    subject_patient_celltype <- as.vector(unlist(subject_patient_meta$cell_type))
    subject_patient_batch <- as.vector(unlist(subject_patient_meta$batch))
    patient_celltype <- c(patient_celltype,subject_patient_celltype)
    patient_batch <- c(patient_batch,subject_patient_batch)
    subject_patient_final <- t(subject_patient)
    patient <- rbind(patient,subject_patient_final)

}

colnames(normal) <- subject_normal_genes
colnames(patient) <- subject_patient_genes
scRNA <- rbind(normal,patient)
colnames(scRNA) <- colnames(normal)
celltype <- c(normal_celltype,patient_celltype)
batch <- c(normal_batch,patient_batch) # unique(batch): CTRL,STIM
# scRNA_all <- cbind(scRNA,as.vector(unlist(celltype)),as.vector(unlist(batch)))
scRNA_all <- cbind(scRNA,celltype,batch)
colnames(scRNA_all) <- c(colnames(normal),"celltype","batch")
scRNA_all <- as.data.frame(scRNA_all)

feat_len <- ncol(scRNA_all) - 2
cat("There are",feat_len,"genes in the simulated scRNA-seq dataset.\n")
data <- subset(scRNA_all,select=(1:feat_len))
data <- as.matrix(data)
class(data) <- "numeric"
label <- as.vector(unlist(scRNA_all$celltype))


celltypes <- unique(label)

for(i in 1:length(celltypes)){
    # train xgboost for the control dataset
    celltype <- celltypes[i]
    cat("This is the identification of cell marker for",celltype,"\n")
    label_binary <- rep(0,length(label))
    celltype_index <- which(label==celltype)
    cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
    label_binary[celltype_index] <- 1
    dtrain <- xgb.DMatrix(data = data, label = label_binary)
    dtrain_model <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
    dtrain_rank <- xgb.importance(model = dtrain_model)
    cat("There are",nrow(dtrain_rank),"identified cellmarkers in the 1st dataset.\n")

    outputfile <- paste0("/home/yinly/Cell_deconvolution/00_Simulation/res_more_confounding_genes/res_simulated_status_ingorant_4celltypes/pbmc_simulated_",celltype,".RData")
    save(celltype,dtrain_rank,file=outputfile)
    cat("The identification of cellmarkers is completed!\n\n")  
}




# library(data.table)
# library(xgboost)
# library(stringr)
# library(SingleCellExperiment)
# library(scuttle)

# load("/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/res/14_addon_simulation/toydata_2.RData")
# normal <- NULL
# normal_celltype <- NULL
# normal_batch <- NULL
# patient <- NULL
# patient_celltype <- NULL
# patient_batch <- NULL
# # We only take the first 4 patients and controls for analysis
# for(i in 1:4){
#     subject_normal <- as.matrix(logcounts(normal_list[[i]]))
#     subject_normal_meta <- colData(normal_list[[i]])
#     subject_normal_genes <- rownames(normal_list[[i]])
#     subject_normal_final <- t(subject_normal)
#     subject_normal_celltype <- as.vector(unlist(subject_normal_meta$cell_type))
#     normal_celltype <- c(normal_celltype,subject_normal_celltype)
#     subject_normal_batch <- as.vector(unlist(subject_normal_meta$batch))
#     normal_batch <- c(normal_batch,subject_normal_batch)
#     normal <- rbind(normal,subject_normal_final)

#     subject_patient <- as.matrix(logcounts(patient_list[[i]]))
#     subject_patient_genes <- rownames(patient_list[[i]])
#     subject_patient_meta <- colData(patient_list[[i]])
#     subject_patient_celltype <- as.vector(unlist(subject_patient_meta$cell_type))
#     subject_patient_batch <- as.vector(unlist(subject_patient_meta$batch))
#     patient_celltype <- c(patient_celltype,subject_patient_celltype)
#     patient_batch <- c(patient_batch,subject_patient_batch)
#     subject_patient_final <- t(subject_patient)
#     patient <- rbind(patient,subject_patient_final)

# }

# colnames(normal) <- subject_normal_genes
# colnames(patient) <- subject_patient_genes
# scRNA <- rbind(normal,patient)
# colnames(scRNA) <- colnames(normal)
# celltype <- c(normal_celltype,patient_celltype)
# batch <- c(normal_batch,patient_batch) # unique(batch): CTRL,STIM
# # scRNA_all <- cbind(scRNA,as.vector(unlist(celltype)),as.vector(unlist(batch)))
# scRNA_all <- cbind(scRNA,celltype,batch)
# colnames(scRNA_all) <- c(colnames(normal),"celltype","batch")
# scRNA_all <- as.data.frame(scRNA_all)

# feat_len <- ncol(scRNA_all) - 2
# cat("There are",feat_len,"genes in the simulated scRNA-seq dataset.\n")
# data <- subset(scRNA_all,select=(1:feat_len))
# data <- as.matrix(data)
# class(data) <- "numeric"
# label <- as.vector(unlist(scRNA_all$celltype))


# celltypes <- unique(label)

# for(i in 1:length(celltypes)){
#     # train xgboost for the control dataset
#     celltype <- celltypes[i]
#     cat("This is the identification of cell marker for",celltype,"\n")
#     label_binary <- rep(0,length(label))
#     celltype_index <- which(label==celltype)
#     cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
#     label_binary[celltype_index] <- 1
#     dtrain <- xgb.DMatrix(data = data, label = label_binary)
#     dtrain_model <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
#     dtrain_rank <- xgb.importance(model = dtrain_model)
#     cat("There are",nrow(dtrain_rank),"identified cellmarkers in the 1st dataset.\n")

#     outputfile <- paste0("/home/yinly/Cell_deconvolution/00_Simulation/res_simulated_status_ignorant/pbmc_simulated_",celltype,".RData")
#     save(celltype,dtrain_rank,file=outputfile)
#     cat("The identification of cellmarkers is completed!\n\n")  
# }


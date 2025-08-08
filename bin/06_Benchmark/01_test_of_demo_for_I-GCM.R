library(data.table)
library(xgboost)
# library(stringr)
library(GeneralisedCovarianceMeasure)
library(pcalg)
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/PCSelect_Parallel_Update.R")
source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/I-GCM.R")

# Load toy data
load("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/res/14_addon_simulation/2celltypes/toydata_batch.RData")
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
scRNA_all$status <- c(rep("CTRL",nrow(normal)),rep("STIM",nrow(patient)))
fwrite(scRNA_all,"/mnt/home/yinly/projects/Cell_deconvolution/dat/Github/Toydata_2celltypes.txt",sep="\t")

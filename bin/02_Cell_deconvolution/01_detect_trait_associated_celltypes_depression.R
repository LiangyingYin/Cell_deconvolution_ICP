library(data.table)
library(biglm)

covariates <- fread("/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_31March2022.txt")
IID <- covariates$FID
Dep <- covariates$Depression

na_index <- which(is.na(Dep))
Dep[na_index] <- 0
Dep_all <- cbind(IID,Dep)
colnames(Dep_all) <- c("IID","Depression")

ukbb_cell_common <- fread("/exeh_3/yinly/Cell_deconvolution/02_Cell_deconvolution/res/UKBB/ukbb_depression_cortex_common.tsv")
cell_common_dat <- merge(ukbb_cell_common,Dep_all,by="IID",sort=FALSE)

ukbb_cell_all <- fread("/exeh_3/yinly/Cell_deconvolution/02_Cell_deconvolution/res/UKBB/ukbb_depression_cortex_all.tsv")
cell_all_dat <- merge(ukbb_cell_all,Dep_all,by="IID",sort=FALSE)
ukbb_cell_ignorant <- fread("/exeh_3/yinly/Cell_deconvolution/02_Cell_deconvolution/res/UKBB/ukbb_depression_cortex_all_status_ignorant.tsv")
cell_ignorant_dat <- merge(ukbb_cell_ignorant,Dep_all,by="IID",sort=FALSE)

feat_len <- ncol(cell_all_dat)-ncol(Dep_all) + 1
res <- matrix(0,nrow=(feat_len-1),ncol=7)
for(i in 2:feat_len){
    celltype <- colnames(cell_common_dat)[i]
    #res[(i-1),1] <- celltype
    cat("This is the analysis for",celltype,"\n")
    # common
    common_cell_prop <- unlist(subset(cell_common_dat,select=i))
    common_obj <- biglm(as.numeric(cell_common_dat$Depression) ~ common_cell_prop,data=cell_common_dat)
    common_res <- summary(common_obj)$mat[2,c(1,5)]
    #common_p_adjust <-p.adjust(common_res[,2],length=nrow())
    # all
    all_cell_prop <- unlist(subset(cell_all_dat,select=i))
    all_obj <- biglm(as.numeric(cell_all_dat$Depression) ~ all_cell_prop,data=cell_all_dat)
    all_res <- summary(all_obj)$mat[2,c(1,5)]
    # ignorant
    ignorant_cell_prop <- unlist(subset(cell_ignorant_dat,select=i))
    ignorant_obj <- biglm(as.numeric(cell_ignorant_dat$Depression) ~ ignorant_cell_prop,data=cell_ignorant_dat)
    ignorant_res <- summary(ignorant_obj)$mat[2,c(1,5)]
    assoc <- c(celltype,common_res,all_res,ignorant_res)
    res[(i-1),] <- assoc
}
colnames(res) <- c("celltype","beta_common","p_common","beta_all","p_all","beta_ignorant","p_ignorant")
common_p_adjust <-p.adjust(res[,3],method ="BH",n=nrow(res))
all_p_adjust <- p.adjust(res[,5],method ="BH",n=nrow(res))
ignorant_p_adjust <- p.adjust(res[,7],method ="BH",n=nrow(res))
res_adjust <- cbind(res,common_p_adjust,all_p_adjust,ignorant_p_adjust)
fwrite(as.data.frame(res_adjust),"/exeh_3/yinly/Cell_deconvolution/02_Cell_deconvolution/res/res_disease_status_assoc/ukbb_depression_frontal_cortex_assoc.txt",sep="\t")


import numpy as np
import scanpy as sc
import pandas as pd
import scipy as sp

def normalize(adata, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):

    if "n_counts" not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata.obs['n_counts'] = adata.obs['total_counts']

    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)

    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata

    if size_factors:
        #sc.pp.normalize_per_cell(adata)
        sc.pp.normalize_total(adata,target_sum=1e6) # get normalized TPM data
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0

    if logtrans_input:
        sc.pp.log1p(adata)

    if normalize_input:
        sc.pp.scale(adata)

    return adata


sample_meta = pd.read_excel('/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/data/GSE158055_sample_metadata.xlsx')
sample_id_selected = sample_meta.iloc[20:304,0].tolist() # get all data
#sample_id_selected = sample_meta.iloc[55:103,0].tolist() # get all data
filepath_prefix = '/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/data/02_patients_splitting/sampleid_'
filepath_suffix = '.h5ad'
sample_id_selected_files = [filepath_prefix + sample_id + filepath_suffix for sample_id in sample_id_selected]
scRNA_list = []
for sample_id_file in sample_id_selected_files:
    scRNA = sc.read_h5ad(sample_id_file)
    scRNA_list.append(scRNA)
    #scRNA_exp_origi = scRNA.X
    print(scRNA.X.shape)
    # scRNA_exp = pd.DataFrame(scRNA.X)
    # scRNA_obs = scRNA.obs 
scRNA_merge = scRNA_list[0]
for adata in scRNA_list[1:]:
    scRNA_merge = scRNA_merge.concatenate(adata)

sc.pp.calculate_qc_metrics(scRNA_merge, inplace=True)
scRNA_merge_norm = normalize(scRNA_merge)
scRNA_merge_exp = pd.DataFrame(scRNA_merge_norm.X)
var_names = scRNA_merge_norm.var.index.tolist()
scRNA_merge_exp.columns = var_names
scRNA_merge_obs = scRNA_merge_norm.obs
row_len = len(scRNA_merge_obs.index.tolist())
# row_indices = list(range(1, (row_len+1)))
row_indices = list(range(0, row_len))
scRNA_merge_obs['ID'] = row_indices
cell_info = scRNA_merge_obs[['orig.ident','patient_ident','batch','ID']]
cell_info = pd.DataFrame(cell_info)
IDs = cell_info['ID'].to_numpy()
cell_info = cell_info.set_index(IDs)
# merge scRNA_exp for Celltypecytes with Celltype_info, notably, the row index indicate the same cell in the two dataframe
scRNA_merge_exp_all = pd.concat([scRNA_merge_exp, cell_info], axis=1) # Celltype_expr(31649, 23878),Celltype_info(31649, 3)
print(scRNA_merge_exp_all.shape)
scRNA_merge_exp_all.to_csv('/exeh_3/yinly/Cell_deconvolution/01_Cellmarker_identification/res/COVID_all.csv',sep='\t', index=False)

import numpy as np
import scanpy as sc
import pandas as pd
import scipy as sp

def normalize(adata, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):

    if "n_counts" not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata.obs['n_counts'] = adata.obs['total_counts']

    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)

    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata

    if size_factors:
        #sc.pp.normalize_per_cell(adata)
        sc.pp.normalize_total(adata,target_sum=1e6) # get normalized TPM data
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0

    if logtrans_input:
        sc.pp.log1p(adata)

    if normalize_input:
        sc.pp.scale(adata)

    return adata


sample_meta = pd.read_excel('/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/GSE158055_sample_metadata.xlsx')
#sample_id_selected = sample_meta.iloc[20:165,0].tolist() # get all data
sample_id_selected = sample_meta.iloc[55:103,0].tolist() # get all data
filepath_prefix = '/home/yujia/project/2023-02-04-Drug_repurposing/cancercn/dat/01_dataset_COVID_single_cell_GSE158055/02_preprocess/02_patients_splitting/sampleid_'
filepath_suffix = '.h5ad'
sample_id_selected_files = [filepath_prefix + sample_id + filepath_suffix for sample_id in sample_id_selected]
scRNA_list = []
for sample_id_file in sample_id_selected_files:
    scRNA = sc.read_h5ad(sample_id_file)
    scRNA_list.append(scRNA)
    #scRNA_exp_origi = scRNA.X
    print(scRNA.X.shape)
    # scRNA_exp = pd.DataFrame(scRNA.X)
    # scRNA_obs = scRNA.obs 
scRNA_merge = scRNA_list[0]
for adata in scRNA_list[1:]:
    scRNA_merge = scRNA_merge.concatenate(adata)

sc.pp.calculate_qc_metrics(scRNA_merge, inplace=True)
scRNA_merge_norm = normalize(scRNA_merge)
scRNA_merge_exp = pd.DataFrame(scRNA_merge_norm.X)
var_names = scRNA_merge_norm.var.index.tolist()
scRNA_merge_exp.columns = var_names
scRNA_merge_obs = scRNA_merge_norm.obs
row_len = len(scRNA_merge_obs.index.tolist())
# row_indices = list(range(1, (row_len+1)))
row_indices = list(range(0, row_len))
scRNA_merge_obs['ID'] = row_indices
cell_info = scRNA_merge_obs[['orig.ident','patient_ident','batch','ID']]
cell_info = pd.DataFrame(cell_info)
IDs = cell_info['ID'].to_numpy()
cell_info = cell_info.set_index(IDs)
# merge scRNA_exp for Celltypecytes with Celltype_info, notably, the row index indicate the same cell in the two dataframe
scRNA_merge_exp_all = pd.concat([scRNA_merge_exp, cell_info], axis=1) # Celltype_expr(31649, 23878),Celltype_info(31649, 3)
print(scRNA_merge_exp_all.shape)
scRNA_merge_exp_all.to_csv('/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/COVID_Wuhan.csv',sep='\t', index=False)




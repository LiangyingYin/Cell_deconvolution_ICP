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

adata = sc.read_h5ad("/mnt/home/yujia/project/2023-04-24-Cell-decomposition/cell_decomposition/dat/09_pbmc_data/01_preprocess_h5ad/5810k_5/5810k_5.h5ad")
sc.pp.calculate_qc_metrics(adata, inplace=True)
scRNA_norm = normalize(adata)
scRNA_exp = pd.DataFrame(scRNA_norm.X)
var_names = scRNA_norm.var.index.tolist()
scRNA_exp.columns = var_names
scRNA_obs = scRNA_norm.obs
row_len = len(scRNA_obs.index.tolist())

scRNA_exp.to_csv('/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/PBMC/PBMC_scRNA_exp_5810k_5.csv',sep='\t', index=True)
scRNA_obs.to_csv('/mnt/home/yinly/projects/Cell_deconvolution/dat/01_Cellmarker_identification/PBMC/PBMC_scRNA_meta_5810k_5.csv',sep='\t', index=True)
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import gc
sc.settings.verbosity = 1

INPUT_H5AD   = "./merge_total.h5ad"
META_CSV     = "./meta_total.csv"
OUTPUT_H5AD  = "./BDdata_processed.h5ad"
MT_PREFIX    = ("MT-", "mt-")
MIN_GENES    = 500
MAX_PCT_MT   = 0.50
MIN_COUNTS   = 1000
N_TOP_HVG    = 2000
N_PCS        = 50
N_NEIGHBORS  = 15
BATCH_KEY    = "patient"
RESOLUTION   = 1.0
RANDOM_STATE = 0

adata = sc.read_h5ad(INPUT_H5AD)
adata.var_names_make_unique()
print(adata)
meta = pd.read_csv(META_CSV, index_col=1)
adata.obs = adata.obs.join(meta, how="left", rsuffix="_meta")
adata.var["mt"] = adata.var_names.str.startswith(MT_PREFIX)
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pl.violin(
    adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0, multi_panel=True,
)
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
adata = adata[
    (adata.obs["pct_counts_mt"] < MAX_PCT_MT)& 
    (adata.obs["total_counts"] > MIN_COUNTS), :
].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(
    adata, n_top_genes=N_TOP_HVG, flavor="seurat", inplace=True
)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS,
                    random_state=RANDOM_STATE)
sc.tl.umap(adata, random_state=RANDOM_STATE)
sc.tl.leiden(adata, resolution=RESOLUTION, random_state=RANDOM_STATE,
             key_added="leiden")
sc.pl.umap(adata, color=["leiden"], legend_loc="on_data", legend_fontsize=8)

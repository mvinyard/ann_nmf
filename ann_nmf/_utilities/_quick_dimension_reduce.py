import scanpy as sc

def _quick_dimension_reduce(adata, n_pcs=50):
    
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

import scanpy as sc

def _quick_norm_log_raw(adata, counts_per_cell_after=10000, return_adata=False, silent=False):
    
    """
    Quickly library-size- and log-normalize raw counts. Saves raw adata. 
    
    Parameters:
    -----------
    adata
        AnnData object containing raw counts of cells x genes in adata.X
        type: anndata._core.anndata.AnnData
        
    counts_per_cell_after
        Passed to sc.pp.normalize_per_cell as the target number of total counts per cell.
        default: 10_000
        type: int
    
    return_adata
        Return the updated AnnData object. 
        default: False
        type: bool
        
    silent
        Run without printing updated AnnData. 
        default: False
        type: bool
        
    Returns:
    --------
    None [adata]
        By default, adata is modified in-place. 
        
    Notes:
    ------
    Quick implementation. May not be ideal for everyone. 
    """

    adata.raw = adata
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=counts_per_cell_after)
    sc.pp.log1p(adata)
    if not silent:
        print(adata)

    if return_adata:
        adata
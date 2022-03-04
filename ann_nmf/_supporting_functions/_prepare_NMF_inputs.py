
import pandas as pd
import scanpy as sc

from ._MessageModule import _MessageModule


def _fetch_highly_variable_genes(
    adata, hv_key="highly_variable", min_mean=0.0125, max_mean=3, min_disp=0.5, **kwargs
):

    """"""
    
    msg = _MessageModule()
    
    if not hv_key in adata.var.columns.tolist():
        msg.scanpy_hv_genes(hv_key)
        sc.pp.highly_variable_genes(
            adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, **kwargs
        )

    print("\n\t**Using {} highly variable genes.".format(sum(adata.var.highly_variable)))
    genes_for_nmf = set(adata.var[adata.var["highly_variable"]].index)

    return genes_for_nmf


def _filter_mito(genes_for_nmf):

    mito_genes = {x for x in genes_for_nmf if x.startswith("MT-")}
    print("\t**Filtering {} mito genes.".format(len(mito_genes)))
    genes_for_nmf -= mito_genes


def _filter_ribo(genes_for_nmf):
    ribo_genes = {
        x for x in genes_for_nmf if x.startswith("RPS") or x.startswith("RPL")
    }
    print("\t**Filtering {} ribo genes.".format(len(ribo_genes)))
    genes_for_nmf -= ribo_genes


def _prepare_NMF_input_genes(
    adata, use_highly_variable=True, filter_mito=True, filter_ribo=True
):

    """"""

    if use_highly_variable:
        genes_for_nmf = _fetch_highly_variable_genes(adata)
    else:
        genes_for_nmf = set(adata.var.index)

    if filter_mito:
        _filter_mito(genes_for_nmf)
    if filter_ribo:
        _filter_ribo(genes_for_nmf)

    return genes_for_nmf

def _fetch_raw_inputs(adata):

    """"""

    X_ = adata.raw.X.toarray()
    idx_ = adata.raw.var_names
    cols_ = adata.raw.obs_names

    return X_, idx_, cols_


def _prepare_NMF_input_matrix(adata, use_raw=True, layer=False, use_X=False):

    """"""

    if use_raw:
        X_, idx_, cols_ = _fetch_raw_inputs(adata)
        return pd.DataFrame(data=X_.T, index=idx_, columns=cols_)
    elif layer:
        assert layer in adata.layers, "Please save input in adata.layers['{}']".format(
            layer
        )
        return pd.DataFrame(
            data=adata.layers[layer].toarray().T,
            index=adata.var_names,
            columns=adata.obs_names,
        )

    elif use_X:
        return pd.DataFrame(
            data=adata.X.toarray().T, index=adata.var_names, columns=adata.obs_names
        )

    else:
        print("Input not properly formatted. Defaulting to adata.X")
        return pd.DataFrame(
            data=adata.X.toarray().T, index=adata.var_names, columns=adata.obs_names
        )
    
def _prepare_NMF_inputs(
    adata,
    use_raw=True,
    layer=False,
    use_X=False,
    use_highly_variable=True,
    filter_mito=True,
    filter_ribo=True,
    **kwargs
):

    """"""

    mtx_df = _prepare_NMF_input_matrix(adata, use_raw=use_raw, layer=layer, use_X=use_X)
    genes_for_nmf = _prepare_NMF_input_genes(
        adata,
        use_highly_variable=use_highly_variable,
        filter_mito=filter_mito,
        filter_ribo=filter_ribo,
        **kwargs
    )

    return mtx_df.loc[genes_for_nmf], genes_for_nmf
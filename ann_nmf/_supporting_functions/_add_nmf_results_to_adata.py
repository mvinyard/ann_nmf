
import os
import pandas as pd
from ._get_best_result import _get_best_result
import signatureanalyzer as sa

def _uns_dict_to_df(adata, key):

    """"""
    return pd.DataFrame(
        data=adata.uns[key]["values"],
        columns=adata.uns[key]["cols"],
        index=adata.uns[key]["rows"],
    )

def _add_consensus_clusters(adata, nmf_result):

    adata.obs = adata.obs.join(_nmf_result["consensus"])
    adata.obs["clusters"] = adata.obs["clusters"].astype("category")
    adata.obs = adata.obs.rename(columns={"clusters": "consensus_cluster"})


def _add_H_to_adata(adata, H):

    adata.obs = adata.obs.join(H)
    adata.obs["max_id"] = adata.obs["max_id"].astype("category")

    adata.uns["signatures"] = list(H)[:-3]
    adata.obsm["X_nmf"] = H.iloc[:, :-3].values


def _add_signatures_to_adata(adata, signatures):

    adata.uns["nmf_genes"] = {}
    adata.uns["nmf_genes"]["cols"] = list(signatures)
    adata.uns["nmf_genes"]["rows"] = signatures.index
    adata.uns["nmf_genes"]["values"] = signatures.values
    adata.uns["nmf_genes"] = _uns_dict_to_df(adata, "nmf_genes")


def _add_markers_to_adata(adata, markers):

    adata.uns["nmf_markers"] = {}
    adata.uns["nmf_markers"]["cols"] = list(markers)
    adata.uns["nmf_markers"]["rows"] = markers.index
    adata.uns["nmf_markers"]["values"] = markers.values
    adata.uns["nmf_markers"] = _uns_dict_to_df(adata, "nmf_markers")


def _add_nmf_results_to_adata(
    adata, nmf_result, h5_path, cut_norm=0, cut_diff=0.1, silent=False, save=True
):

    """

    Parameters:
    -----------
    adata
    nmf_result
    h5_path
    cut_norm
    cut_diff
    silent
    save

    Returns:
    --------
    nmf_genes
        gene x signature table
        type: pandas.DataFrame

    nmf_markers
        cell x marker table
        type: pandas.DataFrame
    """

    aggr, best_run, h5_best, max_k_iter = _get_best_result(h5_path, silent, save)

    markers, signatures = sa.utils.select_markers(
        nmf_result["X"],
        nmf_result["run{}".format(best_run)]["W"],
        nmf_result["run{}".format(best_run)]["H"],
        cut_norm=cut_norm,
        cut_diff=cut_diff,
    )
    
    for result in nmf_result.values():
        if "consensus" in list(result.keys()):
            _add_consensus_clusters(adata, result)

    adata.obs = adata.obs.rename(columns={"max_id": "nmf_id"})
    _add_H_to_adata(adata, nmf_result["run{}".format(best_run)]["H"])
    _add_signatures_to_adata(adata, signatures)
    _add_markers_to_adata(adata, markers)
    
    return adata, aggr, best_run, h5_best
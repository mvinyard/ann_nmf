
import os
from ._select_best_result import _select_best_result
import signatureanalyzer as sa

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


def _add_markers_to_adata(adata, markers):

    adata.uns["nmf_markers"] = {}
    adata.uns["nmf_markers"]["cols"] = list(markers)
    adata.uns["nmf_markers"]["rows"] = markers.index
    adata.uns["nmf_markers"]["values"] = markers.values


def _add_nmf_results_to_adata(
    adata, nmf_result, results_path, cut_norm=0, cut_diff=0.1
):

    """"""

    best = _select_best_result(os.path.dirname(results_path))

    markers, signatures = sa.utils.select_markers(
        nmf_result["X"],
        nmf_result[best]["W"],
        nmf_result[best]["H"],
        cut_norm=cut_norm,
        cut_diff=cut_diff,
    )

    if "consensus" in list(nmf_result.keys()):
        _add_consensus_clusters(adata, nmf_result)

    adata.obs = adata.obs.rename(columns={"max_id": "nmf_id"})
    _add_H_to_adata(adata, nmf_result["run{}/H".format(best)])
    _add_signatures_to_adata(adata, signatures)
    _add_markers_to_adata(adata, markers)
    
    return adata
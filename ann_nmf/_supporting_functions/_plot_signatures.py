


from signatureanalyzer.plotting import marker_heatmap
from signatureanalyzer.utils import select_markers

import pandas as pd
import vinplots

def _linearize_axes(AxesDict):

    axes = []
    for ax_i in AxesDict.keys():
        for ax_j in AxesDict[ax_i].keys():
            axes.append(AxesDict[ax_i][ax_j])

    return axes


def _build_signature_plot(n_signatures, ncols):

    fig = vinplots.Plot()
    fig.construct(nplots=n_signatures, ncols=ncols)
    fig.modify_spines(ax="all", spines_to_delete=["top", "right"])
    axes = _linearize_axes(fig.AxesDict)

    return fig, axes


def _plot_signatures_umap(adata, ncols=4):

    umap = adata.obsm["X_umap"]
    signatures = adata.uns["signatures"]
    n_signatures = len(signatures)
    if n_signatures < ncols:
        ncols = n_signatures

    fig, axes = _build_signature_plot(n_signatures, ncols)

    for n, ax in enumerate(axes):
        sig = signatures[n]
        ax.scatter(umap[:, 0], umap[:, 1], c=adata.obs[sig], s=10)
        ax.set_title("Signature: {}".format(sig))
        ax.set_xticks([])
        ax.set_yticks([])

def _plot_signatures_heatmap(h5_out, best_h5, cut_norm=0, cut_diff=0.1):

    """"""

    X = pd.read_hdf(h5_out, "X")
    H = pd.read_hdf(best_h5, "H")
    W = pd.read_hdf(best_h5, "W")

    markers, signatures = select_markers(X, W, H, cut_norm=cut_norm, cut_diff=cut_diff)
    figure = marker_heatmap(X, signatures, H.sort_values("max_id").max_id)

    return figure
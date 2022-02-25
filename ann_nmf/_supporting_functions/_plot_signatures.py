
import pandas as pd

from signatureanalyzer.plotting import marker_heatmap
from signatureanalyzer.utils import select_markers


def _plot_signatures(h5_out, best_h5, cut_norm=0, cut_diff=0.1):

    """"""

    X = pd.read_hdf(h5_out, "X")
    H = pd.read_hdf(best_h5, "H")
    W = pd.read_hdf(best_h5, "W")

    markers, signatures = select_markers(X, W, H, cut_norm=cut_norm, cut_diff=cut_diff)
    figure = marker_heatmap(X, signatures, H.sort_values("max_id").max_id)

    return markers, signatures, figure
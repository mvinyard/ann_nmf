import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import signatureanalyzer


def _consensus_clustering(output_filepath, best_output_filepath):

    """
    Consensus clustering of ARD-NMF results.

    Parameters:
    -----------
    output_filepath
        type: str

    best_output_filepath
        type: str

    Returns:
    --------
    df
        consensus matrix from results
        type: pandas.DataFrame
    
    assign_p
        assignment probability for selected cluster
        type: pd.Series
    
    Notes:
    ------
    Modified from the original implementation of
    SignatureAnalyzer to account for both the saved "best" and saved "all". It's not necessarily
    better, just a band-aid for this implementation.
    """
    
    niter = signatureanalyzer.utils.get_nruns_from_output(output_filepath)
    H_selected = pd.read_hdf(best_output_filepath, "H")

    x = np.vstack(
        [
            pd.read_hdf(output_filepath, "run{}/H".format(i)).loc[:, "max_id"].values
            for i in range(niter)
        ]
    )
    consensus_matrix = np.vstack(
        [(x[:, [y]] == x[:]).sum(0) for y in range(x.shape[1])]
    )

    df = pd.DataFrame(
        consensus_matrix, index=H_selected.index, columns=H_selected.index
    )
    df = df.loc[
        H_selected.sort_values("max_id").index, H_selected.sort_values("max_id").index
    ]

    assign_p = pd.concat(
        [
            df.loc[
                H_selected[H_selected["max_id"] == x].index,
                H_selected[H_selected["max_id"] == x].index,
            ].mean(1)
            for x in set(H_selected["max_id"])
        ]
    )

    assign_p.name = "assignment"

    return df, assign_p
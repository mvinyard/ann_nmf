import os
import pandas as pd
import signatureanalyzer as sa

import warnings

warnings.filterwarnings("ignore")

from tqdm.notebook import tqdm

from ._ARDNMF_messages import _ARDNMF_messages


def _store_lam(nmf_result, h5_store, iteration):

    """"""

    lam = pd.DataFrame(data=nmf_result["lam"], columns=["lam"])
    lam.index.name = "K0"
    h5_store["run{}/lam".format(iteration)] = lam

    return h5_store


def _store_results(nmf_result, h5_store, iteration):

    """"""

    _store_lam(nmf_result, h5_store, iteration)

    for key in ["H", "W", "Hraw", "Wraw", "markers", "signatures", "log"]:
        h5_store["run{}/{}".format(iteration, key)] = nmf_result[key]

    return h5_store


def _run_ARD_NMF(mtx_df, n_runs=10, verbose=False, outdir="./", **nmf_kwargs):

    """"""

    h5_out = os.path.join(outdir, "nmf_output.h5")
    h5_store = pd.HDFStore(h5_out, "w")
    msg = _ARDNMF_messages(h5_out)
    msg.running()

    NMF_results = {}
    NMF_results["X"] = h5_store["X"] = mtx_df

    for iteration in tqdm(range(n_runs)):

        tag = "\t{}/{}: ".format(iteration, n_runs - 1)
        NMF_results["run{}".format(iteration)] = sa.ardnmf(
            mtx_df, tag=tag, verbose=verbose, **nmf_kwargs
        )
        h5_store = _store_results(NMF_results["run{}".format(iteration)], h5_store, iteration)
        if iteration == n_runs - 1:
            msg.saving_to()
            h5_store.close()

    return NMF_results, h5_out
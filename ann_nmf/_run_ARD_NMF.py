import os
import pandas as pd
import signatureanalyzer as sa

import warnings
warnings.filterwarnings("ignore")

from tqdm.notebook import tqdm

from ._ARDNMF_messages import _ARDNMF_messages

def _store_results(mtx_df, nmf_result, h5_store):

    """"""

    iteration = nmf_result["iter"]

    for key, value in nmf_result.items():
        if key is "lam":
            value = pd.DataFrame(data=value, columns=[key])
            value.index.name = "K0"
        elif key is "objective" or key is "iter":
            value = pd.Series(nmf_result[key])
        h5_store["run{}/{}".format(iteration, key)] = value
        
    return h5_store


def _run_ARD_NMF(mtx_df, n_runs, verbose, outdir, **nmf_kwargs):

    """"""
    
    msg = _ARDNMF_messages(outdir)
    msg.running()
    
    h5_store = pd.HDFStore(msg._h5_out, "w")
    NMF_results = {}
    NMF_results['X'] = h5_store["X"] = mtx_df
    
    for iteration in tqdm(range(n_runs)):
        
        tag = "\t{}/{}: ".format(iteration, n_runs - 1)
        nmf_result = sa.ardnmf(mtx_df, tag=tag, verbose=verbose, **nmf_kwargs)
        nmf_result["iter"] = iteration
        
        NMF_results['run{}/W'.format(iteration)] = nmf_result['W']
        NMF_results['run{}/H'.format(iteration)] = nmf_result['H']
        h5_store = _store_results(mtx_df, nmf_result, h5_store)
        NMF_results[iteration] = nmf_result
        if iteration == n_runs-1:
            msg.saving_to()
            h5_store.close()
            
    return NMF_results
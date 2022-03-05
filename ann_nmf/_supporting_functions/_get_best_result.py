import os
import pandas as pd
import signatureanalyzer as sa

from ._MessageModule import _MessageModule

def _select_best_result(h5_out, msg):

    aggr = sa.utils.get_nlogs_from_output(h5_out)
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr["K"] == max_k].shape[0]
    best_run = int(aggr[aggr["K"] == max_k].obj.idxmin())
    msg.identified_best_run(best_run, max_k_iter, aggr)
    
    return aggr, best_run, max_k_iter

def _store_best_run(h5_out, aggr, best_run, msg):
    
    h5_out_best_file = ".best.".join(os.path.basename(h5_out).split("."))
    h5_out_best = os.path.join(os.path.dirname(h5_out), h5_out_best_file)
    best_run_store = pd.HDFStore(h5_out_best,'a')
    msg.saving_best_run(h5_out_best)
    
    for key in ["H", "W", "lam", "Hraw", "Wraw", "markers", "signatures", "log"]:
        best_run_store[key] = pd.read_hdf(h5_out, "run{}/{}".format(best_run, key))
    
    best_run_store["aggr"] = aggr
    best_run_store.close()
    
    return h5_out_best
    
def _get_best_result(h5_out, silent=False, save=True):
    
    """
    
    Parameters:
    -----------
    h5_out
    silent
    save
    
    Returns:
    --------
    aggr
    
    best_run
    
    h5_out_best
    """
    
    msg = _MessageModule(h5_out, silent)    
    aggr, best_run, max_k_iter = _select_best_result(h5_out, msg)
    if save:
        h5_out_best = _store_best_run(h5_out, aggr, best_run, msg)
    
    return aggr, best_run, h5_out_best, max_k_iter
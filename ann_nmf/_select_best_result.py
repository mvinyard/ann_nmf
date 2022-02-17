import os
import signatureanalyzer as sa

def _select_best_result(outpath):

    aggr = sa.utils.get_nlogs_from_output(os.path.join(outpath, "nmf_output.h5"))
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr["K"] == max_k].shape[0]
    best_run = int(aggr[aggr["K"] == max_k].obj.idxmin())
    print(
        "Run {} had lowest objective with parameters:\n  n = {:g}\n  K = {:g}".format(
            best_run, max_k_iter, aggr.loc[best_run]["K"]
        )
    )

    return best_run

import licorice
import os
import pandas as pd

def _scanpy_highly_variable_genes_message(hv_key):

    """"""

    bold_txt1 = licorice.font_format(hv_key, ["BOLD"])
    bold_txt2 = licorice.font_format("sc.pp.highly_variable_genes", ["BOLD"])
    print(
        "hv_key: {} not present in adata.var. Computing highly-variable genes using {}.".format(
            bold_txt1, bold_txt2
        )
    )

class _MessageModule:
    def __init__(self, h5_out=False, silent=False):

        """"""
        
        self._h5_out = h5_out
        self._saving = "\n*Saving ARD-NMF outputs to: {}".format(self._h5_out)
        self._running = "\n*Running ARD-NMF...\n"
        self._silent = silent
        

    def saving_to(self):
        if not self._silent:
            print(self._saving)

    def running(self):
        if not self._silent:
            print(self._running)
        
    def saving_best_run(self, best_run_h5):
        if not self._silent:
            best_str = licorice.font_format("best", ['BOLD'])
            self._saving_best_run = "\n\t*Saving {} ARD-NMF outputs to: {}".format(best_str, best_run_h5)
            print(self._saving_best_run)

    def identified_best_run(self, best_run, max_k_iter, aggr):
        
        best_run_msg = best_run + 1
        
        if not self._silent:
            print(
            "\nRun {} had lowest objective with parameters:\n  n = {:g}\n  K = {:g}".format(
                best_run_msg, max_k_iter, aggr.loc[best_run]["K"]
            )
        )
            
            
    def scanpy_hv_genes(self, hv_key):
        _scanpy_highly_variable_genes_message(hv_key)
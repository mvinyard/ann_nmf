
import licorice
import os
import pandas as pd


class _ARDNMF_messages:
    def __init__(self, h5_out, silent=False):

        """"""
        
        self._h5_out = h5_out
        self._saving = "\t*Saving ARD-NMF outputs to: {}".format(self._h5_out)
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
            self._saving_best_run = "\t*Saving {} ARD-NMF outputs to: {}".format(best_str, best_run_h5)
            print(self._saving_best_run)

    def identified_best_run(self, best_run, max_k_iter, aggr):
        
        if not self._silent:
            print(
            "\nRun {} had lowest objective with parameters:\n  n = {:g}\n  K = {:g}".format(
                best_run, max_k_iter, aggr.loc[best_run]["K"]
            )
        )
import os
import pandas as pd


class _ARDNMF_messages:
    def __init__(self, outidr="./"):

        """"""

        self._h5_out = os.path.join(outidr, "nmf_output.h5")
        self._saving = "\t*Saving ARD-NMF outputs to: {}".format(self._h5_out)
        self._running = "\n*Running ARD-NMF...\n"

    def saving_to(self):
        print(self._saving)

    def running(self):
        print(self._running)
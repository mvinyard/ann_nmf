
from ._prepare_NMF_inputs import _prepare_NMF_inputs
from ._run_ARD_NMF import _run_ARD_NMF
from ._select_best_result import _select_best_result
from ._ARDNMF_messages import _ARDNMF_messages
from ._add_nmf_results_to_adata import _add_nmf_results_to_adata

class _ARD_NMF_GEX_wrapper:
    
    def __init__(
        self,
        adata,
        outdir="./",
        use_raw=True,
        layer=False,
        use_X=False,
        use_highly_variable=True,
        filter_mito=True,
        filter_ribo=True,
        **kwargs,
    ):

        """"""
        
        self._adata = adata
        
        self._mtx_df, self._genes_for_nmf = _prepare_NMF_inputs(
            self._adata,
            use_raw=use_raw,
            layer=layer,
            use_X=use_X,
            use_highly_variable=use_highly_variable,
            filter_mito=filter_mito,
            filter_ribo=filter_ribo,
            **kwargs,
        )
        
        self._outdir = outdir
        self._msg = _ARDNMF_messages(self._outdir)
        
    def run(self, n_runs=10, cut_norm=0, cut_diff=0.1, verbose=False, **nmf_kwargs):
                
        self._nmf_result = _run_ARD_NMF(self._mtx_df, n_runs, verbose, self._outdir, **nmf_kwargs)
        self._adata = _add_nmf_results_to_adata(self._adata, 
                                                self._nmf_result, 
                                                results_path=self._msg._h5_out, 
                                                cut_norm=cut_norm, 
                                                cut_diff=cut_diff
                                               )
        
    def get_best(self):
        
        """ """
        
        self._best = _select_best_result(self._outdir)
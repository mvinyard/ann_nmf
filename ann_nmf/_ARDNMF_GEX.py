
from ._supporting_functions._prepare_NMF_inputs import _prepare_NMF_inputs
from ._supporting_functions._run_ARD_NMF import _run_ARD_NMF
from ._supporting_functions._get_best_result import _get_best_result
from ._supporting_functions._add_nmf_results_to_adata import _add_nmf_results_to_adata
from ._supporting_functions._consensus_clustering import _consensus_clustering
from ._supporting_functions._plot_signatures import _plot_signatures

from signatureanalyzer.plotting import consensus_matrix
from signatureanalyzer.utils import get_nlogs_from_output

import pydk

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
        silent=False,
        save=True,
        **kwargs,
    ):

        """
        adata
        
        outdir
        
        use_raw
        
        layer
        
        use_X
        
        use_highly_variable
        
        filter_mito
        
        filter_ribo
        
        silent
        
        save
        """
        
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
        self._silent = silent
        self._save = save
        self._aggr = False
        
        pydk.mkdir_flex(self._outdir)
        
    def run(self, n_runs=10, cut_norm=0, cut_diff=0.1, verbose=False, **nmf_kwargs):
        
        """
        n_runs
        
        cut_norm
        
        cut_diff
        
        verbose
        
        **nmf_kwargs
        
        Run ARD-NMF from SignatureAnalyzer. 
        Finally, determines the best result post-run.
        
        """
                
        self._nmf_result, self._h5_out = _run_ARD_NMF(self._mtx_df, n_runs, verbose, self._outdir, **nmf_kwargs)
        self._adata, self._aggr, self._best_run, self._best_h5 = _add_nmf_results_to_adata(
            self._adata,
            self._nmf_result,
            h5_path=self._h5_out,
            cut_norm=cut_norm,
            cut_diff=cut_diff,
        )
        
    def get_best(self):
        
        """
        Determine the best result post-run.
        """
        
        self._aggr, self._best_run, self._best_h5 = _get_best_result(self._h5_out, self._silent, self._save)
        
    def cluster(self):
        
        """
        
        Parameters:
        -----------
        None
        
        Returns:
        --------
        Consensus cluster figure and `d`. Updates class with `self._ConsensusFigure` and `self._d`
        """
        
        
        self._cluster_df, self._assign_p = _consensus_clustering(self._h5_out, self._best_h5)
        if not type(self._aggr) == bool:
            self._aggr, self._best_run, self._best_h5 = _get_best_result(self._h5_out, self._silent, self._save)            
        self._max_k = self._aggr.groupby("K").size().idxmax()
        self._max_k_iter = self._aggr[self._aggr["K"] == self._max_k]['K'][0]
        self._ConsensusFigure, self._d = consensus_matrix(self._cluster_df, n_clusters=int(self._max_k_iter))
        
    def signatures(self, cut_norm=0, cut_diff=0.1):
        
        """
        Obtain and plot clustered gene signatures.
        
        Parameters:
        -----------
        cut_norm
            default: 0
            type: int
            
        cut_diff
            default: 0.1
            type: float
            
        Returns:
        --------
        Modifies the class in-place with the following:        
            self.NMF_markers
            self.NMF_genes
            self._GeneSignatureFigure
            
        """
        
        self._GeneSignatureFigure = _plot_signatures(self._h5_out, self._best_h5, cut_norm=cut_norm, cut_diff=cut_diff)
        self.nmf_markers = self._adata.uns["nmf_markers"]
        self.nmf_genes = self._adata.uns["nmf_genes"]
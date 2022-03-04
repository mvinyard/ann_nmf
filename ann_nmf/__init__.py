

# -------------------------------------------------------------------- #

from ._utilities._sklearn_compatibility import _sklearn_compatibility
_sklearn_compatibility()

# -------------------------------------------------------------------- #

from . import _utilities as ut

from ._ARDNMF_GEX import _ARD_NMF_GEX_wrapper as NMF

# -------------------------------------------------------------------- #
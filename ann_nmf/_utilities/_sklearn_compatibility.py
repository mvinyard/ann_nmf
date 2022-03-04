import sys
import sklearn
import sklearn.neighbors

def _sklearn_compatibility():
    
    """
    Solution from: https://stackoverflow.com/questions/60145652/no-module-named-sklearn-neighbors-base
    """
    
    sys.modules['sklearn.neighbors.base'] = sklearn.neighbors._base

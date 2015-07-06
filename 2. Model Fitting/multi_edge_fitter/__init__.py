import fitter_grav
import fitter_E
import fitter_pij
import fitter_k
import fitter_s
import fitter_sk

try:
	import fitter_s_CVXOPT
except ImportError:
	print "Please download and install CVXOPT suite for enhanced capabilities"
	CVXOPT_available = False
		
__all__ = [fitter_grav,fitter_E,fitter_pij,fitter_k,fitter_s,fitter_sk]


""" 
    Requires:
        Numpy
        Numba
        Scipy
        Pyproj
"""


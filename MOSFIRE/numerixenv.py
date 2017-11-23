from __future__ import division # confidence high
import os

def check_input(xxx):
    """Check if input is a Numarray Array."""
    try:
        import numarray
        return isinstance(xxx,numarray.numarraycore.NumArray)    
    except ImportError:
        pass

def check():
    """Check for running numarray version of pyfits with numpy code."""
    try:
        import pyfits
        if pyfits.__version__ < '1.1':
            raise EnvironmentError("Pyfits 1.1 or later required, pyfits version %s detected\n" % pyfits.__version__)
    except ImportError:
        pass
    try:
        if os.environ['NUMERIX'] == 'numarray':
            raise EnvironmentError("NUMERIX/numarray environment detected; numpy environment required")
    except KeyError:
        pass

    


'''
===================
MOSFIRE Options
===================

'''

import getpass
import os
__version__ = '2014.06.10'

npix = 2048
path_bpm = os.path.join(os.path.dirname(__file__), "data", "badpix_10sep2012.fits")

flat = {
        "version": 1, 
        "edge-order": 4, # Polynomial order for edge of slit 
        "edge-fit-width": 20,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}

wavelength = {
        "datadir" : os.path.join(os.path.dirname(__file__), "data"),
        "version": 2,
        "fractional-wavelength-search": 0.99935, # used in determining oned wavelength solutions
        "chebyshev-degree": 5, # polynomial order for fitting wavelengths
}


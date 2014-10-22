# Help, bugs to: http://mosfire.googlecode.com

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, \
    Rectify
from MOSFIRE import Wavelength

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

#Modify that maskname and band                                                                                                                           
# e.g. maskname = goodsnorth1                                                                                                                            
# e.g. band = H                                                                                                                                          
maskname = 'maskname'
band = 'band'

flatops = Options.flat
waveops = Options.wavelength

#Update the list of names for the offset files below.
obsfiles = ['Offset_1.5.txt', 'Offset_-1.5.txt']

#Flats.handle_flats('Flat.txt', maskname, band, flatops)

#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, obsfiles,waveops)
#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,waveops)
#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)

#update the "lambda_solution_wave_stack_K*.fits" file name                                                                                            
#  to the output file name from the apply_lambda process above.                                                                                 
#Background.handle_background(obsfiles,
    #'lambda_solution_wave_stack_*.fits',
    #maskname, band, waveops)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
#update the "lambda_solution_wave_stack_K*.fits" file name                                                                                            
#  to the output file name from the apply_lambda process above.                                                                                 
# Update the name of the first file in the offset file (use the full path name.                                                                          
#   e.g.    "/Users/user1/MOSFIRE/DRP_CODE/DATA/2014may08/m130114_0451.fits",                                                                        
#Rectify.handle_rectification(maskname, redfiles,
#    "lambda_solution_wave_stack_*.fits",
#    band, 
#    "/scr2/npk/mosfire/2013apr29/m130429_0224.fits",
#    waveops)
#


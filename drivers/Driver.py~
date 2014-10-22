# Help, bugs to: http://mosfire.googlecode.com

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, \
    Rectify
from MOSFIRE import Wavelength

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = ''
band = ''

flatops = Options.flat
waveops = Options.wavelength

obsfiles = ['Offset_1.5.txt', 'Offset_-1.5.txt']

#Flats.handle_flats('Flat.txt', maskname, band, flatops)
#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, obsfiles,
    #waveops)
#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,
    #waveops)

#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)
#Background.handle_background(obsfiles,
    #'lambda_solution_wave_stack_H_m130429_0224-0249.fits',
    #maskname, band, waveops)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
#Rectify.handle_rectification(maskname, redfiles,
#    "lambda_solution_wave_stack_H_m130429_0224-0249.fits",
#    band, 
#    "/scr2/npk/mosfire/2013apr29/m130429_0224.fits",
#    waveops)
#

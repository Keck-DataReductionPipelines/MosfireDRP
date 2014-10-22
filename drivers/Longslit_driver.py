
# Help, bugs to: http://mosfire.googlecode.com
#
# Instructions
#   1. edit band = '' to band = 'Y' or 'J' or 'H' or 'K'
#       e.g. band = 'J'
#   2. edit [709, 1350] to be the pixel values at the beginning and end 
#       of the long slit. Look at the raw data.
#   3. edit row_position to be a location where the standard star is not.
#   4. Decide if you want to use sky lines or Neon lamps for lambda calibration
#   5. Uncomment one line at a time and run mospy on the driver file
#

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify
from MOSFIRE import Wavelength, Longslit

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = 'longslit'
band = 'FIXME'

flatops = Options.flat
waveops = Options.wavelength

longslit = {'yrange': [709, 1350], 'row_position': 1158}

obsfiles = ['Offset_0.txt', 'Offset_13.33333.txt', 'Offset_6.666667.txt', 'Offset_-13.33333.txt',  'Offset_-6.666667.txt']

#Flats.handle_flats('Flat.txt', maskname, band, flatops)

# SKY LINES
# Use sky lines for wavelength calibration
# Use either Sky lines or Neon lines.
#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, obsfiles, waveops, longslit=longslit)
#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles, waveops, longslit=longslit)
#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops, longslit=longslit, smooth=True)

# Update the following line after the apply_lambda_simple step
#Longslit.go(maskname, band, obsfiles ,
#    'lambda_solution_wave_stack_....fits',
#    waveops, longslit)


# NEON
# Use neon for wavelength calibrations
#Wavelength.imcombine('Ne.txt', maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, 'Ne.txt', waveops, longslit=longslit, neon=True)
#Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt', waveops, longslit=longslit)
#Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops, longslit=longslit, smooth=True)

# Update the following line after the apply_lambda_simple step
#Longslit.go(maskname, band, obsfiles,
#    'lambda_solution_wave_stack_....fits',
#    waveops, longslit)




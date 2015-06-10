
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

maskname = 'LONGSLIT-46x0.7'
band = 'J'

flatops = Options.flat
waveops = Options.wavelength


longslit = {'yrange': [900, 1150], 'row_position': 1000, 'mode':'longslit'}

obsfiles = ['Offset_15_HD2050_wise233.txt', 'Offset_-15_HD2050_wise233.txt']

#Flats.handle_flats('Flat.txt', maskname, band, flatops, longslit = longslit)

# SKY LINES
# Use sky lines for wavelength calibration
# Use either Sky lines or Neon lines.
#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, obsfiles, waveops, longslit=longslit)
#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles, waveops, longslit=longslit)
#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops, longslit=longslit, smooth=True)


#Background.handle_background(obsfiles,
#    'lambda_solution_wave_stack_J_m150609_0654-0661.fits',
#    maskname, band, waveops)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
#update the "lambda_solution_wave_stack_K*.fits" file name
#  to the output file name from the apply_lambda process above.
# Update the name of the first file in the offset file (use the full path name.
#   e.g.    "/Users/user1/MOSFIRE/DRP_CODE/DATA/2014may08/m130114_0451.fits",

#Rectify.handle_rectification(maskname, redfiles,
#    "lambda_solution_wave_stack_J_m150609_0654-0661.fits",
#    band, 
#    obsfiles,
#    waveops)


# NEON
# Use neon for wavelength calibrations
#Wavelength.imcombine('Ne.txt', maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, 'Ne.txt', waveops, longslit=longslit, neon=True)
#Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt', waveops, longslit=longslit)
#Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops, longslit=longslit, smooth=True)
#print redfiles
#obsfiles = ['eps_off_-10.txt', 'eps_off_10.txt']
# Update the following line after the apply_lambda_simple step
#Longslit.go(maskname, band, obsfiles,
#    'lambda_solution_wave_stack_H_m150112_0199-0201.fits',
#    waveops, longslit, extension='/Volumes/PromiseRAID/MOSFIRE/DRP_CODE/DATA/2015jan12/m150112_0199.fits')




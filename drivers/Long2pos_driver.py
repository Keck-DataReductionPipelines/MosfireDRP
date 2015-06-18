
import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify
from MOSFIRE import Wavelength, Longslit

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = 'long2pos'
band = 'H'

flatops = Options.flat
waveops = Options.wavelength

# Note: for long2pos, the row position is ignored, and the middle point of the slit is used
longslit = {'yrange': [[1062,1188],[887,1010]], 'row_position': 0, 'mode':'long2pos'}


# SETUP FILES FOR DATES before June 10, 2015

obsfiles_posC_narrow = ['Offset_-21_HIP85871_7.25_PosC.txt', 'Offset_-7_HIP85871_7.25_PosC.txt']
targetCnarrow = "HIP85871_posC_narrow"

obsfiles_posA_narrow = ['Offset_7_HIP85871_7.25_PosA.txt', 'Offset_21_HIP85871_7.25_PosA.txt']
targetAnarrow = "HIP85871_posA_narrow"

obsfiles_posC_wide = ['Offset_-14_HIP85871_7.25_PosC.txt','Offset_-7_HIP85871_7.25_PosC.txt']
targetCwide = "HIP85871_posC_wide"

obsfiles_posA_wide = ['Offset_14_HIP85871_7.25_PosA.txt','Offset_21_HIP85871_7.25_PosA.txt']
targetAwide = "HIP85871_posA_wide"

# SETUP FILES for DATES after June 10,2015 

obsfiles_posC_narrow = ['Offset_7_FS134_posC.txt','Offset_-7_FS134_PosC.txt']
targetCnarrow = "FS134_posC_narrow"
obsfiles_posA_narrow = ['Offset_7_FS134_posA.txt','Offset_-7_FS134_PosA.txt']
targetAnarrow = "FS134_posA_narrow"
obsfiles_posC_wide = ['Offset_0_FS134_posC.txt','Offset_-7_FS134_PosC.txt']
targetCwide = "FS134_posC_wide"
obsfiles_posA_wide = ['Offset_0_FS134_posA.txt','Offset_-7_FS134_PosA.txt']
targetAwide = "FS134_posA_wide"

# Argon files
argon = ['Ar.txt']

Flats.handle_flats('Flat.txt', maskname, band, flatops, longslit = longslit)
# if you have thermal flats in K...
#Flats.handle_flats('Flat.txt', maskname, band, flatops,lampOffList='FlatThermal.txt', longslit=longslit)

# Uses the argon calibration taken in the afternoon with long2pos for the wavelength calibration

Wavelength.imcombine(argon, maskname, band, waveops)
Wavelength.fit_lambda_interactively(maskname, band, argon, waveops, longslit=longslit, argon=True)
Wavelength.fit_lambda(maskname, band, argon, argon, waveops, longslit=longslit)
Wavelength.apply_lambda_simple(maskname, band, argon, waveops, longslit=longslit, smooth=True)


# UPDATE THE RESULTING WAVELENGTH SOLUTION HERE!
wavelength_file = 'lambda_solution_wave_stack_H_m150428_0091-0091.fits'

########### NARROW SLITS ############

# Long2POS: NARROW SLITS:  use the -7 and -21 for posC and 7 and 21 for posA

obsfiles = obsfiles_posC_narrow
target = targetCnarrow

#IO.fix_long2pos_headers(obsfiles)

#Background.handle_background(obsfiles,
#   wavelength_file,
#    maskname, band, waveops,  target=target)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]

#Rectify.handle_rectification(maskname, redfiles,
#    wavelength_file,
#    band, 
#    obsfiles,
#    waveops, target=target)


obsfiles = obsfiles_posA_narrow
target = targetAnarrow

#IO.fix_long2pos_headers(obsfiles)

#Background.handle_background(obsfiles,
#        wavelength_file,
#    maskname, band, waveops,  target=target)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]

#Rectify.handle_rectification(maskname, redfiles,
#    wavelength_file,
#    band, 
#    obsfiles,
#   waveops, target=target)

########### WIDE SLITS ############

# SPECTROPHOTOMETRIC Long2POS: THIS SECTION IS FOR THE REDUCTION OF THE WIDE SLITS.

#obsfiles = obsfiles_posC_wide
#target = targetCwide

#IO.fix_long2pos_headers(obsfiles)

#Background.handle_background(obsfiles,
#    wavelength_file,
#    maskname, band, waveops,  target=target)

#obsfiles = obsfiles_posA_wide
#target = targetAwide

#IO.fix_long2pos_headers(obsfiles)

#Background.handle_background(obsfiles,
#    wavelength_file,
#    maskname, band, waveops,  target=target)






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




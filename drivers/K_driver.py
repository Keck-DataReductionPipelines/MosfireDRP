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
# e.g. band = K
maskname = 'maskname'
band = 'K'

flatops = Options.flat
waveops = Options.wavelength

#Flats.handle_flats('Flat.txt', maskname, band, flatops,lampOffList='FlatThermal.txt')

#Update the list of names for the offset files below.
obsfiles = ['Offset_2.txt', 'Offset_-2.txt']
#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.imcombine('Ne.txt', maskname, band, waveops)
#Wavelength.imcombine('Ar.txt', maskname, band, waveops)

#Wavelength.fit_lambda_interactively(maskname, band, obsfiles,     waveops, bypass=False)
#Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ne.txt', neon=True)
#Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ar.txt', argon=True)

#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles, waveops)
#Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt', waveops, wavenames2='Ar.txt')
#LROI = [[21000, 22800]] * 1
#LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles, 'Ne.txt', LROI, waveops)

#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)
#Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops)
#Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ne.txt', LROIs, waveops, neon=True)

#update the "merge_lambda_solution_wave_stack_K*.fits" file name
#  to the output file name from the apply_lambda_sky_and_arc process above.
#Background.handle_background(obsfiles,
#    'merged_lambda_solution_wave_stack_K_*.fits',
#    maskname,
#    band,
#    waveops)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
#update the "merge_lambda_solution_wave_stack_K*.fits" file name
#  to the output file name from the apply_lambda_sky_and_arc process above.
# Update the name of the first file in the offset file (use the full path name.
#   e.g.    "/Users/user1/MOSFIRE/DRP_CODE/DATA/2014may08/m130114_0451.fits",
#Rectify.handle_rectification(maskname, redfiles,
#    'merged_lambda_solution_wave_stack_K_*.fits',
#   band, 
#    "full_path_to_first_file_in_offsets.fits",
#   waveops)



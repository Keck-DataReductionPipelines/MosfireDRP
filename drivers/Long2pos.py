import os, time, logging
import MOSFIRE
from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength
from MOSFIRE.MosfireDrpLog import info, debug, warning, error
logger = logging.getLogger(__name__)
import numpy as np
from matplotlib import pyplot as pl
import pyfits as pf
np.seterr(all='ignore')
flatops = Options.flat
waveops = Options.wavelength

#Driver file automatically generated on Wed Jul 29 15:04:02 2015
#For questions and comments, email mosfiredrp@gmail.com, submit a ticket on the ticketing system, or contact Luca Rizzi @ WMKO

# THIS DRIVER IS SETUP TO REDUCE A STAR CALLED HIP85871_7.25, change this string on the entried below



maskname = 'long2pos_specphot (align)'
band = 'H'

#Set bypass to True to autofit wavelenth solution instead of manually fitting.
bypassflag=False

# these are the narrow slits
obsfiles_posCnarrow = ['Offset_-21_HIP85871_7.25_PosC.txt', 'Offset_-7_HIP85871_7.25_PosC.txt']
target_posCnarrow = "HIP85871_7.25_POSC_NARROW"
IO.fix_long2pos_headers(obsfiles_posCnarrow)
obsfiles_posAnarrow = ['Offset_7_HIP85871_7.25_PosA.txt', 'Offset_21_HIP85871_7.25_PosA.txt']
target_posAnarrow = "HIP85871_7.25_POSA_NARROW"
IO.fix_long2pos_headers(obsfiles_posAnarrow)
# these are the wide slits, comment out if you are not using specphot
obsfiles_posCwide = ['Offset_-14_HIP85871_7.25_PosC.txt', 'Offset_-7_HIP85871_7.25_PosC.txt']
target_posCwide = "HIP85871_7.25_POSC_WIDE"
IO.fix_long2pos_headers(obsfiles_posCwide)
obsfiles_posAwide = ['Offset_14_HIP85871_7.25_PosA.txt', 'Offset_21_HIP85871_7.25_PosA.txt']
target_posAwide = "HIP85871_7.25_POSA_WIDE"
IO.fix_long2pos_headers(obsfiles_posAwide)

# Note: for long2pos, the row position is ignored, and the middle point of the slit is used
longslit = {'yrange': [[1062,1188],[887,1010]], 'row_position': 0, 'mode':'long2pos'}
Flats.handle_flats('Flat.txt', maskname, band, flatops,longslit=longslit)

# in this case, we are using the argon lines.
# replace this with neon=['Ne.txt'] if you prefer to use Ne, and edit the following lines accordingly
argon = ['Ar.txt']
Wavelength.imcombine(argon, maskname, band, waveops)
Wavelength.fit_lambda_interactively(maskname, band, argon,waveops,longslit=longslit, argon=True, bypass=bypassflag)
Wavelength.fit_lambda(maskname, band, argon,argon,waveops,longslit=longslit)
Wavelength.apply_lambda_simple(maskname, band, argon, waveops, longslit=longslit, smooth=True)

# make sure you use the correct wavelength file generated before
Wavelength_file = 'lambda_solution_wave_stack_H_m150428_0091-0091.fits'

# narrow
Background.handle_background(obsfiles_posAnarrow,Wavelength_file,maskname,band,waveops, target=target_posAnarrow)
Background.handle_background(obsfiles_posCnarrow,Wavelength_file,maskname,band,waveops, target=target_posCnarrow)
# wide
Background.handle_background(obsfiles_posAwide,Wavelength_file,maskname,band,waveops, target=target_posAwide)
Background.handle_background(obsfiles_posCwide,Wavelength_file,maskname,band,waveops, target=target_posCwide)

# narrow
redfiles = ["eps_" + file + ".fits" for file in obsfiles_posAnarrow]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_posAnarrow,waveops, target=target_posAnarrow)
redfiles = ["eps_" + file + ".fits" for file in obsfiles_posCnarrow]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_posCnarrow,waveops, target=target_posCnarrow)
# wide
redfiles = ["eps_" + file + ".fits" for file in obsfiles_posAwide]
redfiles = [redfiles[0]]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_posAwide,waveops, target=target_posAwide)
redfiles = ["eps_" + file + ".fits" for file in obsfiles_posCwide]
redfiles = [redfiles[0]]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_posCwide,waveops, target=target_posCwide)


import os, time, logging
import MOSFIRE
from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength
from MOSFIRE.MosfireDrpLog import info, debug, warning, error
logger = logging.getLogger(__name__)
import numpy as np, pylab as pl, pyfits as pf
np.seterr(all='ignore')
flatops = Options.flat
waveops = Options.wavelength

#Driver file automatically generated on Sat Jul 25 17:46:43 2015
#For questions and comments, email mosfiredrp@gmail.com, submit a ticket on the ticketing system, or contact Luca Rizzi @ WMKO

maskname = 'maskname'
band = 'band'

#Set bypass to True to autofit wavelenth solution instead of manually fitting.
bypassflag=False
obsfiles=['Offset_1.25.txt','Offset_-1.25.txt']

Flats.handle_flats('Flat.txt', maskname, band, flatops,lampOffList='FlatThermal.txt')

Wavelength.imcombine(obsfiles, maskname, band, waveops)
# if you have Ar
Wavelength.imcombine('Ar.txt', maskname, band, waveops)
# if you have Ne
Wavelength.imcombine('Ne.txt', maskname, band, waveops)
Wavelength.fit_lambda_interactively(maskname, band, obsfiles,waveops, bypass=bypassflag)

Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ar.txt', argon=True)
Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ne.txt', neon=True)

Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,waveops)
Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt',waveops, wavenames2='Ar.txt')
LROI = [[21000,22800]]*1
LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles, 'Ne.txt', LROI, waveops)

Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)
Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ne.txt', LROIs, waveops)

Wavelength_file = 'merged_lambda_solution_wave_stack_K_*.fits'

Background.handle_background(obsfiles,Wavelength_file,maskname,band,waveops)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles,waveops)


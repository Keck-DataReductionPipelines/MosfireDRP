import os, time, logging
import MOSFIRE
from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength
from MOSFIRE.MosfireDrpLog import info, debug, warning, error
logger = logging.getLogger(__name__)
import numpy as np
from matplotlib import pyplot as pl
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
np.seterr(all='ignore')
flatops = Options.flat
waveops = Options.wavelength

#Driver file automatically generated on Wed Jul 29 15:11:22 2015
#For questions and comments, email mosfiredrp@gmail.com, submit a ticket on the ticketing system, or contact Luca Rizzi @ WMKO

maskname = 'LONGSLIT-3x0.7'
band = 'H'

#Set bypass to True to autofit wavelenth solution instead of manually fitting.
bypassflag=False
# modify the target name to match your observations
obsfiles=['Offset_5_HIP17971.txt','Offset_-5_HIP17971.txt']
target="HIP17971"


# modify the yrange to match the size of your longslit
# row position is the extraction line used for the initial wavelength solution. It should be away from your target
longslit = {'yrange':[968,1100],'row_position':1034,'mode':'longslit'}
Flats.handle_flats('Flat.txt', maskname, band, flatops,longslit=longslit)

Wavelength.imcombine(obsfiles, maskname, band, waveops)
Wavelength.fit_lambda_interactively(maskname, band, obsfiles,waveops,longslit=longslit, bypass=bypassflag)
Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,waveops,longslit=longslit)
Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops,longslit=longslit)

# make sure you use the file generated on the previous step
Wavelength_file = 'lambda_solution_wave_stack_H_m121227_0162-0311.fits'

Background.handle_background(obsfiles,Wavelength_file,maskname,band,waveops,target=target)

redfiles = ["eps_" + file + ".fits" for file in obsfiles]
Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles,waveops, target=target)


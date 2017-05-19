"""

==================================
MOSFIRE Wavelength Calibrations

This module is responsible for determining wavelength solutions.

The wavelength is represented as a physical model of MOSFIRE using the 
grating equation, the size of a pixel (18 micron), and the focal length
of the camera (250 mm), the lines per mm of the grating (110.5), and the
order of the grating (Y: 6, J: 5, H: 4, K: 3).

formally, a fifth-order chebyshev polynomial is fit as the wavelength solution,
though the order is specified through Options.py

 scale is (pixel size) / (camera focal length)

-- Helper functions also exist for determining the on-order region of 
a spectrum --


control flow

    1. fit a spatially central pixel on a spectrum interactively. Goal is to
    achieve RMS of 0.1 Angstrom or better. The entry point is a class called
    InteractiveSolution called from a wrapper called fit_lambda_interactively.
    Interactive solution is a visual wrapper around the functions
    find_known_lines & fit_chebyshev_to_lines. Currently an old optical model
    is used to guestimate the wavelength solution but this could be updated in
    the future. 
        -> Produces lambda_center_coeffs_maskname.npy

    2. based on interactive fits, perform automated fits "outwards" from the
    central pixel in the spectrum. These fits are performed using the final
    linelist from the interactive fit. The entry point is a function called
    fit_lambda that iteratively calls fit_outwards_refit.

        -> Produces lambda_2d_coeffs_maskname.npy

    3. Apply these lambda fits to produce a full wavelength solution

INPUT:

OUTPUT:

npk Apr/May  2012 - Significant enhancements w/ first light data
npk April 26 2011
npk   May  4 2011

"""

import logging as log
from multiprocessing import Pool
import os
import itertools
import time

import numpy as np
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
from matplotlib import pyplot as pl

from scipy.interpolate import interp1d
from scipy import signal 
from scipy import optimize
from matplotlib.widgets import Button
from numpy.polynomial import chebyshev as CV


from MOSFIRE import CSU, Fit, IO, Options, Filters, Detector
from MosfireDrpLog import debug, info, warning, error

import pdb

__version__ = "1May2012"

MADLIMIT = 0.1

try:
    __IPYTHON__
    reload(Options)
    reload(CSU)
    reload(IO)
    reload(Fit)
    reload(Filters)
except:
    pass

#
# Glue code
#

def filelist_to_wavename(files, band, maskname, options):
    start = files[0].split('/')[-1].rstrip(".fits")
    end = files[-1].split('/')[-1].rstrip(".fits").split("_")[1]
    name = "wave_stack_{0}_{1}-{2}.fits".format(band, start, end)

    return name


def grating_results(band):
    '''returns the dlambda/dpixel in angstrom for a band'''
    orders = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = orders[band]
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length
    dlambda = scale * d / order * 10000

    return dlambda


def filelist_to_path(files, band, maskname, options):
    outf = filelist_to_wavename(files, band, maskname, options)

    return outf

def imcombine(files, maskname, bandname, options, extension=None):
    ''' This version of imcombine is used to create the wave_stack file
    which is used only by the wavelength fitting routine. This imcombine does
    not produce science results.

    Args:
        files: list of strings with file names
        maskname: string with name of the mask
        bandname: string with name of band
        options: passed across from Options file
        extension: path to file that contains a well formated fits header
            this should be used only when the detector server fails
            to write the full FITS header
    
    Results:
        writes a median combined image in electron. It is called
            wave_stack_[bandname]_[filename range].fits'''
    
    pixelflat_file = "pixelflat_2d_{0}.fits".format(bandname)
    flat = IO.readfits(pixelflat_file, use_bpm=True)[1]
    flat = flat.filled(1.0)

    files = IO.list_file_to_strings(files)

    info("combining Wavelength files")
    for file in files:
        debug(str(file))
    debug("list complete")

    ADUs = np.zeros((len(files), 2048, 2048))
    prevssl = None
    prevmn = None
    patternid = None
    header = None


    for i in xrange(len(files)):
        fname = files[i]
        thishdr, data, bs = IO.readmosfits(fname, options, extension=extension)
        info("Checking maskname and filter for {} {}/{}".format(fname, maskname, thishdr['filter']))
        ADUs[i,:,:] = data.filled(0)

        if thishdr["aborted"]:
            raise Exception("Img '%s' was aborted and should not be used" %
                    fname)

        if prevssl is not None:
            if len(prevssl) != len(bs.ssl):
                # todo Improve these checks
                info( "Reading file "+str(fname))
                error("This file contains "+str(len(bs.ssl))+" slits instead of "+str(len(prevssl)))
                raise Exception("The stack of input files seems to be of "
                        "different masks")
        prevssl = bs.ssl

        if maskname is not None:
            if maskname != thishdr["maskname"]:
                warning("Maskname specified ({0}) does not match header maskname "
                    " ({1}).".format(maskname, thishdr["maskname"]))

        if thishdr["BUNIT"] != "ADU per coadd":
            error("The units of '%s' are not in ADU per coadd and "
                    "this violates an assumption of the DRP. Some new code " 
                    "is needed in the DRP to handle the new units of "
                    "'%s'." % (fname, thishdr["BUNIT"]))
            raise Exception("The units of '%s' are not in ADU per coadd and "
                    "this violates an assumption of the DRP. Some new code " 
                    "is needed in the DRP to handle the new units of "
                    "'%s'." % (fname, thishdr["BUNIT"]))

        ''' Construct Header'''
        if header is None:
            header = thishdr

        header.set("imfno%2.2i" % (i), fname)

        for key in header.keys():
            try: val = header[key]
            except KeyError: 
                warning("Header should have key '%s' but does not" % key)
                print("Header should have key '%s' but does not" % key)
            
            if key in thishdr:
                if val != thishdr[key]:
                    newkey = "hierarch " + key + ("_img%2.2i" % i)
                    try: header.set(newkey.rstrip(), thishdr[key])
                    except ValueError: pass

        ''' Now handle error checking'''

        if maskname is not None:
            if thishdr["maskname"] != maskname:
                warning("File %s uses mask '%s' but the stack is of '%s'" %
                    (fname, thishdr["maskname"], maskname))

        for key in header.keys():
            val = header[key]

            if key in thishdr:
                if val != thishdr[key]:
                    newkey = "hierarch " + key + ("_img%2.2i" % i)
                    try: header.set(newkey.rstrip(), thishdr[key])
                    except ValueError: pass

        ''' Now handle error checking'''

        if maskname is not None:
            if thishdr["maskname"] != maskname:
                warning("File %s uses mask '%s' but the stack is of '%s'" %
                    (fname, thishdr["maskname"], maskname))

        info('Done.')

    wavename = filelist_to_wavename(files, bandname, maskname, options)
    info('Combining images to make {}'.format(wavename))
    header.set("frameid", "median")
    electrons = np.median(np.array(ADUs) * Detector.gain, axis=0)
    IO.writefits(electrons, maskname, wavename, options, overwrite=True,
            header=header)
    info("Done")
    
def fit_lambda(maskname, 
        bandname, 
        wavenames, 
        guessnames, 
        options,
        longslit=None,
        neon=None,
        extension=None,
        wavenames2=None):

    """Fit the two-dimensional wavelength solution to each science slit

    Inputs:
        bandname: The mask name as a string
        wavenames: List of wavelength files
        options: Dictionary of wavelength options
        longlist: True if a longslit
        neon: path to neon image [2k x 2k frame]
        extension: path to file that contains a well formated fits header
            this should be used only when the detector server fails
            to write the full FITS header

    Prints:
        This step prints out lines like
            resid ang S01 @ p####: 0.25 rms 0.15 mad [shift10]

        which means that the residual error in angstrom units for slit #1,
        at pixel location #### is 0.25 angstrom RMS and 0.15 median absolute
        deviation. The shift refers to the average amount of pixel shift from
        the solution determined during the interactive fitting step. 

    Results:
        
        lambda_coeffs_...npy: Coefficients file containing an array of 
        dictionaries:
            {"slitno": The 0-indexed slit number into the barset
            "lines": The list of fitted emission lines
            "center_sol": The Chebyshev coefficients in the center of the 
                slit
            "2d": A dictionary containing:
                {'positions': The rows that comprise the slit
                'delts:' The mean standard deviation of features. 
                    This list is the length of 'lines'
                'lambdaRMS': RMS of delts
                'lambdaMAD': MAD of delts
                'coeffs': The Chebyshev coefficients for each position}
            }

                
        lambda_solution....fits: A fits file with wavelength [Ang] for each
            pixel in the image
        sigs_solution...fits: A fits file with the standard deviation per row
            of the wavelength solution
        rectified_wave_stack: The tilted spectra are interpolated onto a common
            wavelength grid.
"""
    global bs, data, lamout, center_solutions, edgedata, data2, center_solutions2 
    np.seterr(all="ignore")

    """ Set defaults for second set of lines to None """
    data2=None
    center_solutions=None    

    wavenames = IO.list_file_to_strings(wavenames)
    debug("WAVENAMES"+ str(wavenames))
    wavename = filelist_to_wavename(wavenames, bandname, maskname,
            options).rstrip(".fits")

    guessnames = IO.list_file_to_strings(wavenames)
    debug("GUESSNAMES"+str(guessnames))
    guessname = filelist_to_wavename(guessnames, bandname, maskname,
            options).rstrip(".fits")

    fn = "lambda_coeffs_{0}.npy".format(wavename)

    info("%s] Writing to: %s" % (maskname, fn))

    wavepath = filelist_to_path(wavenames, bandname, maskname,
            options)
    drop, data = IO.readfits(wavepath, use_bpm=True)
    header, drop,bs = IO.readmosfits(wavenames[0], options, extension=extension)

    fnum = guessname
    center_solutions = IO.load_lambdacenter(fnum, maskname, options)
    edgedata, metadata = IO.load_edges(maskname, bandname, options)

    """ This neon flag looks like it may be removed -MK 2014 June 10 """
    if neon is not None:
        drop, Neon = IO.readfits(neon, use_bpm=True)
        data += Neon

    if wavenames2 is not None:
        wavenames2 = IO.list_file_to_strings(wavenames2)
        debug("WAVENAMES Second Set:"+str(wavenames2))
        wavename2 = filelist_to_wavename(wavenames2, bandname, maskname,
                                         options).rstrip(".fits")

        guessnames2 = IO.list_file_to_strings(wavenames2)
        debug("GUESSNAMES Second Set:"+str(guessnames2))
        guessname2 = filelist_to_wavename(guessnames2, bandname, maskname,
                                         options).rstrip(".fits")

        wavepath2 = filelist_to_path(wavenames2, bandname, maskname,
                                    options)
        drop, data2 = IO.readfits(wavepath2, use_bpm=True)
        header2, drop,bs2 = IO.readmosfits(wavenames2[0], options, extension=extension)

        fnum2 = guessname2
        center_solutions2 = IO.load_lambdacenter(fnum2, maskname, options)

    if longslit is not None and longslit['mode'] is "longslit":
        info("*** Longslit mode ***  Slitedges set to:")
        info("Bottom: "+str(edgedata[0]["yposs_bot"][0]))
        info("Top:    "+str(edgedata[0]["yposs_top"][0]))

    solutions = []
    lamout = np.zeros(shape=(2048, 2048), dtype=np.float32)

    tock = time.time()

    multicore = False
    if multicore:
        p = Pool()
        solutions = p.map(fit_lambda_helper, range(len(bs.ssl)))
        p.close()
    else:
        solutions = map(fit_lambda_helper, range(len(bs.ssl)))

    tick = time.time()

    info("-----> Mask took %i" % (tick-tock))

    try: os.remove(fn)
    except: pass
    np.save(fn, solutions)

    return solutions


def fit_lambda_helper(slitno):
    """This helper function exists for multiprocessing suport"""
    
    global bs, data, lamout, center_solutions, edgedata, data2, center_solutions2

    slitidx = slitno-1

    tick = time.time()

    slitedges = edgedata

    sol_1d = center_solutions[slitidx]["sol_1d"]
    edge = slitedges[slitidx]
    linelist = center_solutions[slitidx]["linelist"]

    start  = center_solutions[slitidx]["extract_pos"]
    bottom = np.ceil(edge["bottom"](1024))+2
    top    = np.ceil(edge["top"](1024))-2

    info(("* Fitting Slit %s from %i to %i" % (bs.ssl[slitno]["Target_Name"],
        bottom, top)))

    if data2 is not None:
        sol_1d2 = center_solutions2[slitidx]["sol_1d"]
        linelist2 = center_solutions2[slitidx]["linelist"]
        sol_2d = fit_outwards_refit(data, bs, sol_1d, linelist, Options.wavelength,
                                    start, bottom, top, slitno, data2=data2, linelist2=linelist2, sol_1d2=sol_1d2)
    else:
        sol_2d = fit_outwards_refit(data, bs, sol_1d, linelist, Options.wavelength,
                                    start, bottom, top, slitno)

    sol = {"slitno": slitno, "center_sol": np.array(sol_1d[1]), "2d":
        sol_2d, "lines": np.array(linelist)}
    
    info("S%2.2i] TOOK: %i s" % (slitno, time.time()-tick))

    return sol

def apply_interactive(maskname, band, options, apply=None, to=None, neon=False,
    argon=False, extension=None,short_exp = False):
    """Fit the one-dimensional wavelength solution to each science slit"""
    np.seterr(all="ignore")

    # Load the guess wavelength solution data
    wavenames = IO.list_file_to_strings(apply)
    wavename = filelist_to_path(wavenames, band, maskname, options)
    fn = "lambda_center_coeffs_{0}.npy".format(wavename.rstrip(".fits"))
    waves = np.load(fn)

    # Load the arc lamp
    to_files = IO.list_file_to_strings(to)
    to_filename = filelist_to_path(to_files, band, maskname, options)
    mfits = IO.readfits(to_filename, use_bpm=True)
    (drop, data) = mfits
    (header, drop, bs) = IO.readmosfits(wavenames[0], options, extension=extension)

    mfits = header, data, bs
    linelist = pick_linelist(header, neon=neon, argon=argon, short_exp = short_exp)

    solutions = []
    pix = np.arange(2048)
    for slitno in xrange(len(waves)):
        info("Slit number = "+str(slitno+1))
        csuslits = bs.scislit_to_csuslit(slitno+1)

        try:
            l = len(csuslits)
            if l > 1:
                csuslit = csuslits[l/2]
            else:
                csuslit = csuslits[0]
        except:
            csuslit = csuslits

        extract_pos = bs.science_slit_to_pixel(slitno+1)
        cfit = waves[slitno]['sol_1d'][1]

        spec = \
            np.ma.mean(data[extract_pos-1:extract_pos+1, :],
                axis=0) # axis = 0 is spatial direction

        STD = np.inf
        n_attempts = 5
        while (n_attempts > 0) and (STD > .3) :
            ll = CV.chebval(pix, cfit)
            [xs, sxs, sigmas] = find_known_lines(linelist, ll, spec, options)
            [deltas, cfit, perror] = fit_chebyshev_to_lines(xs, sxs, linelist, options)

            ok  = np.isfinite(deltas)
            STD = np.std(deltas[ok])
            MAD = np.median(np.abs(deltas[ok]))
            n_attempts -= 1

        solutions.append({"linelist": linelist, "MAD": MAD, 
                "foundlines": xs, "foundlinesig": sxs,
                "sol_1d": [deltas, cfit, sigmas], "STD":
                STD, "slitno": slitno, "extract_pos":
                extract_pos})

        info("slitno %2.0i STD: %1.2f MAD: %1.2f" % (slitno+1, STD, MAD))


    # Output filename
    info(str(to_filename))
    outfn = "lambda_center_coeffs_{0}.npy".format(to_filename.rstrip(".fits"))

    np.save(outfn, solutions)

def check_wavelength_roi(maskname, band, skyfiles, arcfiles, LROI, options, no_check=False):
    '''The purpose of this function is to help the user selection a wavelength
        range of interest over which to normalize the arcs versus sky solutions.
        '''
    skyfiles = IO.list_file_to_strings(skyfiles)
    skyfilename = filelist_to_path(skyfiles, band, maskname, options)
    fn = "lambda_center_coeffs_{0}.npy".format(skyfilename.rstrip(".fits"))
    skysols = np.load(fn)

    # Load the arc wavelength solution data
    arcfiles = IO.list_file_to_strings(arcfiles)
    arcfilename = filelist_to_path(arcfiles, band, maskname, options)
    fn = "lambda_center_coeffs_{0}.npy".format(arcfilename.rstrip(".fits"))
    arcsols = np.load(fn)

    if len(skysols) != len(arcsols): 
        error("Number of slits in sky (%i) and arcs (%i) is different" % ( len(skysols) , len(arcsols)))
        raise Exception("Number of slits in sky (%i) and arcs (%i) is different" % ( len(skysols) , len(arcsols)))

    pix = np.arange(2048)
    pl.figure(1)

    if len (LROI) == 1: LROI *= len(skysols)
    if len (LROI) != len(skysols): 
        error("Number of solutions is not equal to the LROI vector (%i!=%i)" % ( len(LROI), len(skysols)))
        raise Exception("Number of solutions is not equal to the LROI vector (%i!=%i)" % ( len(LROI), len(skysols)))

    MeanDiffs = []
    for i in xrange(len(skysols)):
        s = skysols[i]
        a = arcsols[i]

        ls = CV.chebval(pix, s['sol_1d'][1])
        la = CV.chebval(pix, a['sol_1d'][1])

        roi = (ls > LROI[i][0]) & (ls < LROI[i][1])
        la -= np.mean( (la-ls)[roi] )

        diff = ls - la
        pl.plot(ls, diff)
        pl.axvline(LROI[i][0], color='blue')
        pl.axvline(LROI[i][1], color='red')

        roi = (ls > 21500) & (ls < 22000)
        MeanDiffs.append(np.mean(ls[roi] - la[roi]))


    if not no_check:
        pl.xlim(19000, 25000)
        pl.ylim(-7,7)
        pl.grid(True)
        pl.title("Close this window to continue")
        pl.xlabel("Sky Wavelength [Angstrom]")
        pl.ylabel("Sky - Arc Wavelength [Angstrom]")

        MeanDiffs =np.array(MeanDiffs)
        info("RMS in 21500 < lambda < 22000 is %2.2f Ang" % (
            np.sqrt(np.mean(MeanDiffs**2))))
        pl.show()

    return LROI


def fit_lambda_interactively(maskname, band, wavenames, options, neon=None,
                             longslit=None,argon=None, extension=None,
                             bypass=False, noninteractive=False, short_exp=False):
    """Fit the one-dimensional wavelength solution to each science slit
    
    Args:
        maskname: The maskname
        band: The spectral band [Y, J, H, K]
        wavenames: List of wavelength standard files
        options: Options dictionary
        neon: Using neon emission lines
        argon: Using argon emission lines
        extension: path to file that contains a well formated fits header
            this should be used only when the detector server fails
            to write the full FITS header
        longslit: Longslit dictionary containing {"yrange": [a,b] and "row_position": YY}
            Note that [a,b] is the range to extract the longslit spectrum over
            row_position is the location to perform the interactive solution over. This
                row should be clean of any contaminatring light

        noninteractive: Bypass the manual fitting and run an autofit routine.
        bypass: Bypass the manual fitting and run an autofit routine.  (This is 
                a duplicate of noninteractive above for backward compatibility).

        """

    ## Set noninteractive mode if either noninteractive or bypass is set
    noninteractive = noninteractive or bypass

    np.seterr(all="ignore")

    
    wavenames = IO.list_file_to_strings(wavenames)
    input_f = filelist_to_path(wavenames, band, maskname, options)
    
    debug("{0} resolves to input files: {1}".format(str(wavenames), str(input_f)))
    info("The wavelength files resolve to input file {0}".format(str(input_f)))
    mfits = IO.readfits(input_f, use_bpm=True)
    (drop, data) = mfits
    (header, drop, bs) = IO.readmosfits(wavenames[0], options, extension=extension)

    mfits = header, data, bs

    name = filelist_to_wavename(wavenames, band, maskname, options)
    fn = "lambda_center_coeffs_{0}.npy".format(name.rstrip(".fits"))

    linelist = pick_linelist(header, neon=neon, argon=argon, short_exp = short_exp)
    
    try: 
        solutions = np.load(fn)
        info( "Solutions loaded from: "+str(fn))
    except IOError: solutions = None

    lamout = np.zeros(shape=(2048, 2048), dtype=np.float32)

    tock = time.time()
    
    outfilename = fn
    if noninteractive is False:
        fig = pl.figure(1,figsize=(16,8))
    else:
        fig = None
    info("Started interactive solution")
    if longslit is not None and longslit['mode'] is "longslit":
        starting_pos = longslit["row_position"]
        info("*** LONGSLIT MODE *** Extract position set to %i" % starting_pos)
    else:
        starting_pos = None
    debug("using line list")
    debug(linelist)
    II = InteractiveSolution(fig, mfits, linelist, options, 1,
        outfilename, solutions=solutions, noninteractive=noninteractive, starting_pos=starting_pos)
    info( "Waiting")

    if noninteractive is False:
        pl.ioff()
        pl.show()
        #pl.draw()
        #pl.show(block=True)

    info("save to: "+str(fn))
    np.save(outfilename, np.array(II.solutions))
    np.save('barset.npy', [II.bs])


def polyfit2d(f, x, y, unc=1.0, orderx=1,ordery=1):
    """Fit a polynomial surface to 2D data, assuming the axes are
    independent of each other.

    Evaluate the fit via:
      f = np.polyval(polyy, y) + np.polyval(polyx, x)

    Usage:
      polyx, polyy, cov = polyfit2d(f, x, y, unc=None, bad=None)

    Input:
      f = an array of values to fit
      x, y = the coordinates for each value of v

    Optional inputs:
      unc = the uncertainties, either one for each value of v or a
            single value for all values; if set to None, then the
            stdev of v is used (an inaccurate measurement of the
            uncertainty, to be sure)
      orderx = the polynomial order (for one dimension) of the fitting
              function
      ordery = same as above
              

    Return value:
     polyx, polyy = the polynomial coefficients (same format as
                    polyfit)
     cov = the covariance matrix

    v1.0.0 Written by Michael S. Kelley, UMD, Mar 2009
    modified 13 aug 2012 by NPK, Caltech
    """

    # the fitting function
    def chi(p, y, x, f, unc, orderx, ordery):
        cy = p[:1+ordery]
        cx = p[1+orderx:]
        model = np.zeros(f.shape) + np.polyval(cy, y) + np.polyval(cx, x)
        chi = (f - model) / unc
        return chi

    # run the fit
    lsq = optimize.leastsq
    guess = np.zeros((orderx + 1)*(ordery+1))
    result = lsq(chi, guess, args=(y, x, f, unc, orderx, ordery), 
        full_output=True)
    fit = result[0]
    cov = result[1]
    cy = fit[:1+ordery]
    cx = fit[1+orderx:]

    return (cx, cy, cov)

def helper_apply_sky_and_arc(positions, coeffs, lambdaMAD):

    global lams, sigs

        

def find_pixel_offset(lam_sky, coeff_arc, LROI):
    '''Find the best pixel offset between the sky wavelength solution and
    Chebyshev fitting polymoials by looping over a pixel range.'''

    roi = (lam_sky > LROI[0]) & (lam_sky < LROI[1])

    dpixs = np.arange(-1, 1, .01)
    RMSs = np.zeros(len(dpixs))
    
    xx = np.arange(2048)
    for rms_cnt in xrange(len(dpixs)):
        dpix = dpixs [rms_cnt]
        dl = lam_sky - CV.chebval(xx - dpix, coeff_arc)
        RMSs[rms_cnt] = np.sqrt(np.mean(dl[roi]**2))

    minix = np.argmin(RMSs)
    best_dpix = dpixs[minix]
    return best_dpix


def apply_lambda_sky_and_arc(maskname, bandname, skynames, arcnames, LROIs,
    options, longslit=None, smooth=True, neon=True, extension=None, short_exp = False):
    
    global lams, sigs
    

    skynames = IO.list_file_to_strings(skynames)
    arcnames = IO.list_file_to_strings(arcnames)

    skyname = filelist_to_wavename(skynames, bandname, maskname,
            options).rstrip(".fits")
    arcname = filelist_to_wavename(arcnames, bandname, maskname,
            options).rstrip(".fits")

    info(str(skyname))
    info(str(arcname))
    drop, data = IO.readfits(skyname+'.fits', use_bpm=True)
    header, drop, bs = IO.readmosfits(skynames[0], options, extension=extension)
    skydata = data.filled(0)

    drop, data = IO.readfits(arcname+'.fits', use_bpm=True)
    header, drop, bs = IO.readmosfits(arcnames[0], options, extension=extension)
    arcdata = data.filled(0)

    slitedges, edgeinfo = IO.load_edges(maskname, bandname, options)
    SkyL = IO.load_lambdadata(skyname, maskname, bandname, options)
    ArcL = IO.load_lambdadata(arcname, maskname, bandname, options)

    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[bandname]

    skylines = pick_linelist(header, short_exp = short_exp)
    arclines = pick_linelist(header, neon=neon)

    # write lambda
    lams = np.zeros((2048, 2048), dtype=np.float32)
    sigs = np.zeros((2048, 2048), dtype=np.float)
    xx = np.arange(2048)
    xsamp = np.array(np.linspace(0, 2047, 10), dtype=np.int)
    
    xypairs = []

    if len(SkyL) != len(ArcL):
        error("Number of lines in Sky file not the same as those in Arc file")
        raise Exception("Number of lines in Sky file not the same as those in Arc file")

    fitpix = np.arange(0,2048,100)
    solutions = []
    for i in xrange(len(SkyL)):
        slp = SkyL[i]["2d"]["positions"].astype(np.int)
        slc = SkyL[i]["2d"]["coeffs"]
        slm = SkyL[i]["2d"]["lambdaMAD"]
        alp = ArcL[i]["2d"]["positions"].astype(np.int)
        alc = ArcL[i]["2d"]["coeffs"]
        alm = ArcL[i]["2d"]["lambdaMAD"]

        info("2d wavelengths: Slit %i/%i" % (i+1, len(SkyL)))

        prev = 0
        dpixels = []
        for j in xrange(len(slp)):

            if (slm[j] < 0.2) and (alm[j] < 0.2):
                
                coeff_sky = slc[j]
                coeff_arc = alc[j]
                lam_sky = CV.chebval(xx, coeff_sky)
                lam_arc = CV.chebval(xx, coeff_arc)

                # minimize pixel solution
                best_dpix = find_pixel_offset(lam_sky, coeff_arc, LROIs[i])
                dpixels.append(best_dpix)

                # Refit the chebyshev
                lambdas = CV.chebval(fitpix - best_dpix, coeff_arc)
                arc_positions = lambdas > np.mean(LROIs[i])
                sky_positions = lambdas <= np.mean(LROIs[i])

                to_fit = CV.chebval(fitpix[sky_positions], coeff_sky)
                to_fit = np.append(to_fit,
                    CV.chebval(fitpix[arc_positions], coeff_arc))
                
                coeffs = CV.chebfit(fitpix, to_fit, options["chebyshev-degree"])


                prev = lams[slp[j],:] = CV.chebval(xx, coeffs)
                prevcoeff = coeffs
                SkyL[i]["2d"]["coeffs"][j,:] = coeffs
            else:
                lams[slp[j],:] = prev

            

        info("Shifted arc by an average of %1.2f pixels" % (np.mean(dpixels)))

        if np.isfinite(np.mean(dpixels)):
            if smooth == True:
                xr = np.arange(len(slp))

                for i in xrange(lams.shape[1]):
                    ff = np.poly1d(Fit.polyfit_clip(xr, lams[slp, i], 3))
                    d = lams[slp,i] - ff(xr)
                    lams[slp, i] = ff(xr)

    info("{0}: writing lambda".format(maskname))

    ### FIX FROM HERE
    fn = "merged_lambda_coeffs_{0}_and_{1}".format(skyname, arcname)
    np.save(fn, SkyL)
    

    header = pf.Header()
    header.set("maskname", maskname)
    header.set("filter", bandname)
    header.set("object", "Wavelengths {0}/{1}".format(maskname, bandname))

    IO.writefits(lams, maskname, "merged_lambda_solution_{0}_and_{1}.fits".format(skyname, arcname), 
            options, overwrite=True, header=header)
                
    info("{0}: rectifying".format(maskname))
    dlam = np.ma.median(np.diff(lams[1024,:]))
    hpp = Filters.hpp[bandname] 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)
    nspec = len(ll_fid)

    rectified = np.zeros((2048, nspec), dtype=np.float32)

    for i in xrange(2048):
        ll = lams[i,:]
        ss = skydata[i,:]

        f = interp1d(ll, ss, bounds_error=False)
        rectified[i,:] = f(ll_fid)

    header.set("object", "Rectified wave FIXME")
    header.set("wat0_001", "system=world")
    header.set("wat1_001", "wtype=linear")
    header.set("wat2_001", "wtype=linear")
    header.set("dispaxis", 1)
    header.set("dclog1", "Transform")
    header.set("dc-flag", 0)
    header.set("ctype1", "AWAV")
    header.set("cunit1", "Angstrom")
    header.set("crval1", ll_fid[0])
    header.set("crval2", 0)
    header.set("crpix1", 1)
    header.set("crpix2", 1)
    header.set("cdelt1", 1)
    header.set("cdelt2", 1)
    header.set("cname1", "angstrom")
    header.set("cname2", "pixel")
    header.set("cd1_1", dlam.item())
    header.set("cd1_2", 0)
    header.set("cd2_1", 0)
    header.set("cd2_2", 1)


    IO.writefits(rectified, maskname, "merged_rectified_{0}_and_{1}.fits".format(skyname, arcname), 
            options, overwrite=True, lossy_compress=True, header=header)
    

def apply_lambda_simple(maskname, bandname, wavenames, options,
        longslit=None, smooth=True, neon=None, short_exp = False):
    """Convert solutions into final output products. This is the function that
    should be used for now."""

    wavenames = IO.list_file_to_strings(wavenames)
    wavename = filelist_to_wavename(wavenames, bandname, maskname,
            options).rstrip(".fits")

    wavepath = filelist_to_path(wavenames, bandname, maskname,
            options)
    drop, data = IO.readfits(wavepath, use_bpm=True)
    header, drop, bs = IO.readmosfits(wavenames[0], options)
    data = data.filled(0)

    slitedges, edgeinfo = IO.load_edges(maskname, bandname, options)
    Ld = IO.load_lambdadata(wavename, maskname, bandname, options)

    if longslit is not None and longslit['mode'] is "longslit":
        info("*** LONGSLIT MODE *** Slit edges set to:")
        info("Bottom: "+str(slitedges[0]["yposs_bot"][0]))
        info("Top:    "+str(slitedges[0]["yposs_top"][0]))
        
    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[bandname]

    lines = pick_linelist(header, neon=neon, short_exp = short_exp)

    # write lambda
    lams = np.zeros((2048, 2048), dtype=np.float32)
    sigs = np.zeros((2048, 2048), dtype=np.float)
    xx = np.arange(2048)
    xsamp = np.array(np.linspace(0, 2047, 10), dtype=np.int)
    
    xypairs = []

    xs = []
    ys = []
    zs = []
    for i in xrange(len(Ld)):
        lp = Ld[i]["2d"]["positions"].astype(np.int)
        lc = Ld[i]["2d"]["coeffs"]
        lm = Ld[i]["2d"]["lambdaMAD"]
        info("Creating 2d wavelength map: Slit %i/%i" % (i+1, len(Ld)))

        prev = 0
        for j in xrange(len(lp)):
            sigs[lp[j],:] = lm[j]

            if lm[j] < 0.18:
                prev = lams[lp[j],:] = CV.chebval(xx, lc[j])
                prevcoeff = lc[j]
                xs.extend(np.ones(len(xsamp)) * lp[j])
                ys.extend(xsamp)
                zs.extend(lams[lp[j], xsamp])
            else:
                lams[lp[j],:] = prev

        if smooth == True:
            xr = np.arange(len(lp))
            for k in xrange(lams.shape[1]):
                ff = np.poly1d(Fit.polyfit_clip(xr, lams[lp, k], 3))
                d = lams[lp,k] - ff(xr)
                lams[lp, k] = ff(xr)


#         if False == True:
#             xs,ys,zs = map(np.array, [xs,ys,zs])
#             info("smoothing")
# 
#             polyx, polyy, cov = polyfit2d(np.array(zs,dtype=np.double),
#                     np.array(ys, dtype=np.double),
#                     np.array(xs, dtype=np.double),
#                     orderx=3,ordery=3)
# 
#             xx, yy = np.array(np.meshgrid(np.arange(2048), lp),
#                     dtype=np.double)
# 
#             M = lams[lp,:] = np.polyval(polyy, yy) + np.polyval(polyx, xx)


    info("writing {} for {}".format("lambda_solution_{0}.fits".format(wavename), maskname))
    header = pf.Header()
    header.set("maskname", maskname)
    header.set("filter", bandname)
    header.set("object", "Wavelengths {0}/{1}".format(maskname, bandname))

    wavename = wavename.rstrip(".fits")
    IO.writefits(lams, maskname, "lambda_solution_{0}.fits".format(wavename), 
            options, overwrite=True, header=header)
                

    info("writing {} for {}".format("sigs_solution_{0}.fits".format(wavename), maskname))
    header.set("object", "Sigmas {0}/{1}".format(maskname, bandname))
    IO.writefits(sigs, maskname, "sigs_solution_{0}.fits".format(wavename), 
            options, overwrite=True, header=header, lossy_compress=True)

    info("writing {} for {}".format("rectified_{0}.fits".format(wavename), maskname))
    dlam = 0
    central_line = 1024
    step = 0
    while dlam==0:
        line = central_line+(10*step)
        dlam = np.ma.median(np.diff(lams[line,:]))
        if dlam==0:
            line = central_line-(10*step)
            dlam = np.ma.median(np.diff(lams[line,:]))
        step=step+1

    # if a masked array, as returned by numpy 1.9, then get the internal representation
    if type(dlam) == np.ma.MaskedArray:
        dlam = dlam.item()
        
    info("Non-empty line found at pixel "+str(line))
    hpp = Filters.hpp[bandname] 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)
    nspec = len(ll_fid)

    rectified = np.zeros((2048, nspec), dtype=np.float32)

    for i in xrange(2048):
        ll = lams[i,:]
        ss = data[i,:]

        f = interp1d(ll, ss, bounds_error=False)
        rectified[i,:] = f(ll_fid)


    header.set("object", "Rectified wave FIXME")
    header.set("wat0_001", "system=world")
    header.set("wat1_001", "wtype=linear")
    header.set("wat2_001", "wtype=linear")
    header.set("dispaxis", 1)
    header.set("dclog1", "Transform")
    header.set("dc-flag", 0)
    header.set("ctype1", "AWAV")
    header.set("cunit1", "Angstrom")
    header.set("crval1", ll_fid[0])
    header.set("crval2", 0)
    header.set("crpix1", 1)
    header.set("crpix2", 1)
    header.set("cdelt1", 1)
    header.set("cdelt2", 1)
    header.set("cname1", "angstrom")
    header.set("cname2", "pixel")
    header.set("cd1_1", dlam)
    header.set("cd1_2", 0)
    header.set("cd2_1", 0)
    header.set("cd2_2", 1)


    IO.writefits(rectified, maskname, "rectified_{0}.fits".format(wavename), 
            options, overwrite=True, lossy_compress=True, header=header)
    
#
# Fitting Methods
#   

# Physical models for instrument
def param_guess_functions(band):
    """Parameters determined from experimentation with cooldown 9 data"""

    alpha_pixel = np.poly1d([-8.412e-16, 3.507e-12, -3.593e-9, 
        6.303e-9, 0.9963]) 

    # Note that these numbers were tweaked by hand by npk on 28 apr
    # they are not reliable. this function should change dramatically.
    if band == 'Y' or band == 'J':
        sinbeta_position = np.poly1d([0.0239, 36.2])
        sinbeta_pixel = np.poly1d([-2.578e-7, 0.00054, -0.2365])
        gamma_pixel = np.poly1d([1.023e-25, -4.313e-22, 7.668e-17, 6.48e-13])
    elif band == 'H' or band == 'K':
        sinbeta_position = np.poly1d([2.331e-2, 38.24])
        sinbeta_pixel = np.poly1d([-2.664e-7, 5.534e-4, -1.992e-1])
        gamma_pixel = np.poly1d([1.033e-25, -4.36e-22, 4.902e-19, -8.021e-17,
            6.654e-13]) 

    delta_pixel = np.poly1d([-1.462e-11, 6.186e-8, -5.152e-5, -0.0396,
        1193])  - 50

    return [alpha_pixel, sinbeta_position, sinbeta_pixel, 
            gamma_pixel, delta_pixel]
    
def dlambda_model(p):
    """Returns an approximate dlambda/dpixel """
    x = 1024
    order = p[4]
    y = p[5]
    (alpha, sinbeta, gamma, delta) = p[0:4]
    sinbeta = np.radians(sinbeta)
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length

    costerm = np.cos(scale * (y-1024))

    return scale/(order/d) * sinbeta / costerm

def pick_linelist(header, neon=False, argon=False, short_exp = False):
    band = header["filter"]
    
    # following linelinests are produced by ccs and can be found in the iraf
    # data bases on ramekin:
    # /scr2/mosfire/mosfire_sky_lines/database

    if band == 'Y':
        lines = np.array([
              9793.6294 , 9874.84889 , 9897.54143 , 9917.43821 , 10015.6207 ,
             10028.0978 , 10046.7027 , 10085.1622 , 10106.4478 , 10126.8684 ,
              10174.623 , 10192.4683 , 10213.6107 , 10289.3707 , 10298.7496 ,
             10312.3406 , 10350.3153 , 10375.6394 , 10399.0957 , 10421.1394 ,
             10453.2888 , 10471.829 , 10512.1022 , 10527.7948 , 10575.5123 ,
             10588.6942 , 10731.6768 , 10753.9758 , 10774.9474 , 10834.1592 ,
             10844.6328 , 10859.5264 , 10898.7224 , 10926.3765 , 10951.2749 ,
             10975.3784 , 11029.8517 , 11072.4773 , 11090.083 , 11140.9467 ,
             11156.0366 , ])


    if band == 'H':
        lines = np.array([
                14605.0225 , 14664.9975 , 14698.7767 , 14740.3346 , 14783.7537 ,
                14833.029 , 14864.3219 , 14887.5334 , 14931.8767 , 15055.3754 ,
                15088.2599 , 15187.1554 , 15240.922 , 15287.7652 ,
                15332.3843 , 15395.3014 , 15432.1242 , 15570.0593 , 15597.6252 ,
                15631.4697 , 15655.3049 , 15702.5101 , 15833.0432 , 15848.0556 ,
                15869.3672 , 15972.6151 , 16030.8077 , 16079.6529 , 16128.6053 ,
                16194.6497 , 16235.3623 , 16317.0572 , 16351.2684 , 16388.4977 ,
                16442.2868 , 16477.849 , 16502.395 , 16553.6288 , 16610.807 ,
                16692.2366 , 16708.8296 , 16732.6568 , 16840.538 , 16903.7002 ,
                16955.0726 , 17008.6989 , 17078.3519 , 17123.5694 , 17210.579 ,
                17248.5646 , 17282.8514 , 17330.8089 , 17386.0403 , 17427.0418 ,
                17449.9205 , 17505.7497 , 17653.0464 , 17671.843 , 17698.7879 ,
                17811.3826 , 17880.341  , 17993.9600 , 18067.9500 ])


    if band == 'J':
        # Removed: 12589.2998 12782.9052 12834.5202
        lines = np.array([
            11538.7582 , 11591.7013 , 11627.8446 , 11650.7735 , 11696.3379 ,
            11716.2294 , 11788.0779 , 11866.4924 , 11988.5382 , 12007.0419 ,
            12030.7863 , 12122.4957 , 12135.8356 , 12154.9582 , 12196.3557 ,
            12229.2777 , 12257.7632 , 12286.964 , 12325.9549 , 12351.5321 ,
            12400.8893 , 12423.349 , 12482.8503 , 12502.43 , 
            12905.5773 , 12921.1364 , 12943.1311 ,
            12985.5595 , 13021.6447 , 13052.818 , 13085.2604 , 13127.8037 ,
            13156.9911 , 13210.6977 , 13236.5414 , 13301.9624 , 13324.3509 ,
             13421.579])

    if (band == 'K') and short_exp:
        # remove: 19518.4784 , 19593.2626 ,  19618.5719 ,19678.046 ,19839.7764 ,20193.1799 ,20499.237 ,21279.1406 ,21580.5093 ,21711.1235 , 21873.507 ,22460.4183 ,22690.1765 ,22985.9156,23914.55, 24041.62,22742.1907 
        lines = np.array([
        19642.4493 , 
        19701.6455 , 19771.9063 , 
        20008.0235 ,  20275.9409 , 20339.697 , 20412.7192 ,
        20563.6072 , 20729.032 , 20860.2122 , 20909.5976 ,
        21176.5323 , 21249.5368 ,  21507.1875 , 21537.4185 ,
        21802.2757 ,  21955.6857 ,
        22125.4484 , 22312.8204 , 22517.9267
        ])
        print("using short exposure line list")
    elif band == 'K':
        #drop: 19751.3895, 19736.4099, 21711.1235, 22
        lines = np.array([
        19518.4784 , 19593.2626 , 19618.5719 , 19642.4493 , 19678.046 ,
        19701.6455 , 19771.9063 , 19839.7764 ,
        20008.0235 , 20193.1799 , 20275.9409 , 20339.697 , 20412.7192 ,
         20499.237 , 20563.6072 , 20729.032 , 20860.2122 , 20909.5976 ,
        21176.5323 , 21249.5368 , 21279.1406 , 21507.1875 , 21537.4185 ,
        21580.5093 , 21711.1235 , 21802.2757 , 21873.507 , 21955.6857 ,
        22125.4484 , 22312.8204 , 22460.4183 , 22517.9267 , 22690.1765 ,
        22742.1907 , 22985.9156, 23914.55, 24041.62])

    if neon:
        # http://www2.keck.hawaii.edu/inst/mosfire/data/MosfireArcs/mosfire_Ne_vac.list
        # Trimmed using PDF of id'd lines
        info("Picking Neon's arc line list")
        if band == 'Y':
            lines = np.array([
                9668.071,
                10298.238,
                10565.303,
                10801.001,
                10847.448,
                11146.072,
                11180.585,
                11393.552,
                11412.257])

        if band == 'J':
            lines = np.array([
                11393.552,
                11412.257,
                11525.900,
                11539.503,
                11617.260,
                11792.271,
                11988.194,
                12069.636,
                12462.799,
                12692.674,
                12915.546])

        if band == 'H':
            lines = np.array([
                14933.886,
                14990.415,
                15144.236,
                15195.083,
                15234.877,
                15352.384,
                15411.803,
                15608.478,
                16027.147,
                16272.797,
                16049.737,
                16479.254,
                16793.378,
                16866.255,
                17166.622])

        if band == 'K':
            lines = np.array([
                #19579.094 ,
                19582.455 ,
                20355.771 ,
                21047.013 ,
                21714.039 ,
                22434.265 ,
                22536.528 ,
                22667.971 ,
                23106.784 ,
                23379.343 ,
                23571.764 ,
                23642.934 ,
                #23918.541 ,
                24168.025 ,
                #24256.224 ,
                #24371.661 ,
                #24378.260 ,
                #24390.011 ,
                #24454.531 ,
                #24459.775 ,
                #24466.068 ,
                #24471.606 ,
                ])

    if argon:
        info("Picking Argon's arc line list")
        if band == 'Y':
            lines = np.array([
                9660.43,
                9787.18,
                10054.81,
                10472.92,
                10480.90,
                10676.49,
                10684.70,
                10883.94,
                11081.90, ])
        if band == 'J':
            lines = np.array([
                11491.25,
                11671.90,
                11722.70,
                11946.55,
                12029.94,
                12115.64,
                12143.06,
                12346.77,
                12406.22,
                12442.73,
                12491.08,
                12736.90,
                12806.24,
                12960.20,
                13011.82,
                13217.61,
                13276.27,
                13316.85,
                13370.77,
                13410.26,
                13507.88,
                13626.38 ])
        if band == 'H':
            lines = np.array([
                14654.35,
                14743.17,
                15050.62,
                15176.84,
                15306.07,
                15333.54,
                15406.85,
                15904.03,
                15993.86,
                16184.44,
                16441.44,
                16524.38,
                16744.65,
                16945.21,
                17449.67,
                17919.61])
        if band == 'K':
            lines = np.array([
                19822.91,
                19971.18,
                20322.56,
                20574.43,
                20621.86,
                20739.22,
                20816.72,
                20991.84,
                21338.71,
                21540.09,
                22045.58,
                22083.21,
                23139.52,
                23851.54,])

    lines = np.array(lines)
    return np.sort(lines)


def guess_wavelength_solution(slitno, header, bs):
    """Given a slit number guess the coefficient values
    return [order, y0, alpha, sinbeta, gamma, delta]
    """

    band = header['filter'].rstrip()
    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]

    y0 = bs.csu_slit_to_pixel(slitno)
    csupos_mm = bs.csu_slit_center(slitno)


    [alpha_pixel, sinbeta_position, sinbeta_pixel, gamma_pixel, 
            delta_pixel] = param_guess_functions(band)

    retv = [alpha_pixel(y0),
            sinbeta_position(csupos_mm) + sinbeta_pixel(y0),
            gamma_pixel(y0),
            delta_pixel(y0),
            order,
            y0, 
            csupos_mm]

    return retv

def refine_wavelength_guess(wave,spec,linelist):
    """Do a cross correlation with the sky lines to get a better 
    guess of the wavelength solution.                            

    INPUTS:                                                      
    ---------                                                    
    wave - wavelength array (needs to be same length as spec)    
    spec - spectrum flux correponding to those wavelengths       
    linelist - a list of line centroids of the reference         

    HISTORY:                                                     
    --------                                                     
    2013-06-27 - T. Do                                           
    """

    # find what is the average peak height to construct the reference
    # spectrum template                                                                                                                          
    peaks = signal.find_peaks_cwt(spec,np.array([2.0,4.0]),noise_perc=5.0,min_snr =1.5)
    avePeak = np.mean(spec[peaks])

    refSpec = np.zeros(len(spec))

    for i in np.arange(len(linelist)):
        refSpec = refSpec + Fit.gaussian([avePeak,linelist[i],2.0,0,0],wave)

    corr = signal.correlate(spec,refSpec,mode='same')
    lags = np.arange(len(spec))-len(spec)/2

    # peak velocity corresponding to the pixel peak                                                                                                  
    peakInd = np.argmax(corr)
    peakLag = lags[peakInd]

    # return the number of pixels that need to be shifted                                                                                             
    return -peakLag

def plot_mask_solution_ds9(fname, maskname, options):
    """makes a ds9 region file guessing the wavelength solution"""


    (header, data, bs) = IO.readmosfits(fname, options)

    linelist = pick_linelist(header)

    ds9 = """# Region file format: DS9 version 4.1
global color=red 
"""
    pix = np.arange(2048)
    colors = ["red", "blue"]
    cidx = 0


    for i in xrange(1,len(bs.ssl)+1):
        slits = bs.scislit_to_csuslit(i)

        info("Guessing: "+str(slits))
        
        cidx = (cidx + 1) % 2
        color = colors[cidx]
        for slitno in slits:
            guess = guess_wavelength_solution(slitno, header, bs)
            ll = wavelength_model(guess, pix)

            if bs.is_alignment_slitno(slitno): color = 'green'
            
            for line in linelist:
                x = np.argmin(np.abs(ll - line))

                ds9 += "circle(%f, %f, 1) # color=%s text={}\n" % (x, guess[5],
                        color)

    path = './'


    fname = fname.rstrip(".fits")
    fn = os.path.join(path, maskname, ("guess_waves_%s.reg" % fname))
    try: os.remove(fn)
    except: pass

    try:
        f = open(fn, 'w')
        f.write(ds9)
        f.close()
    except:
        error("Could not write %s" % fn)


def estimate_half_power_points(slitno, header, bs):
    """This helper function is used to determine the filter half-power points.
    This function is primarily used by the flat-field code to determine the 
    on order regions of an image.  """

    band = header['filter'].rstrip()
    parguess = guess_wavelength_solution(slitno, header, bs)
    pix = np.arange(2048.)
    ll = wavelength_model(parguess, pix)

    hpp = Filters.hpp[band]
    return [ np.argmin(np.abs(ll-hpp[0])), np.argmin(np.abs(ll-hpp[1])) ]



def xcor_known_lines(lines, ll, spec, spec0, options):
    """
    lines[N]: list of lines in wavelength units
    ll[2048]: lambda vector
    spec[2048]: spectrum vector (as function of lambda)
    options: wavelength options
    """
    inf = np.inf
    dxs = []
    sigs = []

    pix = np.arange(2048.)

    for lam in lines:
        f = options["fractional-wavelength-search"]
        roi = np.where((f*lam < ll) & (ll < lam/f))[0]



        if not roi.any():
            dxs.append(inf)
            sigs.append(inf)
            continue
        
        lags = np.arange(-len(roi)/2, len(roi)/2)
        cors = Fit.xcor(spec[roi], spec0[roi], lags)

        fit = Fit.mpfitpeak(lags, cors)

        if (fit.perror is None) or (fit.status < 0):
            dxs.append(inf)
            sigs.append(inf)
            continue

        dxs.append(fit.params[1])
        sigs.append(fit.params[2])

    return map(np.array, [dxs, sigs])


def find_known_lines(lines, ll, spec, options):
    """
    lines[N]: list of lines in wavelength units
    ll[2048]: lambda vector
    spec[2048]: spectrum vector (as function of lambda)
    options: wavelength options
    """
    inf = np.inf
    xs = []
    sxs = []
    sigmas = []

    pix = np.arange(len(spec))

    for lam in lines:
        f = options["fractional-wavelength-search"]
        roi = (f*lam < ll) & (ll < lam/f)

        if not roi.any():
            xs.append(0.0)
            sxs.append(inf)
            continue

        istd = 1/np.sqrt(np.abs(spec[roi].data))

        lsf = Fit.mpfitpeak(pix[roi], spec[roi].data,
                error=istd)

        if (lsf.perror is None) or (lsf.status < 0):
            xs.append(0.0)
            sxs.append(inf)
            continue

        mnpix = np.min(pix[roi])
        mxpix = np.max(pix[roi])

        if (mnpix + 4) > lsf.params[1] < (mxpix-4):
            xs.append(0.)
            sxs.append(inf)
            continue

        if mnpix < 7:
            xs.append(0.0)
            sxs.append(inf)
            continue

        if mxpix > 2040:
            xs.append(0.0)
            sxs.append(inf)
            continue
        
        xs.append(lsf.params[1])
        sxs.append(lsf.perror[1])
        sigmas.append(lsf.params[2])
        
    return map(np.array, [xs, sxs, sigmas])

def fit_chebyshev_to_lines(xs, sxs, lines, options):
    """Fit a chebyshev function to the best fit determined lines.
    Note the best fits may fail and the first part of this function culls
    bad fits, while the second part of the function has all the chebyshev
    action.  For reference the numpy.polynomial.chebyshev package is imported
    as CV """

    ok = np.isfinite(sxs)
    L = len(xs[ok])
    badfit = np.zeros(options["chebyshev-degree"]+1)
    baddelt = np.ones(L) * 9999.0

    if L < 6:
        return [baddelt, badfit, lines[ok]]

    if np.median(lines) < 1000:
        error("Units fed to this function are likely in micron but "
                "should be in angstrom")
        raise Exception("Units fed to this function are likely in micron but "
                "should be in angstrom")

    cfit = CV.chebfit(xs[ok], lines[ok], options["chebyshev-degree"])
    delt = CV.chebval(xs[ok], cfit) - lines[ok]

    if cfit is None:
        return [baddelt, badfit, lines[ok]]

    return [np.array(delt), np.array(cfit), np.array(lines[ok])]



def fit_model_to_lines(xs, sxs, lines, parguess, options, fixed):

    ok = np.isfinite(sxs)

    if len(np.where(ok)[0]) < 3:
        return [[np.inf], parguess, None]

    slambda = sxs * dlambda_model(parguess)

    parinfo = [
        {'fixed': 0, 'value': parguess[0], 'parname': 'alpha', 'step': 1e-5,
            'limited': [0,0], 'limits': [0,0]},
        {'fixed': 0, 'value': parguess[1], 'parname': 'sinbeta', 
            'step': 1e-7, 'limited': [0,0], 'limits': [0, 0]},
        {'fixed': fixed, 'value': parguess[2], 'parname': 'gamma','step': 1e-12,
            'limited': [1,1], 'limits': [-50e-13, 50e-13]},
        {'fixed': fixed, 'value': parguess[3], 'parname': 'delta', 'step': 1,
            'limited': [1,1], 'limits': [0, 2048]},
        {'fixed': 1, 'value': parguess[4], 'parname': 'order', 
            'limited': [0,0], 'limits': [3, 7]},
        {'fixed': 1, 'value': parguess[5], 'parname': 'Y',
            'limited': [0,0], 'limits': [0, 2048]}
    ]

    merit_function = Fit.mpfit_residuals(wavelength_model)
    
    lsf = Fit.mpfit_do(merit_function, xs[ok], lines[ok], 
            parinfo, error=slambda[ok])

    delt = np.abs(wavelength_model(lsf.params, xs[ok]) - lines[ok])

    xsOK = xs[ok]
    linesOK = lines[ok]


    return [ delt, lsf.params, lsf.perror]

def guesslims(spec):
    """Guess the spectral limits"""
    f = 1.1

    s = spec.copy()
    s.sort()
    return [-500, s[-10]*f]

def polyfit_err(x,y,order):

    coef = np.polyfit(x,y,order)
    fun = np.poly1d(coef)

    return fun, coef, np.std(y - fun(x))



class InteractiveSolution:

    header = None
    data = None
    bs = None
    parguess = None
    linelist0 = None
    foundlines = None
    options = None
    slitno = None   # the science slit number.
    spec = None
    good_solution = False
    done = False
    starting_pos = None # This overrides the extract position for longslits

    ll = None
    pix = None
    xlim = None
    ylim = None
    MAD = None
    STD = None

    first_time = True

    def __init__(self, fig, mfits, linelist, options, slitno, outfilename, 
            solutions=None, noninteractive=False, starting_pos=None):
        self.header = mfits[0]
        self.data = mfits[1]
        self.bs = mfits[2]
        self.options = options
        self.linelist0 = linelist
        self.slitno = slitno
        self.fig = fig
        self.outfilename = outfilename
        self.starting_pos = starting_pos
        self.noninteractive = noninteractive
        self.pix = np.arange(2048)
        band = self.header["filter"].rstrip()
        self.xrng = Filters.hpp[band][:]
        self.band = band
        self.xrng[0] *= 0.99
        self.xrng[1] /= 0.99
        self.sigma_clip = False
        self.xlim = self.xrng
        if solutions is None:
            self.solutions = range(len(self.bs.ssl))
        else:
            self.solutions = solutions

        if self.noninteractive:
            self.setup()
            self.fit_event(0,0)
            #self.nextobject(0,0)  ### the call to next object is built-in in fit_event

            while self.done is False:
                self.fit_event(0,0)
                self.nextobject(0,0)
        else:
            # follow line prevents window from going full screen when the
            # 'f'it button is pressed.
            pl.rcParams['keymap.fullscreen'] = ''
            self.cid = self.fig.canvas.mpl_connect('key_press_event', self)
            self.setup()
            self.fit_event(0,0)

    
    def setup(self):
        csuslits = self.bs.scislit_to_csuslit(self.slitno)
        try:
            l = len(csuslits)
            if l > 1:
                csuslit = csuslits[l/2]
            else:
                csuslit = csuslits[0]
        except:
            csuslit = csuslits


#         info(str(csuslits)+" "+str(csuslit))
        info('CSU slits {} acting as slit number {}'.format(str(csuslits), str(csuslit)))
        self.linelist = self.linelist0
        if self.starting_pos is None:
            self.extract_pos = self.bs.science_slit_to_pixel(self.slitno)
        else:
            info("LONGSLIT mode: forced longslit center line")
            '''This is used in longslits to handle a forced start position'''
            self.extract_pos = self.starting_pos

        info("Extracting at %i " % self.extract_pos)

        S = self.solutions[self.slitno-1]
        if type(S) is not int: # previously setup
            self.MAD = S["MAD"]
            self.STD = S["STD"]
            self.linelist = S["linelist"]
            self.foundlines = S["foundlines"]
            self.foundlinesig = S["foundlinesig"]
            self.extract_pos = S["extract_pos"]
            self.cfit = S["sol_1d"][1]
        else:
            self.MAD = self.STD = self.foundlines = self.linesig = None
            tx = np.arange(0,2048,100)
            parguess = guess_wavelength_solution(csuslit, self.header,
                self.bs)
            ll = wavelength_model(parguess, tx)

            refineGuess= True
            if refineGuess:
                # do an inital fit 
                cfit = CV.chebfit(tx, ll, self.options["chebyshev-degree"])
                inputWave = CV.chebval(self.pix, cfit)

                self.spec = \
                    np.ma.mean(self.data[self.extract_pos-1:self.extract_pos+1, :],
                        axis=0) # axis = 0 is spatial direction

                # refine the wavelength guess using the reference lines
                refShift = refine_wavelength_guess(inputWave,self.spec,self.linelist)

                # add the shift to the wavelength solution                           
                ll = wavelength_model(parguess, tx+refShift)

            self.cfit = CV.chebfit(tx, ll, self.options["chebyshev-degree"])

# I don't understand why there is this offset added, but I am removing it 
#    because the refineGuess method results in more accurate positioning.      
#            if self.band != 'J':
#                self.cfit[0] -= 16.0

        self.spec = \
            np.ma.mean(self.data[self.extract_pos-1:self.extract_pos+1, :],
                axis=0) # axis = 0 is spatial direction

        self.ll = CV.chebval(self.pix, self.cfit)

        if self.noninteractive:
            pass
        else:
            info("Launching graphics display.    ")
            self.redraw()

    def toggle_noninteractive(self,x,y):
        info("############ NON INTERACTIVE MODE ENABLED ###########")
        info("# From now on, the fit will proceed automatically   #")
        info("#####################################################")
        self.noninteractive = True
        self.quit(0,0)
        self.nextobject(0,0)
        

    def draw_found_lines(self):
        pl.subplot(2,1,1)

        pl.grid(True)
        xmin, xmax, ymin, ymax = pl.axis()
        if self.foundlines is not None:
            foundlams = CV.chebval(self.foundlines, self.cfit)
            ok = np.isfinite(self.foundlinesig) 

            for i in xrange(len(self.linelist)):
                if not ok[i]: continue
                D = (foundlams[i] - self.linelist[i])
                pl.axvline(foundlams[i], color='orange', ymax=.75, ymin=.25,
                        linewidth=1.5)
                pl.text(foundlams[i], 1500, "%1.2f" % D, rotation='vertical',
                        size=10)

            pl.subplot(2,1,2)
            pl.xlim(self.xlim)
            pl.grid(True)
            #pl.axhline(0.1)
            #pl.axhline(-0.1)
            pl.axhline(self.STD)
            pl.axhline(-1*self.STD)

            if self.STD < 0.1: fmt = 'go'
            else: fmt = 'bo'
            pl.plot(self.linelist[ok], (foundlams[ok] - self.linelist[ok]), 
                    fmt)
            pl.xlim(self.xlim)



    def draw_done(self):
        if not self.done: 
            return

        mid = np.mean(self.xlim)*.99

        pl.subplot(2,1,1)
        pl.text(mid, 0, 'Done!', size=32, color='red')

    def draw_vertical_line_marks(self):
        pl.subplot(2,1,1)
        xmin, xmax, ymin, ymax = pl.axis()
        i = 0
        for line in self.linelist:
            pl.axvline(line, color='red', linewidth=.5)

            pl.text(line, ymax*.75, "%5.1f" % (line), 
                    rotation='vertical', color='black')

            i = i+1
            fwl = self.options['fractional-wavelength-search']
            pl.plot([line*fwl,line/fwl], [0,0], linewidth=2)

    def redraw(self):
        pl.ion()
        pl.clf()

        pl.subplot(2,1,1)
        pl.subplots_adjust(left=.1,right=.95,bottom=.1,top=.90)
        pl.plot(self.ll, self.spec, linestyle='steps-mid')

        if self.MAD is None:
            pl.title("[%i] Press 'z' to zoom, 'x' to unzoom, 'c' to shift, " 
                    "'f' to fit, 'k' to toggle sigma clipping. 'h' for help" % self.slitno)
        else:
            name = self.bs.ssl[self.slitno-1]["Target_Name"]

            pl.title(u"[%i,%s, p%i] Best fit STD: %0.2f $\AA$, MAD: %0.2f $\AA$: " \
                     % (self.slitno, name, self.extract_pos,self.STD, self.MAD))

        pl.ioff()
        self.draw_vertical_line_marks()
        self.draw_found_lines()
        pl.ion()

        pl.subplot(2,1,1)
        xmin, xmax, ymin, ymax = pl.axis()
        pl.xlim(self.xlim)
        if self.band == 'Y': pl.ylim([-100, 1000])
        else: pl.ylim([-1000, ymax*.8])

        if np.max(self.spec) < 200:
            pl.ylim([-100,500])
        
        self.draw_done()


    def shift(self, x, y):
        """Shift the observed spectrum"""
        theline = np.argmin(np.abs(x - self.linelist))

        delt = x - self.linelist[theline] 
        self.ll -= delt
        self.redraw()

    def drop_point(self, x, y):
        """Drop point nearest in x from set"""
        theline = np.argmin(np.abs(x - self.linelist))
        self.linelist = np.delete(self.linelist, theline)
        if self.foundlines is not None:
            self.foundlines = np.delete(self.foundlines, theline)
            self.foundlinesig = np.delete(self.foundlinesig, theline)

        self.redraw()
        
    def unzoom(self, x, y):
        """Show the full spectrum"""
        self.xlim = self.xrng
        pl.ion()
        pl.subplot(2,1,1) ; pl.xlim(self.xlim)
        pl.subplot(2,1,2) ; pl.xlim(self.xlim)

    def zoom(self, x, y):
        """Zoom/pan the view"""
        self.xlim = [x*.988,x/.988]
        pl.ion()
        pl.subplot(2,1,1) ; pl.xlim(self.xlim)
        pl.subplot(2,1,2) ; pl.xlim(self.xlim)

    def fastforward(self, x, y):
        """Fast forward to next uncalib obj """
        for i in xrange(self.slitno+1, len(self.solutions)):
            if type(self.solutions[i]) is int:
                self.slitno = i-1
                self.setup()
                break


    def nextobject(self, x, y):
        """Go to the next object"""
        self.slitno += 1
        self.done = False
        if self.slitno > len(self.bs.ssl): 
            self.done = True
            self.slitno -= 1
            if self.noninteractive is False:
                self.draw_done()

        info("Saving to: "+str(self.outfilename))
        np.save(self.outfilename, np.array(self.solutions))

        if self.done is False:
            self.setup()
            self.fit_event(0,0)

    def prevobject(self, x, y):
        """Go to the previous object"""
        self.slitno -= 1
        self.done = False
        if self.slitno < 1: 
            warning("first limit")
            self.slitno = 1
        self.setup()

    def quit(self, x, y):
        """Quit and save the results """
        info("Closing figure")
        self.fig.canvas.mpl_disconnect(self.cid)
        pl.close(self.fig)


    def reset(self, x, y):
        """Reset the fitting performed on this object """
        self.MAD = None
        self.solutions[self.slitno-1] = self.slitno
        self.setup()

    def savefig(self, x, y):
        """Save the figure to disk"""
        pass

    def toggle_sigma_clip(self,x,y):
        if self.sigma_clip is True: 
            self.sigma_clip = False
            info("Sigma clipping disabled")
        else:
            self.sigma_clip = True
            info("Sigma clipping enabled")
            self.fit_event(0,0)
    
    def fit_event(self, x, y):
        """Fit Chebyshev polynomial to predicted line locations """

        [xs, sxs, sigmas] = find_known_lines(self.linelist, self.ll,
            self.spec, self.options)
        self.foundlines = xs
        self.foundlinesig = sxs
        mask = (np.isfinite(sxs))
        local_linelist=self.linelist[mask]
        xs = xs[mask]
        sxs = sxs[mask]
        [deltas, cfit, perror] = fit_chebyshev_to_lines(xs, sxs,
            local_linelist, self.options)
        self.cfit = cfit
        self.ll = CV.chebval(self.pix, self.cfit)

        # Calculate current std error
        error = np.std(deltas[np.isfinite(deltas)])
        if self.sigma_clip is True or self.noninteractive:
            # prepare a sigma tolerance (reject values of deltas > tolerance * sigma)
            tolerance = 3
            # if the std error is > 0.10, iteratively reject lines
            while error>0.10:
#                 info("#####################################################")
                warning("Large error detected. Iterating with sigma clipping")
                warning("Current error is "+str(error))
                warning("Current tolerance is "+str(tolerance)+" sigmas")
                warning("Number of lines used for fit: "+str(len(xs)))
                warning("Filtering with rms = "+str(np.std(deltas[np.isfinite(deltas)])))
                mask = (abs(deltas)<tolerance*np.std(deltas[np.isfinite(deltas)]))
                warning("Number of rejected lines: "+str(len(xs)-len(xs[mask])))
                local_linelist = local_linelist[mask]
                xs=xs[mask]
                sxs=sxs[mask]
                info("Fitting again...")
                [deltas, cfit, perror] = fit_chebyshev_to_lines(xs, sxs,
                                                                local_linelist, self.options)
                error = np.std(deltas[np.isfinite(deltas)])
                if error<=0.10:
                    info("The error is now {}.  The error is acceptable, continuing...".format(error))
                else:
                    info("The error is now {}".format(error))
#                 info("#####################################################")
                self.cfit = cfit
                self.ll = CV.chebval(self.pix, self.cfit)
                tolerance = tolerance - 0.2
                self.foundlines = xs
                self.foundlinesig = sxs
                self.linelist = local_linelist
        ok = np.isfinite(deltas)
        self.STD = np.std(deltas[ok])
        self.MAD = np.median(np.abs(deltas[ok]))

        debug("STD: %1.2f MAD: %1.2f" % (self.STD, self.MAD))
        debug(str(self.cfit))


        self.solutions[self.slitno-1] = {"linelist": self.linelist, "MAD":
                self.MAD, "foundlines": self.foundlines, "foundlinesig":
                self.foundlinesig, "sol_1d": [deltas, cfit, sigmas], "STD":
                self.STD, "slitno": self.slitno, "extract_pos":
                self.extract_pos}

        info('Stored slit number: {}'.format(str(self.solutions[self.slitno-1]['slitno'])))

        if self.noninteractive is False:
            self.redraw()
        else:
            self.nextobject(0,0)

    def __call__(self, event):
        kp = event.key
        x = event.xdata
        y = event.ydata

        info( str(kp)+" "+str(x)+" "+str(y))

        actions_mouseless = {".": self.fastforward, "n": self.nextobject, "p":
                self.prevobject, "q": self.quit, "r": self.reset, "f":
                self.fit_event, "k": self.toggle_sigma_clip, "\\": self.fit_event, "b": self.toggle_noninteractive}

        actions = { "c": self.shift, "d": self.drop_point,
                "z": self.zoom, "x": self.unzoom, "s": self.savefig}

        if (kp == 'h') or (kp == '?'):
            info("Commands Desc")
            for key, value in actions.items():
                info("%8s %s" % (key, value.__doc__))
            for key, value in actions_mouseless.items():
                info("%8s %s" % (key, value.__doc__))

        if kp in actions_mouseless:
            actions_mouseless[kp](x, y)

        if x is None: return
        if y is None: return

        if kp in actions:
            actions[kp](x, y)

def fit_wavelength_solution(data, parguess, lines, options, 
        slitno, search_num=145, fixed=False):
    """Tweaks the guessed parameter values and provides 1d lambda solution
    
    """

    pix = np.arange(2048.)
    MAD = np.inf

    y0 = parguess[5]
    spec = np.ma.median(data[y0-1:y0+1, :], 
        axis=0) # axis = 0 is spatial direction

    d = 0.1
    dsinbetas = np.sort(np.abs(np.linspace(-d/2., d/2., search_num)))

    sinbetadirection = 1.0

    iteration = 0

    DRAW = False
    if DRAW:
        pl.ion()
        pl.figure(2, figsize=(16,5))
        pl.xlim([2.03,2.3])

    #print "iter  dsb      MAD"
    for dsinbeta in dsinbetas:
        dsinbeta *= sinbetadirection
        sinbetadirection *= -1


        pars = parguess
        pars[1] = parguess[1] + dsinbeta
        ll = wavelength_model(parguess, pix)
        [xs, sxs, sigmas] = find_known_lines(lines, ll, spec, options)
        [deltas, params, perror] = fit_model_to_lines(xs, sxs, lines, 
                pars, options, fixed)

        if DRAW:
            pl.figure(2)
            pl.xlim([1.94,2.1])
            ll2 = wavelength_model(params, pix)
            pl.plot(ll2, spec)
            for line in lines:
                pl.axvline(line ,color='red')
            pl.draw()

        MAD = np.ma.median(deltas)
        iteration += 1
        #print "%2.2i] %3.0i %1.4f %1.4f" % (slitno, iteration, dsinbeta, MAD)


        if MAD > MADLIMIT:
            continue
        else: 
            #print "%i] found: %3i %+1.5f %3.6f" % (slitno, iteration, dsinbeta, MAD)
            break


    if MAD <= MADLIMIT:
        #print("%3i: %3.5f %4.3f %3.3e %4.1f %1.4f" % (slitno, params[0], params[1],
            #params[2], params[3], MAD))

        return [deltas, params, perror, sigmas]
    else:
        warning("%3i: Could not find parameters" % slitano)
        return [[], parguess, None, []]


def construct_model(slitno):
    """Given a matrix of Chebyshev polynomials describing the wavelength
    solution per pixel, return a 2-d function returning the wavelength
    solution.

    For a set of Chebyshev coefficients c_1 to c_6 there is a strong
    correlation between c_n and c_1 in pixel space. This correlation is
    used to produce the 2-d wavelength solution.
    """

    global coeffs, data, linelist, options

    info("Constructing model on %i" % slitno)

    pix = np.arange(2048)

    cfits   = coeffs[slitno-1]['2d']['coeffs']
    delts   = coeffs[slitno-1]['2d']['delts']
    sds     = coeffs[slitno-1]['2d']['lambdaRMS']
    mads    = coeffs[slitno-1]['2d']['lambdaMAD']
    positions= coeffs[slitno-1]['2d']['positions']


    ok = (sds < 0.2) & (mads < 0.1)
    c0coeff = Fit.polyfit_sigclip(positions[ok], cfits[ok,0], 1, nmad=3)
    c0fun = np.poly1d(c0coeff)

    cfit_coeffs = [c0coeff]
    cfit_funs = [c0fun]

    for i in xrange(1, cfits.shape[1]):
        if i < 3: order = 1
        else: order = 1
        ci_coeff = Fit.polyfit_sigclip(cfits[ok,0], cfits[ok,i], order, nmad=3)
        cfit_coeffs.append(ci_coeff)

        ci_fun = np.poly1d(ci_coeff)
        cfit_funs.append(ci_fun)

    lambdaRMS = []
    lambdaMAD = []

    cpolys = []
    # Check the fits now
    if True:
        for i in xrange(len(positions)):
            pos = positions[i]
            c0 = c0fun(pos)
            cs = [c0]
            cs.extend([f(c0) for f in cfit_funs[1:]])

            ll_now = CV.chebval(pix, cs)
            spec_here = np.ma.median(data[pos-1:pos+1,:], axis=0)

            [xs, sxs, sigmas] = find_known_lines(linelist,
                ll_now, spec_here, options)

            [delt, cfit, lines] = fit_chebyshev_to_lines(xs, sxs,
                linelist, options)

            rms = np.std(delt)
            rmsR = np.median(np.abs(lines/delt))
            if rms  > .2:
                pre = bcolors.FAIL
            elif rms > .1:
                pre = bcolors.OKBLUE
            else:
                pre = bcolors.ENDC


            info(pre + "2d model S%2.2i p%4.4i: residual %1.3f Angstrom " \
                    "RMS / %3.1e lam/dlam / %1.3f Angstrom MAD" % (slitno, pos,
                            rms, rmsR, np.median(np.abs(delt))) + bcolors.ENDC)


            lambdaRMS.append(np.std(delt))
            lambdaMAD.append(np.median(np.abs(delt)))
            cpolys.append(cs)


        """The return function takes a pixel position and returns wavelength
        in angstrom"""

        return {"functions": np.array(cfit_coeffs), "cpolys": np.array(cpolys),
                "RMS": np.array(lambdaRMS), "MAD": np.array(lambdaMAD),
                "positions": positions[ok]}

    else:
        info("No fit stats")

        return {"functions": np.array(cfit_coeffs), "cpolys": [], "RMS":
                np.zeros(len(positions[ok])), "MAD":
                np.zeros(len(positions[ok])), "positions": positions[ok]}


#
# Two dimensional wavelength fitting
#
def fit_outwards_refit(data, bs, sol_1d, lines, options, start, bottom, top,
        slitno, linelist2=None, data2=None, sol_1d2=None):
    '''Fit Chebyshev polynomials across a slit bounded by bottom to top.

    Args:
        data: 2k x 2k data frame
        bs: The barset
        sol_1d: The guess solution taken at the 'start' position
        lines: The list of lines (wavelength in Angstrom)
        options: Options dictionary
        start: The start position (should be where the interactive wavelength
            fit occurred.
        bottom/top: The bottom/top of the slit in pixels

    Optional Arguments:
        linelist2 - second line list used for when both the neon and argon 
                    lines with be fit simultaneously.
        data2     - is the spectra data set that goes along with the argon
                    line list when argon and neon are used simultaneously.
    Returns:
        Npixel length array of dictionaries containing:
        'coeffs': Chebyshev coefficients
        'delts': The offset of the line position to the measured position
            in Angstrom
        'lambdaRMS': The standard deviaiton of delts
        'lambdaMAD': The MAD of delts
        'positions': Line positions in pixel units'

        The wavelength solution at some pixel is computed by

        sol = fit_outwards_refit(...)
        wave_vector = CV.chebval(arange(2048), sol[pixelnum]['coeffs'])

        where wave_vector holds the wavelength in angstroms of a 
            row.
        
    Modification history: 
        2014 June 17   MK  - Added code to handle two line lists simultaneously
                             This will allow observers to fit both the Argon 
                             and the Neon lines. Added optional parameters data2 
                             linelist2.

    '''
    lags = np.arange(-30,30)
    pix = np.arange(2048.)
    linelist = lines

    def fit_parameters(yhere):
        """
        Return chebyshev fit to a pixel column 
        2014 June 17 MK- Added a second set a variables to indicate that there
        are two sets of lines we want to fit. This is used 
        for the arc line data. We want to fit both Ne and
        Argon simultaneously.
        
        """        

        cfit = sol_1d[1]
        spec_here = np.ma.median(data[int(yhere)-2:int(yhere)+2, :], axis=0)
        shift = Fit.xcor_peak(spec_here, spec0, lags)
        ll_here = CV.chebval(pix - shift, cfit)
        [xs, sxs, sigmas] = find_known_lines(linelist,
                                             ll_here, spec_here, options)

        if data2 is not None:
            cfit2 = sol_1d2[1]
            spec_here2 = np.ma.median(data2[yhere-2:yhere+2, :], axis=0)
            shift2 = Fit.xcor_peak(spec_here2, spec2, lags)
            ll_here2 = CV.chebval(pix - shift2, cfit2)

            [xs2, sxs2, sigmas2] = find_known_lines(linelist2,
                                  ll_here2, spec_here2, options)

        """ Fit a chebyshev to the measured positions """
        if data2 is not None:
            "fit both line lists"
            """combine the line lists"""
            clinelist= np.concatenate([linelist,linelist2])
            cxs  = np.concatenate([xs, xs2])
            csxs = np.concatenate([sxs, sxs2])
             
            """combine the measured xs and sxs arrays that have the measured
            line positions"""
            [delt, cfit, lines] = fit_chebyshev_to_lines(cxs, csxs,
                                                         clinelist, options)
        else:
            [delt, cfit, lines] = fit_chebyshev_to_lines(xs, sxs,
                                                         linelist, options)

        #if np.std(delt) < .01: pdb.set_trace()
        debug("resid ang S%2.2i @ p%4.0i: %1.2f rms %1.2f mad [shift%2.0f]" % \
                (slitno+1, yhere, np.std(delt), np.median(np.abs(delt)),
                    shift))

        return cfit, delt

    def sweep(positions):
        ret = []
        cfits = []
        sds = []
        mads = []

        for position in positions:
            cfit, delt = fit_parameters(position)

            cfits.append(cfit)
            sds.append(np.std(delt))
            mads.append(np.median(np.abs(delt)))


        cfits, sds, mads = map(np.array, [cfits, sds, mads])
        #model = construct_model(cfits, positions, sds)

        assert(len(positions) == len(cfits))
        return {'coeffs': cfits, 'delts': delt, 'lambdaRMS':
                sds, 'lambdaMAD': mads, "positions": np.array(positions)}

    """ Start of main section of fit_outwards """
    pix = np.arange(2048.)

    positions = np.concatenate((np.arange(start, top, 1), 
        np.arange(start-1,bottom,-1)))
    positions = np.arange(bottom, top, 1)
    info("Computing 0 spectrum at %i" % start)
    spec0 = np.ma.median(data[start-1:start+1, :], axis=0)
    if data2 is not None:
            spec2 = np.ma.median(data2[start-1:start+1, :], axis=0)
    params = sweep(positions)

    return params

class NoSuchFit(Exception):
    pass

def fit_to_coefficients(fits, pos, slitno=None):
    """Given a fit structure from fit_outwards_refit and a pixel y position
    return the coefficients of a Chebyshev Polynomial"""


    if slitno is None:
        slitno = 1
        found = False
        for fit in fits:
            if type(fit) is int: continue
            fitpos = fit["positions"]
            mn = min(fitpos) ; mx = max(fitpos)
            if (pos >= mn) and (pos <= mx):
                found = True
                break
            slitno += 1
        if not found:
            warning("Position %i does not have a fitted wavelength " \
            " solution" % pos)
            raise NoSuchFit("Position %i does not have a fitted wavelength " \
            " solution" % pos)
            return np.zeros(5)

        fit = fits[slitno-1]
        
    else:
        fit = fits[slitno-1]
        fitpos = fit["positions"]
        mn = min(fitpos) ; mx = max(fitpos)
        if not ((pos >= mn) and (pos <= mx)):
            warning("Slitno %i has a pixel range of %i-%i but " \
                    "position %i was requested" % (slitno, mn, mx, pos))
            raise Exception("Slitno %i has a pixel range of %i-%i but " \
                    "position %i was requested" % (slitno, mn, mx, pos))
            return np.zeros(5)

    funs = fit["functions"]

    cs = np.zeros(funs.shape[0])
    cs[0] = np.poly1d(funs[0])(pos)
    for i in xrange(1,len(cs)):
        cs[i] = np.poly1d(funs[i])(cs[0])

    return np.array(cs)

    
def fit_to_lambda(fits, pix_lambda, pix_spatial, slitno=None):
    """Given a fit structure from fit_outwards_refit and a pixel x,y position
    return the wavelength value"""

    cs = fit_to_coefficients(fits, pix_spatial, slitno=slitno)

    return CV.chebval(pix_lambda, cs) 


# Model Functions

def wavelength_model(p, x):
    """Returns wavelength [um] as function of pixel (x)
    
    The parameter list, p, contains
    p[0:3] -- alpha, beta, gamma, delta, model parameters
    p[4] -- the grating order.
    p[5] -- the pixel y position on detector [pix]

    x -- the x pixel position (dispersion direction)

    returns wavelength [micron]
    """
    order = p[4]
    y = p[5]
    (alpha, sinbeta, gamma, delta) = p[0:4]
    sinbeta = np.radians(sinbeta)
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length

    costerm = np.cos(scale * (y-1024))

    return (alpha/(order/d) * 1/costerm * \
            (np.sin(scale * (x-1024)) + sinbeta) + \
            gamma * (x - delta)**3)*1e4

def mask_model(p, xs):
    """Fit a continuous smooth function to parameters in the mask.

    parameters:
        linear model is:
        x = xs - p[3]         2         3
        p[0] + p[1] x + p[2] x  + p[3] x  + discontinuity

        p[4:] -- [N] list of discontinuities
        """

    cpix = p[0]
    cpar = p[1]
    radius_pix = p[2]
    radius_par = p[3]
    coeffs = p[4:]

    vals = []
    for i in xrange(len(coeffs)):
        x = np.array(xs[i]) - p[3]
        c = coeffs[i]
        y = p[0] + p[1] * x + p[2] * x*x + c
        vals.extend(y)

    return np.array(vals).ravel() * 1e4


def plot_mask_fits(maskname, fname, options):

    from matplotlib.backends.backend_pdf import PdfPages

    mfits = IO.readmosfits(fname, options)
    header, data, bs = mfits

    band = header['filter'].rstrip()

    fname = fname.rstrip(".fits")
    Lmask = IO.load_lambdamodel(fname, maskname, band, options)
    Lslit = IO.load_lambdadata(fname, maskname, band, options)

    assert(len(Lmask) == len(Lslit))
    outname = os.path.join(path, "mask_fit_%s.pdf" % fname)
    pp = PdfPages(outname)

    for i in xrange(len(Lmask)):
        print i
        ll = Lmask[i]
        ls = Lslit[i]
        coeffs  = ls['2d']['coeffs']
        sds     = ls['2d']['lambdaRMS']
        posc    = ls['2d']['positions']
        fits    = ll['functions']
        pos     = ll['positions']



        #N = fits.shape[1]
        N = 6
        ny = 2.0
        nx = np.ceil(N/ny)

        pl.clf()
        
        c0 = np.poly1d(fits[0])
        c0s = c0(pos)

        for i in xrange(1,N):
            pl.subplot(int(nx),int(ny),i)
            f = np.poly1d(fits[i])
            pl.plot(pos, f(c0s), color='orange')
            ylim = pl.ylim()
            pl.scatter(posc, coeffs[:, i])
            pl.ylim(ylim)

            #newfit = np.poly1d(Fit.polyfit_clip(c0(posc), coeffs[:,i], 0))
            #pl.plot(pos, newfit(c0(pos)))


        pp.savefig()

    pp.close()



def plot_sky_spectra(maskname, fname, options, short_exp = False):

    from matplotlib.backends.backend_pdf import PdfPages

    fp = os.path.join(options['indir'], fname)
    mfits = IO.readmosfits(fp)
    header, data, bs = mfits

    
    band = header["filter"].rstrip()
    fname = fname.rstrip(".fits")
    solutions = IO.load_lambdadata(fname, maskname, band, options)

    outname = os.path.join(path, "sky_spectra_%s.pdf" % fname)

    pp = PdfPages(outname)
    band = header['filter'].rstrip()

    # determine region to cutoff spectra for xlims
    linelist = pick_linelist(header, short_exp = short_exp)
    hpps = Filters.hpp[band]

    # Pick top 95% of flux for ylims
    sdata = np.sort(data, None)
    ymax = sdata[-15000]


    pix = np.arange(2048)
    for solution in solutions:
        slitno = solution["slitno"]
        parameters = solution["center_sol"][0]
        info("Slit: {0}".format(slitno))

        parguess = guess_wavelength_solution(slitno, header, bs)
        y0 = parguess[-2]

        ll = wavelength_model(parameters, pix)
        measured = data[y0, :]

        pl.clf()
        pl.title("Slit {0}".format(solution["slitno"]))
        pl.plot(ll, measured, linewidth=.2)
        pl.xlim(hpps)
        pl.ylim(-30, ymax)

        for line in linelist:
            pl.axvline(line, color='red', linewidth=.1)

        pp.savefig()

    pp.close()



def plot_data_quality(maskname, fname, options):

    from matplotlib.backends.backend_pdf import PdfPages
    path = os.path.join(options["indir"])
    if not os.path.exists(path):
        error("Output directory '%s' does not exist. This " 
            "directory should exist." % path)
        raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

    fp = os.path.join(path, fname)
    mfits = IO.readmosfits(fp)
    header, data, bs = mfits

    
    fname = fname.rstrip(".fits")
    path = './'
    solname = os.path.join(path, "lambda_coeffs_%s.npy" % fname)
    solutions = np.load(solname)
    solname = os.path.join(path, "mask_solution_%s.npy" % fname)
    masksol = np.load(solname)[0]


    outname = os.path.join(path, "wavelength_fits_%s.pdf" % fname)

    pp = PdfPages(outname)

    filter_fun = (lambda x:
            (x[1] is not None) and
            (x[1][0] < 1e-5) and
            (x[2] < .2) and
            (x[3] == True))

    all_pix = []
    all_alphas = []
    all_betas = []
    all_gammas = []
    all_deltas = []
    for solution in solutions:
        sol_2d = solution["2d"]
        info("Slit: {0}".format(solution["slitno"]))
        ff = filter(filter_fun, sol_2d)
        ar = np.array(map(lambda x: x[0], ff))

        if len(ar) == 0: continue

        pixels = ar[:,5]

        alphas = ar[:,0]
        betas = ar[:,1]
        gammas = ar[:,2]
        deltas = ar[:,3]
        sds = ar[:,4]


        all_pix.extend(pixels)
        all_alphas.extend(alphas)
        all_betas.extend(betas)
        all_gammas.extend(gammas)
        all_deltas.extend(deltas)

        alphamodel = np.poly1d(np.polyfit(pixels, alphas, 1))
        betamodel  = np.poly1d(np.polyfit(pixels, betas, 1))
        gammamodel = np.poly1d(np.polyfit(pixels, gammas, 1))
        deltamodel = np.poly1d(np.polyfit(pixels, deltas, 1))

        info("Scatters: {0:3.5} {1:3.5} {2:3.5} {3:3.5}".format(
                np.std(alphas-alphamodel(pixels)),
                np.std(betas-betamodel(pixels)),
                np.std(gammas-gammamodel(pixels)),
                np.std(deltas-deltamodel(pixels)),
                ))

        pl.clf()
        pl.subplot(2,2,1)
        pl.title("Slit {0}".format(solution["slitno"]))
        pl.scatter(pixels, alphas)
        pl.plot(pixels, alphamodel(pixels))
        
        pl.ylim([.993,1/.993])
        pl.xticks(rotation=90)
        pl.ylabel(r'$\alpha$')

        pl.subplot(2,2,2)
        pl.scatter(pixels, betas)
        pl.plot(pixels, betamodel(pixels))
        pl.xticks(rotation=90)
        pl.ylabel(r'$\beta$')

        pl.subplot(2,2,3)
        pl.scatter(pixels, gammas)
        pl.plot(pixels, gammamodel(pixels))
        pl.ylim([0,1e-12])
        pl.xticks(rotation=90)
        pl.ylabel(r'$\gamma$')

        pl.subplot(2,2,4)
        pl.scatter(pixels, deltas)
        pl.plot(pixels, deltamodel(pixels))
        pl.xticks(rotation=90)
        pl.ylabel(r'$\delta$')


        pp.savefig()

    band = header['filter'].rstrip()
    [alpha_pixel, sinbeta_position, sinbeta_pixel, gamma_pixel, 
            delta_pixel] = param_guess_functions(band)
    pl.clf()
    pl.subplot(1,1,1)
    pl.scatter(all_pix, all_alphas, c=all_deltas)
    pl.plot(all_pix, alpha_pixel(all_pix), 'r')

    ff = np.poly1d(np.polyfit(all_pix, all_alphas, 4))
    pl.plot(all_pix, ff(all_pix))
    info("Alpha: "+str(ff))
    pl.ylabel(r'$\alpha$')
    pp.savefig()

    pl.clf()
    delts = all_alphas - ff(all_pix)
    pl.scatter(all_pix, delts, c=all_gammas)
    pl.ylabel(r'$\Delta \alpha$')
    info("Scatter is {0} pixels".format(np.std(delts)*2048))
    pp.savefig()

    pl.clf()
    pl.scatter(all_pix, all_betas, s=.1)
    pl.ylabel(r'$\beta$')
    pp.savefig()

    pl.clf()
    pl.scatter(all_pix, all_gammas, c=all_gammas)
    pl.plot(all_pix, gamma_pixel(all_pix), 'r')
    ff = np.poly1d(np.polyfit(all_pix, all_gammas, 4))
    info("Gamma: "+str(ff))

    pl.plot(all_pix, ff(all_pix), 'b')
    pl.ylabel(r'$\gamma$')
    pp.savefig()

    pl.clf()
    delta_pixel = np.poly1d([4.284e-5, -0.1145, 1219])
    pl.scatter(all_pix, all_deltas, c=all_gammas)
    pl.plot(all_pix, delta_pixel(all_pix), 'r')
    
    ff = np.poly1d(np.polyfit(all_pix, all_deltas, 4))
    info("Delta: "+str(ff))

    pl.ylabel(r'$\delta$')
    pp.savefig()

    pp.close()

if __name__ == "__main__":
    np.set_printoptions(precision=3)


    cc = np.load("lambda_coeffs_m120507_0230.npy")


    px = []
    ys = []
    for c in cc: 
        px.append(c['2d']['positions'])
        ys.append(c['2d']['coeffs'][:,0])







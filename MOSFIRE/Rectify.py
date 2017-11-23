
import os
import sys
import time

import numpy as np
from matplotlib import pyplot as pl
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf

from multiprocessing import Pool
import scipy as sp
import scipy.ndimage
from scipy import interpolate as II
import warnings

import pdb

import MOSFIRE
from MOSFIRE import Background, CSU, Fit, IO, Options, Filters, Detector, Wavelength
from MOSFIRE.MosfireDrpLog import debug, info, warning, error


def handle_rectification(maskname, in_files, wavename, band_pass, files, options,
        commissioning_shift=3.0, target='default'):
    '''Handle slit rectification and coaddition.

    Args:
        maskname: The mask name string
        in_files: List of stacked spectra in electron per second. Will look
            like ['electrons_Offset_1.5.txt.fits', 'electrons_Offset_-1.5.txt.fits']
        wavename: path (relative or full) to the wavelength stack file, string
        band_pass: Band pass name, string
        barset_file: Path to a mosfire fits file containing the full set of
            FITS extensions for the barset. It can be any file in the list
            of science files.
    Returns:
        None

    Writes files:
        [maskname]_[band]_[object name]_eps.fits --
            The rectified, background subtracted, stacked eps spectrum
        [maskname]_[band]_[object name]_sig.fits --
            Rectified, background subtracted, stacked weight spectrum (STD/itime)
        [maskname]_[band]_[object_name]_itime.fits
            Rectified, CRR stacked integration time spectrum
        [maskname]_[band]_[object_name]_snrs.fits
            Rectified signal to noise spectrum
    '''

    global edges, dats, vars, itimes, shifts, lambdas, band, fidl, all_shifts
    band = band_pass

    
    dlambda = Wavelength.grating_results(band)

    hpp = Filters.hpp[band]
    fidl = np.arange(hpp[0], hpp[1], dlambda)

    lambdas = IO.readfits(wavename, options)

    if np.any(lambdas[1].data < 0) or np.any(lambdas[1].data > 29000):
        info("***********WARNING ***********")
        info("The file {0} may not be a wavelength file.".format(wavename))
        info("Check before proceeding.")
        info("***********WARNING ***********")

    edges, meta = IO.load_edges(maskname, band, options)
    shifts = []

    posnames = []
    postoshift = {}
    
    for file in in_files:

        info(":: "+str(file))
        II = IO.read_drpfits(maskname, file, options)

        off = np.array((II[0]["decoff"], II[0]["raoff"]),dtype=np.float64)
        if "yoffset" in II[0]:
            off = -II[0]["yoffset"]
        else:
            # Deal with data taken during commissioning
            if II[0]["frameid"] == 'A': off = 0.0
            else: off = commissioning_shift

        try: off0
        except: off0 = off

        shift = off - off0

        shifts.append(shift)
        posnames.append(II[0]["frameid"])
        postoshift[II[0]['frameid']] = shift
    
        info("Position {0} shift: {1:2.2f} as".format(off, shift))
    # this is to deal with cases in which we want to rectify one single file
    if len(set(posnames)) is 1:
        plans = [['A']]
    else:
        plans = Background.guess_plan_from_positions(set(posnames))

    all_shifts = []
    for plan in plans:
        to_append = []
        for pos in plan:
            to_append.append(postoshift[pos])

        all_shifts.append(to_append)

    # Reverse the elements in all_shifts to deal with an inversion
    all_shifts.reverse()

    theBPM = IO.badpixelmask()

    all_solutions = []
    cntr = 0

    if target is 'default':
        outname = maskname
    else:
        outname = target

    for plan in plans:
        if len(plan) is 1:
            p0 = 'A'
            p1 = 'B'
        else:
            p0 = plan[0].replace("'", "p")
            p1 = plan[1].replace("'", "p")
        suffix = "%s-%s" % (p0,p1)
        info("Handling plan %s" % suffix)
        fname = "bsub_{0}_{1}_{2}.fits".format(outname,band,suffix)
        EPS = IO.read_drpfits(maskname, fname, options)
        EPS[1] = np.ma.masked_array(EPS[1], theBPM, fill_value=0)

        fname = "var_{0}_{1}_{2}.fits".format(outname, band, suffix)
        VAR = IO.read_drpfits(maskname, fname, options)
        VAR[1] = np.ma.masked_array(VAR[1], theBPM, fill_value=np.inf)

        fname = "itime_{0}_{1}_{2}.fits".format(outname, band, suffix)
        ITIME = IO.read_drpfits(maskname, fname, options)
        ITIME[1] = np.ma.masked_array(ITIME[1], theBPM, fill_value=0)


        dats = EPS
        vars = VAR
        itimes = ITIME

        EPS[0]["ORIGFILE"] = fname

        tock = time.time()
        sols = range(len(edges)-1,-1,-1)

        shifts = all_shifts[cntr]
        cntr += 1
        p = Pool()
        solutions = p.map(handle_rectification_helper, sols)
        p.close()

        all_solutions.append(solutions)

    tick = time.time()
    info("-----> Mask took %i. Writing to disk." % (tick-tock))


    output = np.zeros((1, len(fidl)))
    snrs = np.zeros((1, len(fidl)))
    sdout= np.zeros((1, len(fidl)))
    itout= np.zeros((1, len(fidl)))


    # the barset [bs] is used for determining object position
    files = IO.list_file_to_strings(files)
    info("Using "+str(files[0])+" for slit configuration.")
    x, x, bs = IO.readmosfits(files[0], options)
    

    for i_slit in range(len(solutions)):
        solution = all_solutions[0][i_slit]
        header = EPS[0].copy()
        obj = header['OBJECT']
        try:
            target_name = str(bs.ssl[-(i_slit+1)]['Target_Name'], 'utf-8')
        except TypeError:
            target_name = bs.ssl[-(i_slit+1)]['Target_Name']
        header['OBJECT'] = target_name

        pixel_dist = np.float(bs.ssl[-(i_slit+1)]['Target_to_center_of_slit_distance'])/0.18

        pixel_dist -= solution['offset']

        ll = solution["lambda"]

        header["wat0_001"] = "system=world"
        header["wat1_001"] = "wtype=linear"
        header["wat2_001"] = "wtype=linear"
        header["dispaxis"] = 1
        header["dclog1"] = "Transform"
        header["dc-flag"] = 0
        header["ctype1"] = "AWAV"
        header["cunit1"] = "Angstrom"
        header["crval1"] = ll[0]
        header["crval2"] = -solution["eps_img"].shape[0]/2 - pixel_dist
        header["crpix1"] = 1
        header["crpix2"] = 1
        header["cdelt1"] = ll[1]-ll[0]
        header["cdelt2"] = 1
        header["cname1"] = "angstrom"
        header["cname2"] = "pixel"
        header["cd1_1"] = ll[1]-ll[0]
        header["cd1_2"] = 0
        header["cd2_1"] = 0
        header["cd2_2"] = 1


        S = output.shape

        img = solution["eps_img"]
        std = solution["sd_img"]
        tms = solution["itime_img"]


        for i_solution in range(1,len(all_solutions)):
            info("Combining solution %i" %i_solution)
            solution = all_solutions[i_solution][i_slit]
            img += solution["eps_img"]
            std += solution["sd_img"]
            tms += solution["itime_img"]

        output = np.append(output, img, 0)
        output = np.append(output, np.nan*np.zeros((3,S[1])), 0)
        snrs = np.append(snrs, img*tms/std, 0)
        snrs = np.append(snrs, np.nan*np.zeros((3,S[1])), 0)
        sdout = np.append(sdout, std, 0)
        sdout = np.append(sdout, np.nan*np.zeros((3,S[1])), 0)
        itout = np.append(itout, tms, 0)
        itout = np.append(itout, np.nan*np.zeros((3,S[1])), 0)

        header['bunit'] = ('electron/second', 'electron power')
        IO.writefits(img, maskname,
            "{0}_{1}_{2}_eps.fits".format(outname, band, target_name), options,
            overwrite=True, header=header, lossy_compress=False)

        header['bunit'] = ('electron/second', 'sigma/itime')
        IO.writefits(std/tms, maskname,
            "{0}_{1}_{2}_sig.fits".format(outname, band, target_name), options,
            overwrite=True, header=header, lossy_compress=False)

        header['bunit'] = ('second', 'exposure time')
        IO.writefits(tms, maskname,
            "{0}_{1}_{2}_itime.fits".format(outname, band, target_name), options,
            overwrite=True, header=header, lossy_compress=False)

        header['bunit'] = ('', 'SNR')
        IO.writefits(img*tms/std, maskname,
            "{0}_{1}_{2}_snrs.fits".format(outname, band, target_name), options,
            overwrite=True, header=header, lossy_compress=False)

    header = EPS[0].copy()
    header["wat0_001"] = "system=world"
    header["wat1_001"] = "wtype=linear"
    header["wat2_001"] = "wtype=linear"
    header["dispaxis"] = 1
    header["dclog1"] = "Transform"
    header["dc-flag"] = 0
    header["ctype1"] = "AWAV"
    header["cunit1"] = ("Angstrom", 'Start wavelength')
    header["crval1"] = ll[0]
    header["crval2"] = 1
    header["crpix1"] = 1
    header["crpix2"] = 1
    header["cdelt1"] = (ll[1]-ll[0], 'Angstrom/pixel')
    header["cdelt2"] = 1
    header["cname1"] = "angstrom"
    header["cname2"] = "pixel"
    header["cd1_1"] = (ll[1]-ll[0], 'Angstrom/pixel')
    header["cd1_2"] = 0
    header["cd2_1"] = 0
    header["cd2_2"] = 1

    header["bunit"] = "ELECTRONS/SECOND"
    info("############ Final reduced file: {0}_{1}_eps.fits".format(outname,band))
    IO.writefits(output, maskname, "{0}_{1}_eps.fits".format(outname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    header["bunit"] = ""
    IO.writefits(snrs, maskname, "{0}_{1}_snrs.fits".format(outname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    header["bunit"] = "ELECTRONS/SECOND"
    IO.writefits(sdout/itout, maskname, "{0}_{1}_sig.fits".format(outname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    header["bunit"] = "SECOND"
    IO.writefits(itout, maskname, "{0}_{1}_itime.fits".format(outname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)


def r_interpol(ls, ss, lfid, tops, top, shift_pix=0, pad=[0,0], fill_value=0.0):
    '''
    Interpolate the data ss(ls, fs) onto a fiducial wavelength vector.
    ls[n_spatial, n_lam] - wavelength array
    ss[n_spatial, n_lam] - corresponding data array
    lfid[n_lam] - wavelength fiducial to interpolate onto
    shift_pix - # of pixels to shift in spatial direction
    pad - # of pixels to pad in spatial direction
    fill_value - passed through to interp1d
    '''

    S = ss.shape

    output = np.zeros((np.int(S[0]+pad[0]+pad[1]), len(lfid)))
    output[:] = np.nan

    L = np.double(len(lfid))
    
    # First interpolate onto a common wavelength grid
    for i in range(S[0]):

        ll = ls[i,:]
        sp = ss[i,:]
        ok = np.where(ll>1000)[0]

        if len(ok) >= 100:
            f = II.interp1d(ll[ok], sp[ok], bounds_error=False, 
                fill_value = fill_value)
            output[i,:] = f(lfid)


    # Now rectify in spatial
    vert_shift = tops-top-shift_pix

    f = II.interp1d(ls[10, :], vert_shift, bounds_error=False, 
        fill_value = fill_value)

    for i in range(output.shape[1]):
        to_shift = f(fidl[i])
        x = np.arange(output.shape[0])
        y = II.interp1d(x, output[:, i], bounds_error=False,
            fill_value=fill_value)

        output[:,i] = y(x + to_shift)

    return output


def handle_rectification_helper(edgeno):
    ''' All the rectification happens in this helper function. This helper function
    is spawned as a separate process in the multiprocessing pool'''

    global edges, dats, vars, itimes, shifts, lambdas, band, fidl,all_shifts

    pix = np.arange(2048)
    
    edge = edges[edgeno]

    info("Handling edge: "+str(edge["Target_Name"]))

    tops = edge["top"](pix)
    bots = edge["bottom"](pix)

    # Length of the slit in arcsecond
    lenas = (tops[1024] - bots[1024]) * 0.18
    mxshift = np.abs(np.int(np.ceil(np.max(all_shifts)/0.18)))
    mnshift = np.abs(np.int(np.floor(np.min(all_shifts)/0.18)))
    
    top = int(min(np.floor(np.min(tops)), 2048))
    bot = int(max(np.ceil(np.max(bots)), 0))

    ll = lambdas[1].data[bot:top, :]
    eps = dats[1][bot:top, :].filled(0.0)
    vv = vars[1][bot:top, :].filled(np.inf)
    it  = itimes[1][bot:top, :].filled(0.0)

    lmid = ll[ll.shape[0]//2,:]
    hpp = Filters.hpp[band]

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]

    epss = []
    ivss = []
    itss = []
    if len(shifts) is 1: sign = 1
    else: sign = -1
    for shift in shifts:
        output = r_interpol(ll, eps, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value = np.nan)
        epss.append(sign * output)

        ivar = 1/vv
        bad = np.where(np.isfinite(ivar) ==0)
        ivar[bad] = 0.0
        output = r_interpol(ll, ivar, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value=np.nan) 
        ivss.append(output)

        output = r_interpol(ll, it, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value=np.nan) 
        itss.append(output)

        sign *= -1
    # the "mean of empty slice" warning are generated at the top and bottom edges of the array
    # where there is basically no data due to the shifts between a and b positions
    # we could pad a little bit less, or accept the fact that the slits have a couple of rows of
    # nans in the results.
    warnings.filterwarnings('ignore','Mean of empty slice')
    it_img = np.nansum(np.array(itss), axis=0)
    eps_img = np.nanmean(epss, axis=0)
    warnings.filterwarnings('always')

    # Remove any NaNs or infs from the variance array
    ivar_img = []
    for ivar in ivss:
        bad = np.where(np.isfinite(ivar) == 0)
        ivar[bad] = 0.0

        ivar_img.append(ivar)

    IV = np.array(ivar_img)
    bad = np.isclose(IV,0)
    IV[bad] = np.inf
    var_img = np.nanmean(1/np.array(IV), axis=0)
    sd_img = np.sqrt(var_img)

    return {"eps_img": eps_img, "sd_img": sd_img, "itime_img": it_img, 
            "lambda": fidl, "Target_Name": edge["Target_Name"], 
            "slitno": edgeno+1, "offset": np.max(tops-top)}


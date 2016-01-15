
import os
import sys
import time
import warnings

import numpy as np
import pylab as pl
import pyfits as pf
from multiprocessing import Pool
import scipy as sp
import scipy.ndimage
from scipy import interpolate as II

import pdb

import MOSFIRE
from MOSFIRE import CSU, Fit, IO, Options, Filters, Detector, Wavelength
from MosfireDrpLog import debug, info, warning, error

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed()

def rem_header_key(header, key):

    try:
        del header[key]
    except:
        return False

    return True

def guess_plan_from_positions(posnames):
    ''' Based on a Set() of position names guess the observing plan.
        e.g., position names Set(["A", "B", "A'", "B'"]) -> "A-B", "A'-B'" '''
    if posnames == set(["A", "B"]): return [["A", "B"]]
    elif posnames == set(["A'", "B'", "A", "B"]): 
        return [["A", "B"], ["A'", "B'"]]
    else:
        raise Exception("Could not get observing plan from positions %s. You must use the plan keyword" % posnames)


def imcombine(files, maskname, options, flat, outname=None, shifts=None,
    extension=None):
    '''
    From a list of files it imcombine returns the imcombine of several values.
    The imcombine code also estimates the readnoise ad RN/sqrt(numreads) so
    that the variance per frame is equal to (ADU + RN^2) where RN is computed
    in ADUs.

    Arguments:
        files[]: list of full path to files to combine
        maskname: Name of mask
        options: Options dictionary
        flat[2048x2048]: Flat field (values should all be ~ 1.0)
        outname: If set, will write (see notes below for details)
            eps_[outname].fits: electron/sec file
            itimes_[outname].fits: integration time
            var_[outname].fits: Variance files
        shifts[len(files)]: If set, will "roll" each file by the 
            amount in the shifts vector in pixels. This argument
            is used when telescope tracking is poor. If you need
            to use this, please notify Keck staff about poor 
            telescope tracking.

    Returns 6-element tuple:
        header: The combined header
        electrons [2048x2048]:  e- (in e- units)
        var [2048x2048]: electrons + RN**2 (in e-^2 units)
        bs: The MOSFIRE.Barset instance
        itimes [2048x2048]: itimes (in s units)
        Nframe: The number of frames that contribute to the summed
            arrays above. If Nframe > 5 I use the sigma-clipping
            Cosmic Ray Rejection tool. If Nframe < 5 then I drop
            the max/min elements.

    Notes:

        header -- fits header
        ADUs -- The mean # of ADUs per frame
        var -- the Variance [in adu] per frame. 
        bs -- Barset
        itimes -- The _total_ integration time in second
        Nframe -- The number of frames in a stack.

        
        Thus the number of electron per second is derived as: 
            e-/sec = (ADUs * Gain / Flat) * (Nframe/itimes)

        The total number of electrons is:
            el = ADUs * Gain * Nframe


    '''

    ADUs = np.zeros((len(files), 2048, 2048))
    itimes = np.zeros((len(files), 2048, 2048))
    prevssl = None
    prevmn = None
    patternid = None
    maskname = None

    header = None

    if shifts is None:
        shifts = np.zeros(len(files))

    warnings.filterwarnings('ignore')
    for i in xrange(len(files)):
        fname = files[i]
        thishdr, data, bs = IO.readmosfits(fname, options, extension=extension)
        itimes[i,:,:] = thishdr["truitime"]

        base = os.path.basename(fname).rstrip(".fits")
        fnum = int(base.split("_")[1])
        
        if shifts[i] == 0:
            ADUs[i,:,:] = data.filled(0.0) / flat
        else:
            ADUs[i,:,:] = np.roll(data.filled(0.0) / flat, np.int(shifts[i]), axis=0)

        ''' Construct Header'''
        if header is None:
            header = thishdr

        header["imfno%3.3i" % (fnum)] =  (fname, "img%3.3i file name" % fnum)

        map(lambda x: rem_header_key(header, x), ["CTYPE1", "CTYPE2", "WCSDIM",
            "CD1_1", "CD1_2", "CD2_1", "CD2_2", "LTM1_1", "LTM2_2", "WAT0_001",
            "WAT1_001", "WAT2_001", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
            "RADECSYS"])

        for card in header.cards:
            if card == '': continue
            key,val,comment = card
            
            if key in thishdr:
                if val != thishdr[key]:
                    newkey = key + ("_img%2.2i" % fnum)
                    try: header[newkey.rstrip()] = (thishdr[key], comment)
                    except: pass

        ''' Now handle error checking'''

        if maskname is not None:
            if thishdr["maskname"] != maskname:
                error("File %s uses mask '%s' but the stack is of '%s'" %
                    (fname, thishdr["maskname"], maskname))
                raise Exception("File %s uses mask '%s' but the stack is of '%s'" %
                    (fname, thishdr["maskname"], maskname))

        maskname = thishdr["maskname"]
            
        if thishdr["aborted"]:
            error("Img '%s' was aborted and should not be used" %
                    fname)
            raise Exception("Img '%s' was aborted and should not be used" %
                    fname)

        if prevssl is not None:
            if len(prevssl) != len(bs.ssl):
                # todo Improve these checks
                error("The stack of input files seems to be of "
                        "different masks")
                raise Exception("The stack of input files seems to be of "
                        "different masks")
        prevssl = bs.ssl

        if patternid is not None:
            if patternid != thishdr["frameid"]:
                error("The stack should be of '%s' frames only, but "
                        "the current image is a '%s' frame." % (patternid, 
                            thishdr["frameid"]))
                raise Exception("The stack should be of '%s' frames only, but "
                        "the current image is a '%s' frame." % (patternid, 
                            thishdr["frameid"]))

        patternid = thishdr["frameid"]


        if maskname is not None:
            if maskname != thishdr["maskname"]:
                error("The stack should be of CSU mask '%s' frames "
                        "only but contains a frame of '%s'." % (maskname,
                        thishdr["maskname"]))
                raise Exception("The stack should be of CSU mask '%s' frames "
                        "only but contains a frame of '%s'." % (maskname,
                        thishdr["maskname"]))

        maskname = thishdr["maskname"]

        if thishdr["BUNIT"] != "ADU per coadd":
            error("The units of '%s' are not in ADU per coadd and "
                    "this violates an assumption of the DRP. Some new code " 
                    "is needed in the DRP to handle the new units of "
                    "'%s'." % (fname, thishdr["BUNIT"]))
            raise Exception("The units of '%s' are not in ADU per coadd and "
                    "this violates an assumption of the DRP. Some new code " 
                    "is needed in the DRP to handle the new units of "
                    "'%s'." % (fname, thishdr["BUNIT"]))

        ''' Error checking is complete'''
        info("%s %s[%s]/%s: %5.1f s,  Shift: %i px" % (fname, maskname, patternid,
            header['filter'], np.mean(itimes[i]), shifts[i]))

    warnings.filterwarnings('always')

    # the electrons and el_per_sec arrays are:
    #   [2048, 2048, len(files)] and contain values for
    # each individual frame that is being combined.
    # These need to be kept here for CRR reasons.
    electrons = np.array(ADUs) * Detector.gain 
    el_per_sec = electrons / itimes

    output = np.zeros((2048, 2048))
    exptime = np.zeros((2048, 2048))

    numreads = header["READS0"]
    RN_adu = Detector.RN / np.sqrt(numreads) / Detector.gain
    RN = Detector.RN / np.sqrt(numreads)

    # Cosmic ray rejection code begins here. This code construction the
    # electrons and itimes arrays.
    standard = True
    new_from_chuck = False
    # Chuck Steidel has provided a modified version of the CRR procedure. 
    # to enable it, modify the variables above.
    
    if new_from_chuck and not standard:
        if len(files) >= 5:
            print "Sigclip CRR"
            srt = np.argsort(electrons, axis=0, kind='quicksort')
            shp = el_per_sec.shape
            sti = np.ogrid[0:shp[0], 0:shp[1], 0:shp[2]]

            electrons = electrons[srt, sti[1], sti[2]]
            el_per_sec = el_per_sec[srt, sti[1], sti[2]]
            itimes = itimes[srt, sti[1], sti[2]]

            # Construct the mean and standard deviation by dropping the top and bottom two 
            # electron fluxes. This is temporary.
            mean = np.mean(el_per_sec[1:-1,:,:], axis = 0)
            std = np.std(el_per_sec[1:-1,:,:], axis = 0)

            drop = np.where( (el_per_sec > (mean+std*4)) | (el_per_sec < (mean-std*4)) )
            print "dropping: ", len(drop[0])
            electrons[drop] = 0.0
            itimes[drop] = 0.0

            electrons = np.sum(electrons, axis=0)
            itimes = np.sum(itimes, axis=0)
            Nframe = len(files) 

        elif len(files) > 5:
            print "WARNING: Drop min/max CRR"
            srt = np.argsort(el_per_sec,axis=0)
            shp = el_per_sec.shape
            sti = np.ogrid[0:shp[0], 0:shp[1], 0:shp[2]]

            electrons = electrons[srt, sti[1], sti[2]]
            itimes = itimes[srt, sti[1], sti[2]]

            electrons = np.sum(electrons[1:-1,:,:], axis=0)
            itimes = np.sum(itimes[1:-1,:,:], axis=0)

            Nframe = len(files) - 2

        else:
            warning( "With less than 5 frames, the pipeline does NOT perform")
            warning( "Cosmic Ray Rejection.")
            # the "if false" line disables cosmic ray rejection"
            if False: 
                for i in xrange(len(files)):
                    el = electrons[i,:,:]
                    it = itimes[i,:,:]
                    el_mf = scipy.signal.medfilt(el, 5)

                    bad = np.abs(el - el_mf) / np.abs(el) > 10.0
                    el[bad] = 0.0
                    it[bad] = 0.0

                    electrons[i,:,:] = el
                    itimes[i,:,:] = it

            electrons = np.sum(electrons, axis=0)
            itimes = np.sum(itimes, axis=0)
            Nframe = len(files) 

    if standard and not new_from_chuck:
        if len(files) >= 9:
            info("Sigclip CRR")
            srt = np.argsort(electrons, axis=0, kind='quicksort')
            shp = el_per_sec.shape
            sti = np.ogrid[0:shp[0], 0:shp[1], 0:shp[2]]

            electrons = electrons[srt, sti[1], sti[2]]
            el_per_sec = el_per_sec[srt, sti[1], sti[2]]
            itimes = itimes[srt, sti[1], sti[2]]

            # Construct the mean and standard deviation by dropping the top and bottom two 
            # electron fluxes. This is temporary.
            mean = np.mean(el_per_sec[2:-2,:,:], axis = 0)
            std = np.std(el_per_sec[2:-2,:,:], axis = 0)

            drop = np.where( (el_per_sec > (mean+std*4)) | (el_per_sec < (mean-std*4)) )
            info("dropping: "+str(len(drop[0])))
            electrons[drop] = 0.0
            itimes[drop] = 0.0

            electrons = np.sum(electrons, axis=0)
            itimes = np.sum(itimes, axis=0)
            Nframe = len(files) 

        elif len(files) > 5:
            warning( "WARNING: Drop min/max CRR")
            srt = np.argsort(el_per_sec,axis=0)
            shp = el_per_sec.shape
            sti = np.ogrid[0:shp[0], 0:shp[1], 0:shp[2]]

            electrons = electrons[srt, sti[1], sti[2]]
            itimes = itimes[srt, sti[1], sti[2]]

            electrons = np.sum(electrons[1:-1,:,:], axis=0)
            itimes = np.sum(itimes[1:-1,:,:], axis=0)

            Nframe = len(files) - 2

        else:
            warning( "With less than 5 frames, the pipeline does NOT perform")
            warning( "Cosmic Ray Rejection.")
            # the "if false" line disables cosmic ray rejection"
            if False: 
                for i in xrange(len(files)):
                     el = electrons[i,:,:]
                     it = itimes[i,:,:]
                     # calculate the median image
                     el_mf = scipy.signal.medfilt(el, 5)
                     el_mf_large = scipy.signal.medfilt(el_mf, 15)
                     # LR: this is a modified version I was experimenting with. For the version 
                     #     written by Nick, see the new_from_chuck part of this code
                     # sky sub
                     el_sky_sub = el_mf - el_mf_large
                     # add a constant value
                     el_plus_constant = el_sky_sub + 100

                     bad = np.abs(el - el_mf) / np.abs(el_plus_constant) > 50.0
                     el[bad] = 0.0
                     it[bad] = 0.0

                     electrons[i,:,:] = el
                     itimes[i,:,:] = it

            electrons = np.sum(electrons, axis=0)
            itimes = np.sum(itimes, axis=0)
            Nframe = len(files) 


    ''' Now handle variance '''
    numreads = header["READS0"]
    RN_adu = Detector.RN / np.sqrt(numreads) / Detector.gain
    RN = Detector.RN / np.sqrt(numreads)

    var = (electrons + RN**2) 

    ''' Now mask out bad pixels '''
    electrons[data.mask] = np.nan
    var[data.mask] = np.inf

    if "RN" in header:
        error("RN Already populated in header")
        raise Exception("RN Already populated in header")
    header['RN'] = ("%1.3f" , "Read noise in e-")
    header['NUMFRM'] = (Nframe, 'Typical number of frames in stack')


    header['BUNIT'] = 'ELECTRONS/SECOND'
    IO.writefits(np.float32(electrons/itimes), maskname, "eps_%s" % (outname),
                 options, header=header, overwrite=True)

    # Update itimes after division in order to not introduce nans
    itimes[data.mask] = 0.0

    header['BUNIT'] = 'ELECTRONS^2'
    IO.writefits(var, maskname, "var_%s" % (outname),
                 options, header=header, overwrite=True, lossy_compress=True)

    header['BUNIT'] = 'SECOND'
    IO.writefits(np.float32(itimes), maskname, "itimes_%s" % (outname),
                options, header=header, overwrite=True, lossy_compress=True)

    return header, electrons, var, bs, itimes, Nframe

def merge_headers(h1, h2):
    """Merge headers h1 and h2 such that h2 has the nod position name
    appended"""
    
    h = h1.copy()

    patternid = h2["frameid"]

    for key, val, comment in h2.cards:
        if "NAXIS" in key: continue
        if "SIMPLE" in key: continue
        if "BITPIX" in key: continue
        if "EXTEND" in key: continue

        if key in h:
            continue
        else:
            try: h[key] = (val, comment)
            except: pass


    return h


def handle_background(filelist, wavename, maskname, band_name, options, shifts=None, plan=None, extension=None, target='default'): 
    '''
    Perform difference imaging and subtract residual background.

    The plan looks something like: [['A', 'B']]
    In this case, the number of output files is equal to the length of the list (1).

    If you choose to use an ABA'B' pattern then the plan will be: [["A", "B"], ["A'", "B'"]]
    the background subtraction code will make and handle two files, "A-B" and "A'-B'".
    '''
    
    global header, bs, edges, data, Var, itime, lam, sky_sub_out, sky_model_out, band

    band = band_name

    flatname = "pixelflat_2d_%s.fits" % band_name
    hdr, flat = IO.readfits("pixelflat_2d_%s.fits" % (band_name), options)

    if np.abs(np.median(flat) - 1) > 0.1:
        error("Flat seems poorly behaved.")
        raise Exception("Flat seems poorly behaved.")

    '''
        This next section of the code figures out the observing plan
        and then deals with the bookeeping of sending the plan
        to the background subtracter.
    '''

    hdrs = []
    epss = {}
    vars = {}
    bss = []
    times = {}
    Nframes = []

    i = 0
    header = pf.Header()
    for i in xrange(len(filelist)):
        fl = filelist[i]
        files = IO.list_file_to_strings(fl)
        info("Combining")
        if shifts is None: shift = None
        else: shift = shifts[i]
        hdr, electron, var, bs, time, Nframe = imcombine(files, maskname,
            options, flat, outname="%s.fits" % (fl),
            shifts=shift, extension=extension)

        hdrs.append(hdr) 
        header = merge_headers(header, hdr)
        epss[hdr['FRAMEID']] = electron/time
        vars[hdr['FRAMEID']] = var
        times[hdr['FRAMEID']] = time
        bss.append(bs)
        Nframes.append(Nframe)

    positions = {}
    i = 0
    for h in hdrs:
        positions[h['FRAMEID']] = i
        i += 1
    posnames = set(positions.keys())
    if plan is None:
        plan = guess_plan_from_positions(posnames)

    num_outputs = len(plan)


    edges, meta = IO.load_edges(maskname, band, options)
    lam = IO.readfits(wavename, options)

    bs = bss[0]

    for i in xrange(num_outputs):
        posname0 = plan[i][0]
        posname1 = plan[i][1]
        info("Handling %s - %s" % (posname0, posname1))
        data = epss[posname0] - epss[posname1]
        Var = vars[posname0] + vars[posname1]
        itime = np.mean([times[posname0], times[posname1]], axis=0)

        p = Pool()
        solutions = p.map(background_subtract_helper, xrange(len(bs.ssl)))
        p.close()

        write_outputs(solutions, itime, header, maskname, band, plan[i], options, target=target)


def write_outputs(solutions, itime, header, maskname, band_name, plan, options, target):
    sky_sub_out = np.zeros((2048, 2048), dtype=np.float)
    sky_model_out = np.zeros((2048, 2048), dtype=np.float)

    p0 = plan[0].replace("'", "p")
    p1 = plan[1].replace("'", "p")
    suffix = "%s-%s" % (p0,p1)
    xroi = slice(0,2048)

    for sol in solutions:
        if not sol["ok"]: 
            continue

        yroi = slice(sol["bottom"], sol["top"])
        sky_sub_out[yroi, xroi] = sol["output"]
        sky_model_out[yroi, xroi] = sol["model"]

    if target is 'default':
        outname = maskname
    else:
        outname = target
        
    header['BUNIT'] = 'SECOND'
    IO.writefits(itime, maskname, "itime_%s_%s_%s.fits" % (outname, band,
        suffix), options, header=header, overwrite=True, lossy_compress=True)


    header['BUNIT'] = 'ELECTRONS/SECOND'
    IO.writefits(data, maskname, "sub_%s_%s_%s.fits" % (outname, band,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    header['BUNIT'] = 'ELECTRONS/SECOND'
    IO.writefits(sky_sub_out, maskname, "bsub_%s_%s_%s.fits" % (outname, band,
        suffix), options, header=header, overwrite=True)

    header['BUNIT'] = 'ELECTRONS'
    IO.writefits(Var, maskname, "var_%s_%s_%s.fits" % (outname, band,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    header['BUNIT'] = 'ELECTRONS/SECOND'
    IO.writefits(sky_model_out, maskname, "bmod_%s_%s_%s.fits" % (outname,
        band, suffix), options, header=header, overwrite=True,
        lossy_compress=True)

    '''Now create rectified solutions'''
    dlam = Wavelength.grating_results(band)
    hpp = np.array(Filters.hpp[band]) 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)
    nspec = len(ll_fid)


    rectified = np.zeros((2048, nspec), dtype=np.float32)
    rectified_var = np.zeros((2048, nspec), dtype=np.float32)
    rectified_itime = np.zeros((2048, nspec), dtype=np.float32)

    from scipy.interpolate import interp1d
    for i in xrange(2048):
        ll = lam[1][i,:]
        ss = sky_sub_out[i,:]

        ok = np.isfinite(ll) & np.isfinite(ss) & (ll < hpp[1]) & (ll >
                hpp[0])

        if len(np.where(ok)[0]) < 100:
            continue
        f = interp1d(ll[ok], ss[ok], bounds_error=False)
        rectified[i,:] = f(ll_fid)

        f = interp1d(ll, Var[i,:], bounds_error=False)
        rectified_var[i,:] = f(ll_fid)

        f = interp1d(ll, itime[i,:], bounds_error=False)
        rectified_itime[i,:] = f(ll_fid)

    header["wat0_001"] = "system=world"
    header["wat1_001"] = "wtype=linear"
    header["wat2_001"] = "wtype=linear"
    header["dispaxis"] = 1
    header["dclog1"] = "Transform"
    header["dc-flag"] = 0
    header["ctype1"] = "AWAV"
    header["cunit1"] = "Angstrom"
    header["crval1"] = (ll_fid[0], "Starting wavelength Angstrom")
    header["crval2"] = 0
    header["crpix1"] = 1
    header["crpix2"] = 1
    header["cdelt1"] = 1
    header["cdelt2"] = 1
    header["cname1"] = "angstrom"
    header["cname2"] = "pixel"
    header["cd1_1"] = (dlam, "Angstrom/pixel")
    header["cd1_2"] = 0
    header["cd2_1"] = 0
    header["cd2_2"] = (1, "pixel/pixel")

    IO.writefits(rectified_itime, maskname,
        "%s_rectified_itime_%s_%s.fits" % (outname, band_name,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    IO.writefits(rectified, maskname, "%s_rectified_%s_%s.fits" % (outname,
        band_name, suffix), options, header=header, overwrite=True,
        lossy_compress=True)

    IO.writefits(rectified_var, maskname, "%s_rectified_var_%s_%s.fits" %
        (outname, band_name, suffix), options, header=header, overwrite=True,
        lossy_compress=True)

    IO.writefits(rectified*rectified_itime/np.sqrt(rectified_var), maskname,
        "%s_rectified_sn_%s_%s.fits" % (outname, band_name,
        suffix), options, header=header, overwrite=True, lossy_compress=True)


def background_subtract_helper(slitno):
    '''

    Background subtraction follows the methods outlined by Kelson (2003). Here
    a background is estimated as a function of wavelength using B-splines and
    subtracted off. The assumption is that background is primarily a function
    of wavelength, and thus by sampling the background across the full 2-d
    spectrum the background is sampled at much higher than the native spectral
    resolution of mosfire. 

    Formally, the assumption that background is only a function of wavelength
    is incorrect, and indeed a "transmission function" is estimated from the 2d
    spectrum. This gives an estimate of the throughput of the slit and divided
    out.

    1. Extract the slit from the 2d image.
    2. Convert the 2d spectrum into a 1d spectrum
    3. Estimate transmission function

    '''

    global header, bs, edges, data, Var, itime, lam, sky_sub_out, sky_model_out, band
    tick = time.time()

    # 1
    top = np.int(edges[slitno]["top"](1024))  
    bottom = np.int(edges[slitno]["bottom"](1024)) 
    info("Background subtracting slit %i [%i,%i]" % (slitno, top, bottom))

    pix = np.arange(2048)
    xroi = slice(0,2048)
    yroi = slice(bottom, top)

    stime = itime[yroi, xroi]
    slit  = data[yroi, xroi] 
    Var[np.logical_not(np.isfinite(Var))] = np.inf

    lslit = lam[1][yroi,xroi]

    # 2
    xx = np.arange(slit.shape[1])
    yy = np.arange(slit.shape[0])

    X,Y = np.meshgrid(xx,yy)

    train_roi = slice(5,-5)
    ls = lslit[train_roi].flatten().filled(0)
    ss = slit[train_roi].flatten()
    ys = Y[train_roi].flatten()

    dl = np.ma.median(np.diff(lslit[lslit.shape[0]/2,:]))[0]
    if dl == 0:
        return {"ok": False}

    sort = np.argsort(ls)
    ls = ls[sort]
    ys = ys[sort]

    hpps = np.array(Filters.hpp[band]) 

    diff = np.append(np.diff(ls), False)

    OK = (diff > 0.001) & (ls > hpps[0]) & (ls < hpps[1]) & (np.isfinite(ls)) \
            & (np.isfinite(ss[sort]))

    if len(np.where(OK)[0]) < 1000:
        warning("Failed on slit "+str(slitno))
        return {"ok": False}

    # 3
    pp = np.poly1d([1.0])
    ss = (slit[train_roi] / pp(Y[train_roi])).flatten()
    ss = ss[sort]

    knotstart = max(hpps[0], min(ls[OK])) + 5
    knotend = min(hpps[1], max(ls[OK])) - 5


    for i in range(3):
        try:
            delta = dl*0.9
            knots = np.arange(knotstart, knotend, delta)
            bspline = II.splrep(ls[OK], ss[OK], k=5, task=-1, t=knots)
        except:
            delta = dl*1.4
            knots = np.arange(knotstart, knotend, delta)
            try:
                bspline = II.splrep(ls[OK], ss[OK], k=5, task=-1, t=knots)
            except:
                warning("Could not construct spline on slit "+str(slitno))
                return {"ok": False}

        ll = lslit.flatten()
        model = II.splev(ll, bspline)

        oob = np.where((ll < knotstart) | (ll > knotend))
        model[oob] = np.median(ss)
        model = model.reshape(slit.shape)

        output = slit - model

        std = np.abs(output)/(np.sqrt(np.abs(model)+1))

        tOK = (std[train_roi] < 10).flatten() & \
                np.isfinite(std[train_roi]).flatten()  
        OK = OK & tOK[sort]

    return {"ok": True, "slitno": slitno, "bottom": bottom, "top": top,
            "output": output, "model": model, "bspline": bspline}


if __name__ == "__main__":
    background_subtract()



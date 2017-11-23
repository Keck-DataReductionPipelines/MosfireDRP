'''

===================
MOSFIRE Flat Fields 
===================




npk April 14th 2011
Modifications:

2013-08-29:  T. Do - added 'edgeThreshold' keyword as well as changed the initial guess for the edge location when using long slits

2014-01-09:  MK - Functions added to subtract dome "lamps off/thermal" flats from the dome "lamps on" flats. 
                  The functions combine both sets of flats using your current method to create a lamps on and 
                  lamps off flat, and then subtracts those two images to remove the contribution of the dome 
                  emission. Files renamed combflat_lamps_on* and combflat_lamps_off*. The final flat has the 
                  same name that you output: combflat_2s_band_.fits .

                  To reduce the thermal flats, the flat functions have an optional keyword for the 
                  "lampOffList" that acts as the reduction trigger. The driver file should include a call 
                  like the example below.
                      Flats.handle_flats('Flat.txt', maskname, band, flatops, lampOffList='FlatThermal.txt')

2014-04-27:  MK - Nick released a new version of the DRP in March. I merged the differences between the code 
                  that I downloaded March 14 (Nick K.'s updated code to use file lists) with the code we 
                  developed to use the thermal flatts.    

2014-06-12:  MK - Nick released a new version of the DRP in June 2013. I merged the differences between the code 
                  modified on April 27 and the 10 June release. Tested with K-band data and appears to work
                  as expected. 
'''

import os
import time
import unittest

import numpy as np
from matplotlib import pyplot as pl
import scipy, scipy.ndimage
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf

import pdb

import MOSFIRE
from MOSFIRE import Fit, IO, Options, CSU, Wavelength, Filters, Detector
from MOSFIRE.MosfireDrpLog import debug, info, warning, error
__version__ = 0.1

#from IPython.Shell import IPShellEmbed
#start_shell = IPShellEmbed()

def handle_flats(flatlist, maskname, band, options, extension=None,edgeThreshold=450,lampOffList=None,longslit=None, verbose=False):
    '''
    handle_flats is the primary entry point to the Flats module.

    handle_flats takes a list of individual exposure FITS files and creates:
    1. A CRR, dark subtracted, pixel-response flat file.
    2. A set of polynomials that mark the edges of a slit

    Inputs:
    flatlist: 
    maskname: The name of a mask
    band: A string indicating the bandceil

    Outputs:

    file {maskname}/flat_2d_{band}.fits -- pixel response flat
    file {maskname}/edges.np
    '''

    tick = time.time()

    # Check
    bpos = np.ones(92) * -1

    #Retrieve the list of files to use for flat creation.
    flatlist = IO.list_file_to_strings(flatlist)
    if len(flatlist) == 0:
        print('WARNING: No flat files found.')
        raise IOError('No flat files found')
    # Print the filenames to Standard-out
    for flat in flatlist:
        if verbose==True:
            print(str(flat))
        debug(str(flat))

    #Determine if flat files headers are in agreement
    for fname in flatlist:

        hdr, dat, bs = IO.readmosfits(fname, options, extension=extension)
        try: bs0
        except: bs0 = bs

        if np.any(bs0.pos != bs.pos):
            print("bs0: "+str(bs0.pos)+" bs: "+str(bs.pos))
            error("Barset do not seem to match")
            raise Exception("Barsets do not seem to match")

        if hdr["filter"] != band:
            error ("Filter name %s does not match header filter name "
                    "%s in file %s" % (band, hdr["filter"], fname))
            raise Exception("Filter name %s does not match header filter name "
                    "%s in file %s" % (band, hdr["filter"], fname))
        for i in range(len(bpos)):
            b = hdr["B{0:02d}POS".format(i+1)]
            if bpos[i] == -1:
                bpos[i] = b
            else:
                if bpos[i] != b:
                    error("Bar positions are not all the same in "
                            "this set of flat files")
                    raise Exception("Bar positions are not all the same in "
                            "this set of flat files")
    bs = bs0

    # Imcombine the lamps ON flats
    info("Attempting to combine files in {}".format(flatlist))
    if verbose==True: print("Attempting to combine files in {}".format(flatlist))
    out = os.path.join("combflat_2d_{:s}.fits".format(band))
    IO.imcombine(flatlist, out, options, reject="minmax", nlow=1, nhigh=1)

    # Imcombine the lamps OFF flats and subtract the off from the On sets
    if lampOffList != None: 
        #Retrieve the list of files to use for flat creation. 
        info("Attempting to combine Lamps off files in {}".format(lampOffList))
        if verbose==True: print("Attempting to combine Lamps off files in {}".format(lampOffList))
        lampOffList = IO.list_file_to_strings(lampOffList)
        for flat in lampOffList:
            debug(str(flat))
        out = os.path.join("combflat_lamps_off_2d_{:s}.fits".format(band))
        IO.imcombine(flatlist, out, options, reject="minmax", nlow=1, nhigh=1)
        file_on = os.path.join("combflat_2d_{:s}.fits".format(band))
        file_off = os.path.join("combflat_lamps_off_2d_{:s}.fits".format(band))
        file_on_save = os.path.join("combflat_lamps_on_2d_{:s}.fits".format(band))
        IO.imarith(file_on, '-', file_off, file_on_save)

    debug("Combined '%s' to '%s'" % (flatlist, maskname))
    if verbose==True: print("Combined '%s' to '%s'" % (flatlist, maskname))
#     info("Combined flats for '%s'" % (maskname))
    path = "combflat_2d_%s.fits" % band
    (header, data) = IO.readfits(path, use_bpm=True)
    info("Flat written to %s" % path)
    if verbose==True: print("Flat written to %s" % path)

    # Edge Trace
    if bs.long_slit:
        info( "Long slit mode recognized")
        if verbose==True: print( "Long slit mode recognized")
        info( "Central row position:   "+str(longslit["row_position"]))
        info( "Upper and lower limits: "+str(longslit["yrange"][0])+" "+str(longslit["yrange"][1]))
        results = find_longslit_edges(data,header, bs, options, edgeThreshold=edgeThreshold, longslit=longslit)
    elif bs.long2pos_slit:
        info( "Long2pos mode recognized")
        if verbose==True: print( "Long2pos mode recognized")
        results = find_long2pos_edges(data,header, bs, options, edgeThreshold=edgeThreshold, longslit=longslit)
    else:
        info('Finding slit edges in {}'.format(path))
        if verbose==True: print('Finding slit edges in {}'.format(path))
        results = find_and_fit_edges(data, header, bs, options,edgeThreshold=edgeThreshold)
    results[-1]["maskname"] = maskname
    results[-1]["band"] = band
    np.save("slit-edges_{0}".format(band), results)
    save_ds9_edges(results, options)

    # Generate Flat
    out = "pixelflat_2d_%s.fits" % (band)
    if lampOffList != None: 
         make_pixel_flat(data, results, options, out, flatlist, lampsOff=True)
    else:
         make_pixel_flat(data, results, options, out, flatlist, lampsOff=False)

    info( "Pixel flat took {0:6.4} s".format(time.time()-tick))
    if verbose==True: print( "Pixel flat took {0:6.4} s".format(time.time()-tick))

    

def make_pixel_flat(data, results, options, outfile, inputs, lampsOff=None):
    '''
    Convert a flat image into a flat field
    '''

    def pixel_min(y): return int(np.floor(np.min(y)))
    def pixel_max(y): return int(np.ceil(np.max(y)))

    def collapse_flat_box(dat):
        '''Collapse data to the spectral axis (0)'''
        v = np.median(dat, axis=0).ravel()

        return v

    flat = np.ones(shape=Detector.npix)

    hdu = pf.PrimaryHDU((data/flat).astype(np.float32))
    hdu.header.set("version", __version__, "DRP version")
    i = 0
    for flatname in inputs:
        nm = flatname.split("/")[-1]
        hdu.header.set("infile%2.2i" % i, nm)
        i += 1

    slitno = 0
    for result in results[0:-1]:
        slitno += 1
        try:
            hdu.header.set("targ%2.2i" % slitno, result["Target_Name"])
        except ValueError: hdu.header.set("targ%2.2i" % slitno, str(result["Target_Name"],  'utf-8'))

        bf = result["bottom"]
        tf = result["top"]
        try:
            hpps = result["hpps"]
        except:
            error( "No half power points for this slit")
            hpps = [0, Detector.npix[0]]

        xs = np.arange(hpps[0], hpps[1])

        top = pixel_min(tf(xs))
        bottom = pixel_max(bf(xs))

        hdu.header.set("top%2.2i" % slitno, top)
        hdu.header.set("bottom%2.2i" % slitno, bottom)

        info( "%s] Bounding top/bottom: %i/%i" % (result["Target_Name"],
                bottom, top))

        v = collapse_flat_box(data[bottom:top,hpps[0]:hpps[1]])

        x2048 = np.arange(Options.npix)
        v = np.poly1d(np.polyfit(xs,v,
            options['flat-field-order']))(xs).ravel()

        for i in np.arange(bottom-1, top+1):
            flat[i,hpps[0]:hpps[1]] = v

    info("Producing Pixel Flat...")
    for r in range(len(results)-1):
        theslit = results[r]

        try:
            bf = theslit["bottom"]
            tf = theslit["top"]
        except:
            pdb.set_trace()

        for i in range(hpps[0], hpps[1]):
            top = int(np.floor(tf(i)))
            bottom = int(np.ceil(bf(i)))
            
            data[top:bottom, i] = flat[top:bottom,i]

    hdu.data = (data/flat).astype(np.float32)
    bad = np.abs(hdu.data-1.0) > 0.5
    hdu.data[bad] = 1.0
    hdu.data = hdu.data.filled(1)
    if os.path.exists(outfile):
            os.remove(outfile)
    hdu.writeto(outfile)
    info("Done.")


def save_ds9_edges(results, options):
    '''
    Create a ds9 file that saves the fit slit edge positions determined
    by find_and_fit_edges
    '''

    ds9 = ''

    W = Options.npix
    delt = Options.npix/30.
    
    S = 1
    for i in range(len(results) - 1):
        res = results[i]

        top = res["top"]
        bottom = res["bottom"]
        for i in np.arange(W/delt):
            x = delt * i
            sx = x + 1
            ex = x + delt + 1

            sy = top(sx) 
            ey = top(ex) 

            smid = (top(sx) - bottom(sx)) / 2. + bottom(sx)
            emid = (top(ex) - bottom(ex)) / 2. + bottom(sx)

            # three quarter point
            stq = (top(sx) - bottom(sx)) * 3./4. + bottom(sx)
            etq = (top(ex) - bottom(ex)) * 3./4. + bottom(sx)
            # one quarter point
            soq = (top(sx) - bottom(sx)) * 1./4. + bottom(sx)
            eoq = (top(ex) - bottom(ex)) * 1./4. + bottom(sx)

            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=red\n" % (sx, smid, ex, emid)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=red\n" % (sx, stq, ex, etq)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=red\n" % (sx, soq, ex, eoq)


            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)


            if i == W/2:
                    ds9 += " # text={S%2.0i (%s)}" % (S, 
                                    res["Target_Name"])

            ds9 += "\n"

            sy = bottom(sx) + 1
            ey = bottom(ex) + 1
            if i == 10: txt=res["Target_Name"]
            else: txt=""

            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=blue text={%s}\n" % (sx, sy, ex, ey, txt)

        # Vertical line indicating half power points
        try:
            hpps = res["hpps"]
            sx = hpps[0] ; ex = hpps[0]
            sy = bottom(sx) ; ey = top(sx)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)

            sx = hpps[1] ; ex = hpps[1]
            sy = bottom(sx) ; ey = top(sx)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)
        except:
            continue
        
    band = results[-1]["band"]
    fn = "slit-edges_%s.reg" % band
    try:
            f = open(fn,'w')
            f.write(ds9)
            f.close()
    except IOError:
            error("IO Error")
            raise
    except:
            raise

def find_edge_pair(data, y, roi_width, edgeThreshold=450):
    '''
    find_edge_pair finds the edge of a slit pair in a flat

    data[2048x2048]: a well illuminated flat field [DN]
    y: guess of slit edge position [pix]

    Keywords:
    
    edgeThreshold: the pixel value below which we should ignore using
    to calculate edges.
    
    Moves along the edge of a slit image
            - At each location along the slit edge, determines
            the position of the demarcations between two slits

    Outputs:
    xposs []: Array of x positions along the slit edge [pix]
    yposs []: The fitted y positions of the "top" edge of the slit [pix]
    widths []: The fitted delta from the top edge of the bottom [pix]
    scatters []: The amount of light between slits


    The procedure is as follows
    1: starting from a guess spatial position (parameter y), march
        along the spectral direction in some chunk of pixels
    2: At each spectral location, construct a cross cut across the
        spatial direction; select_roi is used for this.
    3: Fit a two-sided error function Fit.residual_disjoint_pair
        on the vertical cross cut derived in step 2.
    4: If the fit fails, store it in the missing list
        - else if the top fit is good, store the top values in top vector
        - else if the bottom fit is good, store the bottom values in bottom
          vector.
    5: In the vertical cross-cut, there is a minimum value. This minimum
        value is stored as a measure of scattered light.

    Another procedure is used to fit polynomials to these fitted values.
    '''

    def select_roi(data, roi_width):
        v = data[int(y-roi_width):int(y+roi_width), int(xp)-2:int(xp)+2]
        v = np.median(v, axis=1) # Axis = 1 is spatial direction

        return v



    xposs_top = []
    yposs_top = []
    xposs_top_missing = []

    xposs_bot = []
    yposs_bot = []
    xposs_bot_missing = []
    yposs_bot_scatters = []

    #1 
    rng = np.linspace(10, 2040, 50).astype(np.int)
    for i in rng:
        xp = i
        #2
        v = select_roi(data, roi_width)
        xs = np.arange(len(v))

        # Modified from 450 as the hard coded threshold to one that
        # can be controlled by a keyword
        if (np.median(v) < edgeThreshold):
            xposs_top_missing.append(xp)
            xposs_bot_missing.append(xp)
            continue

        #3
        ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
        fit_ok = 0 < ff[4] < 4

        if fit_ok:
            (sigma, offset, mult1, mult2, add, width) = ff[0]

            xposs_top.append(xp)
            yposs_top.append(y - roi_width + offset + width)

            xposs_bot.append(xp)
            yposs_bot.append(y - roi_width + offset)

            between = offset + width/2
            if 0 < between < len(v)-1:
                start = int(np.max([0, between-2]))
                stop = int(np.min([len(v),between+2]))
                yposs_bot_scatters.append(np.min(v[start:stop])) # 5

                if False:
                    pl.figure(2)
                    pl.clf()
                    tmppix = np.arange(y-roi_width, y+roi_width)
                    tmpx = np.arange(len(v))
                    pl.axvline(y - roi_width + offset + width, color='red')
                    pl.axvline(y - roi_width + offset, color='red')
                    pl.scatter(tmppix, v)
                    pl.plot(tmppix, Fit.fit_disjoint_pair(ff[0], tmpx))
                    pl.axhline(yposs_bot_scatters[-1])
                    pl.draw()

            else:
                yposs_bot_scatters.append(np.nan)

        else:
            xposs_bot_missing.append(xp)
            xposs_top_missing.append(xp)
            info("Skipping wavelength pixel): %i" % (xp))

    
    return map(np.array, (xposs_bot, xposs_bot_missing, yposs_bot, xposs_top,
        xposs_top_missing, yposs_top, yposs_bot_scatters))

def fit_edge_poly(xposs, xposs_missing, yposs, order):
    '''
    fit_edge_poly fits a polynomial to the measured slit edges.
    This polynomial is used to extract spectra.

    fit_edge_poly computes a parabola, and fills in missing data with a 
    parabola

    input-
    xposs, yposs [N]: The x and y positions of the slit edge [pix]
    order: the polynomial order
    '''

    # First fit low order polynomial to fill in missing data
    fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, 2))

    xposs = np.append(xposs, xposs_missing)
    yposs = np.append(yposs, fun(xposs_missing))

    # Remove any fits that deviate wildly from the 2nd order polynomial
    ok = np.abs(yposs - fun(xposs)) < 1
    if not ok.any():
            error("Flat is not well illuminated? Cannot find edges")
            raise Exception("Flat is not well illuminated? Cannot find edges")

    # Now refit to user requested order
    fun = np.poly1d(Fit.polyfit_clip(xposs[ok], yposs[ok], order))
    yposs_ok = yposs[ok]
    res = fun(xposs[ok]) - yposs[ok]
    sd = np.std(res)
    ok = np.abs(res) < 2*sd


    # Check to see if the slit edge funciton is sane, 
    # if it's not, then we fix it.
    pix = np.arange(2048)
    V = fun(pix)
    if np.abs(V.max() - V.min()) > 10:
        info ("Forcing a horizontal slit edge")
        print("Forcing a horizontal slit edge")
        tmp = yposs_ok[ok]
        fun = np.poly1d(np.median(tmp))
        #fun = np.poly1d(np.median(yposs[ok]))


    return (fun, res, sd, ok)

def find_long2pos_edges(data, header, bs, options, edgeThreshold=450,longslit=None):


    y = 2034
    DY = 44.25


    toc = 0
    ssl = bs.ssl

    slits = []

    top = [0., np.float(Options.npix)]

    start_slit_num = int(bs.msl[0]['Slit_Number'])-1
    if start_slit_num > 0:
        y -= DY * start_slit_num
    # if the mask is a long slit, the default y value will be wrong. Set instead to be the middle
    if bs.long_slit:
        try:
            y=longslit["yrange"][1]
        except:
            error ("Longslit reduction mode is specified, but the row position has not been specified. Defaulting to "+str(y))
            print("Longslit reduction mode is specified, but the row position has not been specified. Defaulting to "+str(y))


    # Count and check that the # of objects in the SSL matches that of the MSL
    # This is purely a safety check
    numslits = np.zeros(len(ssl))
    for i in range(len(ssl)):
        slit = ssl[i]
        M = np.where(slit["Target_Name"] == bs.msl["Target_in_Slit"])

        numslits[i] = len(M[0])
    numslits = np.array(numslits)
    info("Number of slits allocated for this longslit: "+str(np.sum(numslits)))

    # now begin steps outline above
    results = []

    for slit in [0,1]:
        
        result = {}

        result["Target_Name"] = ssl[slit]["Target_Name"]
        # 1 Defines a polynomial of degree 0, which is a constant, with the value of the top of the slit
        result["top"] = np.poly1d([longslit["yrange"][slit][1]]) # 1
        topfun = np.poly1d([longslit["yrange"][slit][1]]) # this is a constant funtion with c=top of the slit # 1
        botfun = np.poly1d([longslit["yrange"][slit][0]]) # this is a constant funtion with c=bottom of the slit # 0

        # xposs_top_this = [10 110 210 .... 1810 1910]
        xposs_top = np.arange(10,2000,100)
        xposs_bot = np.arange(10,2000,100)
        # yposs_top_this = [1104 1104 ... 1104 1104], it's the constant polynomium calculated at the X positions
        yposs_top = topfun(xposs_top)
        yposs_bot = botfun(xposs_bot)

        ''' Deal with the current slit '''
        target=slit
        hpps = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(target+1)[0], header, bs)

        ok = np.where((xposs_top > hpps[0]) & (xposs_top < hpps[1]))

        xposs_bot = xposs_bot[ok]
        yposs_bot = yposs_bot[ok]
        xposs_top = xposs_top[ok]
        yposs_top = yposs_top[ok]

        if len(xposs_bot) == 0:
            error ("The slit edges specifications appear to be incorrect.")
            raise Exception("The slit edges specifications appear to be incorrect.")

        # bot is the polynomium that defines the shape of the bottom of the slit. In this case, we set it to a constant.
        bot = botfun.c.copy() 
        top = topfun.c.copy()


        #4
        result = {}
        result["Target_Name"] = ssl[target]["Target_Name"]
        result["xposs_top"] = xposs_top
        result["yposs_top"] = yposs_top
        result["xposs_bot"] = xposs_bot
        result["yposs_bot"] = yposs_bot
        result["top"] = np.poly1d(top)
        result["bottom"] = np.poly1d(bot)
        result["hpps"] = hpps
        result["ok"] = ok
        results.append(result)
        #print("And the top is"+str(result["top"]))

    results.append({"version": options["version"]})
    return results


def find_longslit_edges(data, header, bs, options, edgeThreshold=450,longslit=None):


    y = 2034
    DY = 44.25


    toc = 0
    ssl = bs.ssl

    slits = []

    top = [0., np.float(Options.npix)]

    start_slit_num = int(bs.msl[0]['Slit_Number'])-1
    if start_slit_num > 0:
        y -= DY * start_slit_num
    # if the mask is a long slit, the default y value will be wrong. Set instead to be the middle
    if bs.long_slit:
        try:
            y=longslit["yrange"][1]
        except:
            error("Longslit reduction mode is specified, but the row position has not been specified. Defaulting to "+str(y))
            print("Longslit reduction mode is specified, but the row position has not been specified. Defaulting to "+str(y))


    # Count and check that the # of objects in the SSL matches that of the MSL
    # This is purely a safety check
    numslits = np.zeros(len(ssl))
    for i in range(len(ssl)):
        slit = ssl[i]
        M = np.where(slit["Target_Name"] == bs.msl["Target_in_Slit"])

        numslits[i] = len(M[0])
    numslits = np.array(numslits)
    info("Number of slits allocated for this longslit: "+str(np.sum(numslits)))

    # now begin steps outline above
    results = []
    result = {}
    result["Target_Name"] = ssl[0]["Target_Name"]

    # 1 Defines a polynomial of degree 0, which is a constant, with the value of the top of the slit
    result["top"] = np.poly1d([longslit["yrange"][1]])
    
    topfun = np.poly1d([longslit["yrange"][1]]) # this is a constant funtion with c=top of the slit
    botfun = np.poly1d([longslit["yrange"][0]]) # this is a constant funtion with c=bottom of the slit

    # xposs_top_this = [10 110 210 .... 1810 1910]
    xposs_top = np.arange(10,2000,100)
    xposs_bot = np.arange(10,2000,100)
    # yposs_top_this = [1104 1104 ... 1104 1104], it's the constant polynomium calculated at the X positions
    yposs_top = topfun(xposs_top)
    yposs_bot = botfun(xposs_bot)

    
    ''' Deal with the current slit '''
    target=0
    hpps = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(target+1)[0], header, bs)

    ok = np.where((xposs_top > hpps[0]) & (xposs_top < hpps[1]))

    xposs_bot = xposs_bot[ok]
    yposs_bot = yposs_bot[ok]
    xposs_top = xposs_top[ok]
    yposs_top = yposs_top[ok]

    if len(xposs_bot) == 0:
        error ("The slit edges specifications appear to be incorrect.")
        raise Exception("The slit edges specifications appear to be incorrect.")

    # bot is the polynomium that defines the shape of the bottom of the slit. In this case, we set it to a constant.
    bot = botfun.c.copy() 
    top = topfun.c.copy()


    #4
    result = {}
    result["Target_Name"] = ssl[target]["Target_Name"]
    result["xposs_top"] = xposs_top
    result["yposs_top"] = yposs_top
    result["xposs_bot"] = xposs_bot
    result["yposs_bot"] = yposs_bot
    result["top"] = np.poly1d(top)
    result["bottom"] = np.poly1d(bot)
    result["hpps"] = hpps
    result["ok"] = ok
    results.append(result)

    results.append({"version": options["version"]})

    return results

def find_and_fit_edges(data, header, bs, options,edgeThreshold=450):
    '''
    Given a flat field image, find_and_fit_edges determines the position
    of all slits.

    The function works by starting with a guess at the location for a slit
    edge in the spatial direction(options["first-slit-edge"]). 
    
    Starting from the guess, find_edge_pair works out in either direction, 
    measuring the position of the (e.g.) bottom of slit 1 and top of slit 2:


    ------ pixel y value = 2048

    Slit 1 data

    ------ (bottom)
    deadband
    ------ (top)

    Slit N pixel data ....

    ------- (bottom) pixel = 0

    --------------------------------> Spectral direction


    1. At the top of the flat, the slit edge is defined to be a pixel value
    2. The code guesses the position of the bottom of the slit, and runs
            find_edge_pair to measure slit edge locations.
    3. A low-order polynomial is fit to the edge locations with
            fit_edge_poly
    4. The top and bottom of the current slit, is stored into the
            result list.
    5. The top of the next slit is stored temporarily for the next
            iteration of the for loop.
    6. At the bottom of the flat, the slit edge is defined to be pixel 4.


    options:
    options["edge-order"] -- The order of the polynomial [pixels] edge.
    options["edge-fit-width"] -- The length [pixels] of the edge to 
            fit over

    '''

    # TODO: move hardcoded values into Options.py
    # y is the location to start
    y = 2034
    DY = 44.25

    toc = 0
    ssl = bs.ssl

    slits = []

    top = [0., np.float(Options.npix)]

    start_slit_num = int(bs.msl[0]['Slit_Number'])-1
    if start_slit_num > 0:
        y -= DY * start_slit_num

    # Count and check that the # of objects in the SSL matches that of the MSL
    # This is purely a safety check
    numslits = np.zeros(len(ssl))
    for i in range(len(ssl)):
        slit = ssl[i]
        M = np.where(slit["Target_Name"] == bs.msl["Target_in_Slit"])

        numslits[i] = len(M[0])
    numslits = np.array(numslits)


    if (np.sum(numslits) != CSU.numslits) and (not bs.long_slit) and (not bs.long2pos_slit):
        error ("The number of allocated CSU slits (%i) does not match "
                " the number of possible slits (%i)." % (np.sum(numslits),
                    CSU.numslits))
        raise Exception("The number of allocated CSU slits (%i) does not match "
                " the number of possible slits (%i)." % (np.sum(numslits),
                    CSU.numslits))

    # if the mask is a long slit, the default y value will be wrong. Set instead to be the middle
    if bs.long_slit:
        y = 1104
        
    # now begin steps outline above
    results = []
    result = {}

    result["Target_Name"] = ssl[0]["Target_Name"]

    # 1
    result["top"] = np.poly1d([y])

    ''' Nomenclature here is confusing:
        
        ----- Edge  -- Top of current slit, bottom of prev slit
        . o ' Data
        ===== Data
        .;.;' Data
        ----- Edge  -- Bottom of current slit, top of next slit
    '''

    topfun = np.poly1d([y])
    xposs_top_this = np.arange(10,2000,100)
    yposs_top_this = topfun(xposs_top_this)

    initial_edges = np.array([2034], dtype=np.int)
    edge = 2034

    # build an array of values containing the lower edge of the slits
    
    for target in range(len(ssl)):        
    # target is the slit number
        edge -= DY * numslits[target]
        initial_edges=np.append(initial_edges,int(edge))

    # collapse the 2d flat along the walenegth axis to build a spatial profile of the slits
    vertical_profile = np.mean(data, axis=1)

    # build an array containing the spatial positions of the slit centers, basically the mid point between the expected
    # top and bottom values of the slit pixels
    spatial_centers = np.array([], dtype=np.int)
    for k in np.arange(0,len(initial_edges)-1):
        spatial_centers = np.append(spatial_centers,(initial_edges[k]+initial_edges[k+1])/2)
    #slit_values=np.array([])
    #for k in np.arange(0, len(spatial_centers)):
    #    slit_values = np.append(slit_values,np.mean(vertical_profile[spatial_centers[k]-3:spatial_centers[k]+3]))
    
    for target in range(len(ssl)):

        y -= DY * numslits[target]
        y = max(y, 1)
        # select a 6 pixel wide section of the vertical profile around the slit center
        threshold_area = vertical_profile[spatial_centers[target]-3:spatial_centers[target]+3]
        # uses 80% of the ADU counts in the threshold area to estimate the threshold to use in defining the slits
        edgeThreshold = np.mean(threshold_area)*0.8
        #if edgeThreshold > 450:
        #    edgeThreshold = 450
        
        info("[%2.2i] Finding Slit Edges for %s ending at %4.0i. Slit "
                "composed of %i CSU slits" % ( target,
                    ssl[target]["Target_Name"], y, numslits[target]))
        info("[%2.2i] Threshold used is %.1f" % (target,edgeThreshold))

        ''' First deal with the current slit '''
        hpps = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(target+1)[0], header, bs)

        if y == 1:
            xposs_bot = [1024]
            xposs_bot_missing = []
            yposs_bot = [4.25]
            botfun = np.poly1d(yposs_bot)
            ok = np.where((xposs_bot > hpps[0]) & (xposs_bot < hpps[1]))
        else:
            (xposs_top_next, xposs_top_next_missing, yposs_top_next, xposs_bot,
                xposs_bot_missing, yposs_bot, scatter_bot_this) = find_edge_pair(
                    data, y, options["edge-fit-width"],edgeThreshold=edgeThreshold)

            ok = np.where((xposs_bot > hpps[0]) & (xposs_bot < hpps[1]))
            ok2 = np.where((xposs_bot_missing > hpps[0]) & (xposs_bot_missing <
                hpps[1]))
            xposs_bot = xposs_bot[ok]
            xposs_bot_missing = xposs_bot_missing[ok2]
            yposs_bot = yposs_bot[ok]
            if len(xposs_bot) == 0:
                botfun = np.poly1d(y-DY)
            else:
                (botfun, bot_res, botsd, botok) =  fit_edge_poly(xposs_bot,
                         xposs_bot_missing, yposs_bot, options["edge-order"])


        bot = botfun.c.copy() 
        top = topfun.c.copy()

        #4
        result = {}
        result["Target_Name"] = ssl[target]["Target_Name"]
        result["xposs_top"] = xposs_top_this
        result["yposs_top"] = yposs_top_this
        result["xposs_bot"] = xposs_bot
        result["yposs_bot"] = yposs_bot
        result["top"] = np.poly1d(top)
        result["bottom"] = np.poly1d(bot)
        result["hpps"] = hpps
        result["ok"] = ok
        results.append(result)

        #5
        if y == 1:
            break
            

        next = target + 2
        if next > len(ssl): next = len(ssl)
        hpps_next = Wavelength.estimate_half_power_points(
                bs.scislit_to_csuslit(next)[0],
                    header, bs)

        ok = np.where((xposs_top_next > hpps_next[0]) & (xposs_top_next <
            hpps_next[1]))
        ok2 = np.where((xposs_top_next_missing > hpps_next[0]) &
            (xposs_top_next_missing < hpps_next[1]))

        xposs_top_next = xposs_top_next[ok]
        xposs_top_next_missing = xposs_top_next_missing[ok2]
        yposs_top_next = yposs_top_next[ok]

        if len(xposs_top_next) == 0:
            topfun = np.poly1d(y)
        else:
            (topfun, topres, topsd, ok) = fit_edge_poly(xposs_top_next,
                xposs_top_next_missing, yposs_top_next, options["edge-order"])

        xposs_top_this = xposs_top_next
        xposs_top_this_missing = xposs_top_next_missing
        yposs_top_this = yposs_top_next

    results.append({"version": options["version"]})

    return results

class FitCheck:
    flat = None
    bs = None
    edges = None
    cutout = None # Is at most 2 edges x 46 slits x 11 pix or 1112 pixels
    
    def __init__(self, maskname, bandname, options, fig):

        self.fig = fig
        self.flat = IO.read_drpfits(maskname, "combflat_2d_%s.fits" % bandname,
                options)

        self.edges, meta = IO.load_edges(maskname, bandname, options)

        self.edgeno=2
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)

        self.draw()


    def draw(self):

        print(self.edgeno)

        pos = 0
        dy = 8
        edgeno = self.edgeno
        edge = self.edges[edgeno]
        edgeprev = self.edges[edgeno-1]
        p = np.round(edge["top"](1024))
        top = min(p+2*dy, 2048)
        bot = min(p-2*dy, 2048)
        self.cutout = self.flat[1][bot:top,:].copy()

        pl.figure(1)
        pl.clf()
        start = 0
        dy = 512
        for i in range(2048/dy):
            pl.subplot(2048/dy,1,i+1)
            pl.xlim(start, start+dy)

            if i == 0: pl.title("edge %i] %s|%s" % (edgeno,
                edgeprev["Target_Name"], edge["Target_Name"]))


            pl.subplots_adjust(left=.07,right=.99,bottom=.05,top=.95)

            pl.imshow(self.flat[1][bot:top,start:start+dy], extent=(start,
                start+dy, bot, top), cmap='Greys', vmin=2000, vmax=6000)

            pix = np.arange(start, start+dy)
            pl.plot(pix, edge["top"](pix), 'r', linewidth=1)
            pl.plot(pix, edgeprev["bottom"](pix), 'r', linewidth=1)
            pl.plot(edge["xposs_top"], edge["yposs_top"], 'o')
            pl.plot(edgeprev["xposs_bot"], edgeprev["yposs_bot"], 'o')


            hpp = edge["hpps"]
            pl.axvline(hpp[0],ymax=.5, color='blue', linewidth=5)
            pl.axvline(hpp[1],ymax=.5, color='red', linewidth=5)

            hpp = edgeprev["hpps"]
            pl.axvline(hpp[0],ymin=.5,color='blue', linewidth=5)
            pl.axvline(hpp[1],ymin=.5,color='red', linewidth=5)


            if False:
                L = top-bot
                Lx = len(edge["xposs"])
                for i in range(Lx):
                    xp = edge["xposs"][i]
                    frac1 = (edge["top"](xp)-bot-1)/L
                    pl.axvline(xp,ymin=frac1)

                for xp in edgeprev["xposs"]: 
                    frac2 = (edgeprev["bottom"](xp)-bot)/L
                    pl.axvline(xp,ymax=frac2)

            start += dy

    def __call__(self, event):
        kp = event.key
        x = event.xdata
        y = event.ydata

        print(kp)

        if kp == 'n': 
            self.edgeno += 1

            if self.edgeno > len(self.edges):
                self.edgeno = len(self.edges)
                print("done")
            else:
                self.draw()

        if kp == 'p': 
            self.edgeno -= 1
            if self.edgeno < 2:
                self.edgeno = 2
                print("Beginning" )
            else:
                self.draw()

    

class TestFlatsFunctions(unittest.TestCase):
    def setUp(self):
            pass

    def test_trace_edge(self):
            (header, data1, targs, ssl, msl, asl) = \
                            IO.readfits_all("/users/npk/desktop/c9/m110326_3242.fits")
            data = data1

            ssl = ssl[ssl["Slit_Number"] != ' ']
            numslits = np.round(np.array(ssl["Slit_length"], 
                    dtype=np.float) / 7.02)

            for i in range(len(ssl)):
                    print(ssl[i]["Target_Name"], numslits[i])

if __name__ == '__main__':
    unittest.main()

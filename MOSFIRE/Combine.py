
import os
import pdb

import numpy as np
try:
    import pyfits as pf
except:
    from astropy.io import fits as pf


import MOSFIRE

from MOSFIRE import Background, IO, Wavelength


def imcombine(filelist, maskname, fname, options, sum_type):
    ''' combine the images in file list into fname.

    Sum type:
        rate -- filelist is in cnt/s
        ivar-rate -- filelist is in s/cnt
        snr-rate -- filelist is in SNR
    '''

    ARR = None
    hdr = None
    i = 1

    itime = 0

    for file in filelist:

        this_hdr, img = IO.readfits(file)
        cards = this_hdr.ascardlist()

        thisitime = this_hdr['truitime']
        itime += thisitime

        if ARR is None: ARR = np.zeros(img.shape)

        if sum_type == 'rate': ARR += img * thisitime
        if sum_type == 'ivar-rate': ARR += thisitime/img
        if sum_type == 'snr-rate': ARR  += img * thisitime


        if hdr is None:
            hdr = this_hdr
        hdr.update("fno%2.2i" % i, file, "--")
        for card in cards:
            key, value, comment = (card.key, card.value, card.comment)

            if hdr.has_key(key) and hdr[key] != value:
                key = key + ("_img%2.2i" % i)

            if len(key) > 8: key = 'HIERARCH ' + key

            try: hdr.update(key, value, comment)
            except ValueError: pass
    
    hdr.update('itime', itime, 'Itime for %i rectified images' % len(filelist))
    if sum_type == 'rate': ARR /= itime 
    if sum_type == 'ivar-rate': ARR = itime/ARR
    if sum_type == 'snr-rate': ARR /= itime

    IO.writefits(ARR, maskname, fname, options, header=hdr, overwrite=True,
            lossy_compress=True)

def get_path(a):
    if os.path.exists(a): return a

    a += ".gz"
    if os.path.exists(a): return a

    raise Exception("No such path: %s" % a)

def gz(name):
    if name[-2:] == 'gz': return '.gz'
    else: return ''

def stack_rectified(wavenames, maskname, band, wavops):
    N = len(wavenames)
    lamnames = []
    suffixs = []

    for i in xrange(N):
        lamnames.append( Wavelength.filelist_to_wavename(wavenames[i], band,
            maskname, wavops).rstrip(".fits"))
        suffixs.append(lamnames[i].lstrip("wave_stack_%s_" % band))
    
    path = os.path.join(wavops["outdir"], maskname)


    recs = []
    ivars = []
    sns = []


    try:
        ls = [get_path(os.path.join(path, "eps_%s_%s_%s.fits") % (maskname,
            suffix, band)) for suffix in suffixs]
        imcombine(ls, maskname, "eps_%s_%s.fits" % (maskname, band), 
                wavops, 'rate')
    except:
        pass

    try:
        ls = [get_path(os.path.join(path, "ivars_%s_%s_%s.fits") % (maskname, 
            suffix, band)) for suffix in suffixs]
        imcombine(ls, maskname, "ivars_%s_%s.fits" % (maskname, 
            band), wavops, 'ivar-rate')
    except:
        pass

    try:
        ls = [get_path(os.path.join(path, "snrs_%s_%s_%s.fits") % (maskname,
            suffix, band)) for suffix in suffixs]
        imcombine(ls, maskname, "snrs_%s_%s.fits" % (maskname, band), 
            wavops, 'snr-rate')
    except:
        pass




def stack_slits(wavenames, maskname, band, wavops):
    pass

def stack_full(wavenames, maskname, band, wavops):
    pass

def stack_files(wavenames, maskname, band, wavops):
    stack_rectified(wavenames, maskname, band, wavops)
    stack_slits(wavenames, maskname, band, wavops)
    stack_full(wavenames, maskname, band, wavops)

def rename_files(wavenames, maskname, band, wavops):

    lamname = Wavelength.filelist_to_wavename(wavenames[0], band, maskname,
            wavops).rstrip(".fits")
    
    suffix = lamname.lstrip("wave_stack_%s_" % band)

    path = os.path.join(wavops["outdir"], maskname)

    fnames = ["rectified_%s%s.fits", "rectified_ivar_%s%s.fits",
            "rectified_sn_%s%s.fits"]

    for fname in fnames:
        try:
            a = get_path(os.path.join(path, fname % (band, "_" + suffix)))
            b = os.path.join(path, fname % (band, "")) + gz(a)
            os.rename(a, b)
        except: 
            print "Ignoring renaming of: ", fname
            pass
    

    edges = IO.load_edges(maskname, band, wavops)
    n_slits = len(edges[0])

    for i in xrange(1, n_slits+1):
        S = "S%2.2i" % (i)
        a = get_path(os.path.join(path, 
            "eps_%s_%s_%s.fits" % (band, suffix, S)))
        
        a_h = pf.open(a)[0].header
        obj = a_h['object']

        b = os.path.join(path, "%s_%s_%s_eps.fits" % (maskname, band, obj)) + \
            gz(a)
        os.rename(a,b)

        a = get_path(os.path.join(path, 
            "ivar_%s_%s_%s.fits" % (band, suffix, S)))
        a_h = pf.open(a)[0].header
        obj = a_h['object']

        b = os.path.join(path, "%s_%s_%s_ivar.fits" % (maskname, band, obj)) + \
            gz(a)
        os.rename(a,b)

    a = get_path(os.path.join(path,
        "eps_%s_%s_%s.fits" % (maskname, suffix, band)))
    b = os.path.join(path,
        "%s_%s_eps.fits" % (maskname, band)) + gz(a)
    os.rename(a,b)

    a = get_path(os.path.join(path,
        "snrs_%s_%s_%s.fits" % (maskname, suffix, band)))
    b = os.path.join(path,
        "%s_%s_snrs.fits" % (maskname, band)) + gz(a)
    os.rename(a, b)

    a = get_path(os.path.join(path,
        "ivars_%s_%s_%s.fits" % (maskname, suffix, band)))
    b = os.path.join(path,
        "%s_%s_ivars.fits" % (maskname, band)) + gz(a)
    os.rename(a, b)

def handle_combine(wavenames, maskname, band, wavops):

    N = len(wavenames)
    assert(N > 0)

    print "Starting"

    if N == 1:  rename_files(wavenames, maskname, band, wavops)
    if N > 1:   stack_files(wavenames, maskname, band, wavops)


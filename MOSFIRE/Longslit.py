
# MOSFIRE Longslit Reductions
# 5 Aug 2012
# Nick Konidaris

import os
import pdb
import numpy as np
import scipy

from MOSFIRE import Detector, IO, Filters, Wavelength

def rectify(dname, lamdat, A, B, maskname, band, wavoptions, 
        longoptions):

    header, data = IO.readfits(dname)
    raw_img = data * Detector.gain / header['TRUITIME']

    dlam = Wavelength.grating_results(band)
    hpp = np.array(Filters.hpp[band]) 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)

    rectified = np.zeros((2048, len(ll_fid)))

    from scipy.interpolate import interp1d

    for i in xrange(2048):
        ll = lamdat[i,:]
        ss = raw_img[i,:]
        ok = np.isfinite(ll) & np.isfinite(ss) & (ll < hpp[1]) & (ll >
                hpp[0])

        if len(np.where(ok)[0]) < 30:
            continue

        f = interp1d(ll[ok], ss[ok], bounds_error=False)
        rectified[i,:] = f(ll_fid)

    header["wat0_001"] = "system=world"
    header["wat1_001"] = "wtype=linear"
    header["wat2_001"] = "wtype=linear"
    header["dispaxis"] = 1
    header["dclog1"] = "Transform"
    header["dc-flag"] = 0
    header["ctype1"] = "AWAV"
    header["cunit1"] = "Angstrom"
    header["crval1"] = ll_fid[0]
    header["crval2"] = 0
    header["crpix1"] = 1
    header["crpix2"] = 1
    header["cdelt1"] = 1
    header["cdelt2"] = 1
    header["cname1"] = "angstrom"
    header["cname2"] = "pixel"
    header["cd1_1"] = dlam
    header["cd1_2"] = 0
    header["cd2_1"] = 0
    header["cd2_2"] = 1


    header["object"] = "rectified [eps]"
    IO.writefits(rectified, maskname, "rectified_%s" % (dname), 
        wavoptions, header=header, overwrite=True, lossy_compress=True)


def imdiff(A, B, maskname, band, header, options):
    s = "[0]"

    targname = A[1]["targname"].rstrip(" ")
    if targname == "":
        objname = A[1]["object"].replace(" ", "_")
    else:
        objname = targname.replace(" ", "_")

    operand1 = A[0] + '[0]'
    operand2 = B[0] + '[0]'

    imnumA = A[0].split('_')[-1].rstrip(".fits")
    imnumB = B[0].split('_')[-1].rstrip(".fits")

    dname = "{0}_{1}_{2}_{3}-{4}_{5}-{6}.fits".format(maskname, objname, band,
        A[1]["frameid"], B[1]["frameid"], imnumA, imnumB)

    try: os.remove(dname)
    except:pass
    print("Data Diff {0}-{1}".format(operand1,operand2))
    IO.imarith(operand1, '-', operand2, dname)

    ''' Now handle variance '''
    numreads = header["READS0"]
    RN_adu = Detector.RN / np.sqrt(numreads) / Detector.gain
    varname = "var_{0}_{1}_{2}_{3}+{4}_{5}+{6}.fits".format(maskname, objname, band,
        A[1]["frameid"], B[1]["frameid"], imnumA, imnumB)

    
    print("Var Sum {0}+{1}".format(operand1,operand2))
    IO.imarith(operand1, '+', operand2, "tmp_" + varname)
    try: os.remove(varname)
    except: pass
    print("Var add RN {0}+{1}".format(operand1,RN_adu**2))
    IO.imarith("tmp_" + varname, '+', RN_adu**2, varname)

    try: os.remove("tmp_" + varname)
    except: pass


    return dname, varname

def apply_flat(scifilename, maskname, band):
    ''' Divides the contents of scifilename by the flat field and
        overwrites scifilename with the same file divided by the flat

        Args:
            scifilename: Path to science file name.
            maskname: The mask name
            band: The filter bands

        Results:
            Overwrites scifilename where the data contents of the file
                are divided by the pixel flat
    '''

    
    pixelflat_file = "pixelflat_2d_{0}.fits".format(band)
    flat = readfits(pixelflat_file, use_bpm=True)
    flat_data = flat[1].filled(1.0)

    header, data = IO.readfits(scifilename)
    
    print("Applying flat to file {0}".format(scifilename))
    IO.writefits(data/flat_data, maskname, scifilename, {}, header=header,
        overwrite=True)



def go(maskname,
        band,
        filenames,
        wavefile,
        wavoptions,
        longoptions,
        use_flat=False):

    '''
    The go command is the main entry point into this module.

    Inputs:
        maskname: String of the mask name
        band: String of 'Y', 'J', 'H', or 'K'
        filenames: List of filenames to reduce
        wavefile: String of path to FITS file with the wavelength solution
        wavoptions: The Wavelength Options dictionary
        longoptions: Dictionary containing:
            {'yrange': The pixel range to extract over
            'row_position': The row to solve the initial wavelength solution on}
        use_flat: Boolean False [default] means to use no flat field
            Boolean True means to divide by the pixelflat
    '''
    wavename = Wavelength.filelist_to_wavename(filenames, band, maskname,
            wavoptions).rstrip(".fits")

    print("Wavefile: {0}".format(wavefile))
    lamhdr, lamdat = IO.readfits(wavefile)

    positions = []
    objname = None
    for listfile in filenames:
        fnames = IO.list_file_to_strings(listfile)
        if len(fnames) != 1:
            raise Exception("I currently expect only one file per position. Remove multiple entries and try again")

        header, data, bs = IO.readmosfits(fnames[0], wavoptions)

        if objname is None:
            objname = header["object"]

        if objname != header["object"]:
            print ("Trying to combine longslit stack of object {0} " 
                    "with object {1}".format(objname, header["object"]))

        print("{0:18s} {1:30s} {2:2s} {3:4.1f}".format(file, header["object"],
            header["frameid"], header["yoffset"]))

        positions.append([fnames[0], header, data, bs])

    print("{0:2g} nod positions found. Producing stacked difference" \
           " image.".format(len(positions)))

    for i in xrange(len(positions)-1):
        A = positions[i]
        B = positions[i+1]

        print("----------- -----".format(A,B))
        
        dname, varname = imdiff(A, B, maskname, band, header, wavoptions)
        if use_flat:
            apply_flat(dname, maskname, band)
            apply_flat(varname, maskname, band)

        rectify(dname, lamdat, A, B, maskname, band, wavoptions,
                longoptions)
        rectify(varname, lamdat, A, B, maskname, band, wavoptions,
                longoptions)
        print dname

    dname, vname = imdiff(B, A, maskname, band, header, wavoptions)
    if use_flat:
        apply_flat(dname, maskname, band)
        apply_flat(vname, maskname, band)
    rectify(dname, lamdat, B, A, maskname, band, wavoptions,
            longoptions)
    rectify(vname, lamdat, B, A, maskname, band, wavoptions,
            longoptions)
    
    if False:
        fname = os.path.join(path, wavename + ".fits")
        B = IO.readfits(fname)
        B = [fname, B[0], B[1]]
        for i in xrange(len(positions)):
            A = positions[i]
            imdiff(A, B, maskname, band, wavoptions)
            rectify(path, dname, lamdat, A, B, maskname, band, wavoptions,
                longoptions)
        imdiff(B, A, maskname, band, wavoptions)
        rectify(path, dname, lamdat, B, A, maskname, band, wavoptions,
                longoptions)



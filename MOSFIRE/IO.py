'''
MOSFIRE Input/Output Utility Code
Written March 2, 2011 by npk

Provides tools to read fits files and parse their headers.
'''

import pyfits as pf
import numpy as np
import unittest
import warnings



import os
import pdb

import MOSFIRE
import CSU
import Options

theBPM = None # the Bad pixel mask

def badpixelmask():
    global theBPM

    path = Options.path_bpm

    if theBPM is None:
        hdulist = pf.open(path)
        header = hdulist[0].header
        theBPM = hdulist[0].data

    return theBPM

def load_edges(maskname, band, options):
    ''' Load the slit edge functions. Returns (edges, metadata) '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "slit-edges_{0}.npy".format(band))

    fn = "slit-edges_{0}.npy".format(band)

    edges = np.load(fn)

    edges,meta = edges[0:-1], edges[-1]

    if meta['maskname'] != maskname:
        print("The maskname for the edge file '%s' does not match "
                "that in the edge file '%s'" % (maskname, meta['maskname']))
        print "Continuing"

    return edges, meta

def load_lambdacenter(fnum, maskname, options):
    ''' Load the wavelength coefficient functions '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_center_coeffs_{0}.npy".format(fnum))

    fn = "lambda_center_coeffs_{0}.npy".format(fnum)
    ld = np.load(fn)

    return ld

def load_lambdadata(wavename, maskname, band, options):
    ''' Load the wavelength coefficient functions '''

    if False:
        fn = os.path.join(options["outdir"], maskname,
            "lambda_coeffs_{0}.npy".format(wavename))


    fn = "lambda_coeffs_{0}.npy".format(wavename)
    ld = np.load(fn)

    return ld

def load_lambdaoutwards(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_outwards_coeffs_{0}.npy".format(fnum))

    fn = "lambda_outwards_coeffs{0}.npy".format(fnum)

    ld = np.load(fn)

    return ld

def load_lambdamodel(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_mask_coeffs_{0}.npy".format(fnum))

    fn = "lambda_mask_coeffs_{0}.npy".format(fnum)

    ld = np.load(fn)
    return ld

def load_flat(maskname, band, options):
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "pixelflat_2d_{0}.fits".format(band))

    fn = "pixelflat_2d_{0}.fits".format(band)

    return readfits(fn, use_bpm=True)


def load_lambdaslit(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_solution_{0}.fits".format(fnum))

    fn = "lambda_solution_{0}.fits".format(fnum)

    print fn

    ret = readfits(fn, options)
    if ret[0]['filter'] != band:
        raise Exception("band name mismatch")

    if ret[0]['maskname'] != maskname:
        rint("The maskname for the edge file '%s' does not match "
                "that in the edge file '%s'" % (maskname, ret[0]['maskname']))
        print "Continuing"

    
    return readfits(fn, options)

def writefits(img, maskname, fname, options, header=None, bs=None,
        overwrite=False, lossy_compress=False):
    '''Convenience wrapper to write MOSFIRE drp-friendly FITS files
    
    Args:
        img: Data array to write to disk
        maskname: Name of the science mask
        fname: Full or relative path to output file
        options: {} Unused
        header: Optional, the header to write
        bs: Optional unused
        overwrite: Force overwrite of file, default False/No.
        lossy_compress: Zero out the lowest order bits of the floats in
            order to make FITS files amenable to compression. The loss is
            at least 10 x less than 5e- which is the lowest reasonable read-
            noise value.

    Results:
        Writes a file to fname with data img and header header.

    '''

    if lossy_compress:
        hdu = pf.PrimaryHDU(floatcompress(img))
    else:
        hdu = pf.PrimaryHDU(img)

    fn = fname

    if header is None: header = {"DRPVER": (MOSFIRE.__version__, "DRP Version Date")}
    else: header["DRPVER"] = (MOSFIRE.__version__, 'DRP Version Date')

    warnings.filterwarnings('ignore')
    if header is not None:
        for k,value, comment in header.cards:
            if hdu.header.has_key(k): continue

            if k == 'COMMENT': continue
            if k == '': continue

            k = k.rstrip()
            hdu.header[k] = (value,comment)

    warnings.filterwarnings('always')
    if overwrite:
        try: 
            os.remove(fn)
            print "Removed old file '{0}'%".format(fn)
        except: pass

    print "Wrote to '%s'" % (fn)
    hdu.writeto(fn)

    if lossy_compress: os.system("gzip --force {0}".format(fn))



def readfits(path, use_bpm=False):
    '''Read a fits file from path and return a tuple of (header, data, 
    Target List, Science Slit List (SSL), Mechanical Slit List (MSL),
    Alignment Slit List (ASL)).'''

    if os.path.exists(path + ".gz"):
        path = path + ".gz"

    if not os.path.exists(path):
        raise Exception("The file at path '%s' does not exist." % path)

    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    if use_bpm:
        theBPM = badpixelmask()
        data = np.ma.masked_array(data, theBPM, fill_value=0)

    return (header, data)



def readheader(path):
    '''Reads a header (only) from a fits file'''

    return pf.getheader(path)

def read_drpfits(maskname, fname, options):
    '''Read a fits file written by the DRP'''

    if os.path.exists(fname): path = fname
    elif os.path.exists(fname + ".gz"): path = fname + ".gz"
    else: path = os.path.join(fname_to_path(fname, options), fname)

    if os.path.exists(path + ".gz"):
        path = path + ".gz"

    if not os.path.exists(path):
        raise Exception("The file at path '%s' does not exist." % path)

    hdulist = pf.open(path)
    output = []

    for hdu in hdulist:
        output.append(hdu.header)

        if "DRPVER" in hdu.header:
            itsver = hdu.header["DRPVER"]
            if itsver != MOSFIRE.__version__:
                raise Exception("The file requested '%s' uses DRP version %f "
                    "but the current DRP version is %f. There might be an "
                    "incompatibility" % (path, itsver, MOSFIRE.__version__))

        else:
            raise Exception("The file requested '%s' does not seem to be "
                    "the result of this DRP. This should never be the "
                    " case.")

        output.append(hdu.data)


    return output

def fname_to_date_tuple(fname):
    '''Take a filename like m120507_0123, return 12may07'''
    months = {"01": "jan", "02": "feb", "03": "mar", "04": "apr", "05": "may",
        "06": "jun", "07": "jul", "08": "aug", "09": "sep", "10": "oct",
        "11": "nov", "12": "dec"}

    if len(fname) != 17:
        raise Exception("The file name '%s' is not of correct length. It "
                "must be of the form mYYmmdd_nnnn.fits" % fname)
    
    try:
        fdate = fname.split("m")[1][0:6]
        yr, mn, dy = "20" + fdate[0:2], fdate[2:4], int(fdate[4:6])
        month = months[mn]
    except:
        print "Could not parse date out of file name: %s" % (fname)
    
    return yr, month, dy

def fname_to_path(fname, options):
    '''Take a filename like m120507_0123, parse date, and return full path'''

    if os.path.isabs(fname): return fname

    yr, month, dy = fname_to_date_tuple(fname)
    path = os.path.join(options["indir"], yr + month + "%2.2i" % dy)
    if not os.path.exists(os.path.join(path, fname)):
        path = os.path.join(options["indir"], yr + month + "%2.2i" % (dy-1))

        if not os.path.exists(path):
            raise Exception("Could not find file '%s' in '%s' out of parsed "
                "%s, %s, %s" % (fname,
                options["indir"], yr, month, dy))

    return path

def list_file_to_strings(fname):
    '''Read the filename in fname and convert to a series of paths.
This emulates IRAF's @file system. However, in addtion, the first line of the file
can be an absolute path. Example:
list.txt
/path/to/files
file1
file2
file3

returns ['/path/to/files/file1', '/path/to/files/file2', '/path/to/files/file3']

whereas
list.txt
file1
file2
file3

returns ['file1', 'file2', 'file3']
'''

    filelist = fname
    if type(fname) == str:
        filelist = [fname]

    if len(fname) == 0:
        return []

    if fname[0][-5:] == '.fits':
        return fname

    output = []

    for fname in filelist:
        print "Loading: %s" % fname
        inputs = np.loadtxt(fname, dtype= [("f", "S100")])

        path = ""
        start_index = 0
        if os.path.isabs(inputs[0][0]):
            path = inputs[0][0]
            start_index = 1

        for i in xrange(start_index, len(inputs)):
            output.append(os.path.join(path, inputs[i][0]))

    return output
    

def readmosfits(fname, options, extension=None):
    '''Read a fits file written by MOSFIRE from path and return a tuple of 
    (header, data, Target List, Science Slit List (SSL), Mechanical Slit 
    List (MSL), Alignment Slit List (ASL)).
    
    Note, the extension is typically not used, only used if the detector server
    does not append slit extension.
    '''

    if os.path.isabs(fname): path = fname
    else: path = os.path.join(fname_to_path(fname, options), fname)

    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    theBPM = badpixelmask()
    data = np.ma.masked_array(data, theBPM)

    if extension is not None:
        hdulist = pf.open(extension)

    try:
        header = hdulist[0].header
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        raise Exception("Improper MOSFIRE FITS File: %s" % path)

    if np.abs(header["REGTMP1"] - 77) > 0.1:
        print "**************************************"
        print ("The temperature of the detector is %3.3f where it "
                "should be 77.000 deg. Please notify Keck support staff." %
                header["REGTMP1"])

    ssl = ssl[ssl.field("Slit_Number") != ' ']
    msl = msl[msl.field("Slit_Number") != ' ']
    asl = asl[asl.field("Slit_Number") != ' ']

    bs = CSU.Barset()
    bs.set_header(header, ssl=ssl, msl=msl, asl=asl, targs=targs)

    return (header, data, bs)

def readscitbl(path):

    print path

    hdulist = pf.open(path)
    header = hdulist[0].header
    try:
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        print "Improper MOSFIRE FITS File: %s" % path

    return header, targs, ssl, msl, asl


def parse_header_for_bars(header):
    '''Parse {header} and convert to an array of CSU bar positions in mm. If 
    the positon is negative it means the barstat is not OK'''

    poss = []
    posfmt = "B%2.2iPOS"
    statfmt = "B%2.2iSTAT"
    for i in range(1,CSU.numbars+1):
        p = posfmt % i
        s = statfmt % i
        pos = np.float32(header[p])
        if (header[s] != 'OK') and (header[s] != 'SETUP'):
            pos *= -1
        poss.append(pos)

    if len(poss) != CSU.numbars:
        raise CSU.MismatchError("Found %i bars instead of %i" % 
                (lens(poss), CSU.numbars))
        

    return np.array(poss)


def floatcompress(data, ndig=14):
    '''Adapted from Finkbeiner IDL routine floatcompress'''

    t = data.dtype
    if not ((t == 'float32') or (t == 'float64')):
        raise Exception("Only works on floating point numbers")

    wzer = np.where(data == 0)
    data[wzer] = 1.0

    log2 = np.ceil(np.log(np.abs(data)) / np.log(2.0))
    mant = np.round(data/2.0**(log2 - ndig))/2.0**ndig
    out = mant*2.0**log2

    out[wzer] = 0.0
    return out

def imarith(operand1, op, operand2, result):
    from pyraf import iraf
    iraf.images()

    pars = iraf.imarith.getParList()
    iraf.imcombine.unlearn()

    print "%s %s %s -> %s" % (operand1, op, operand2, result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

    iraf.imarith.setParList(pars)

def imcombine(filelist, out, options, bpmask=None, reject="none", nlow=None,
        nhigh=None):

    '''Convenience wrapper around IRAF task imcombine

    Args:
        filelist: The list of files to imcombine
        out: The full path to the output file
        options: Options dictionary
        bpmask: The full path to the bad pixel mask
        reject: none, minmax, sigclip, avsigclip, pclip
        nlow,nhigh: Parameters for minmax rejection, see iraf docs
    
    Returns:
        None

    Side effects:
        Creates the imcombined file at location `out'
    '''

    #TODO: REMOVE Iraf and use python instead. STSCI Python has
    # A builtin routine.
    from pyraf import iraf
    iraf.images()


    filelist = [("%s[0]" % f) for f in filelist]
    pars = iraf.imcombine.getParList()
    iraf.imcombine.unlearn()

    path = "flatcombine.lst"
    f = open(path, "w")
    for file in filelist:
        f.write(file + "\n")
    f.close()

    s = ("%s," * len(filelist))[0:-1]
    s = s % tuple(filelist)

    if reject == 'minmax':
        t = iraf.imcombine("@%s" % path, out, #Stdin=filelist, Stdout=1,
            reject=reject, nlow=nlow, nhigh=nhigh)
    else:
        t = iraf.imcombine(s, out, Stdin=filelist, Stdout=1,
            reject=reject)

    iraf.imcombine.setParList(pars)


    

class TestIOFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_readfits(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()

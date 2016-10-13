'''
MOSFIRE Input/Output Utility Code
Written March 2, 2011 by npk

Provides tools to read fits files and parse their headers.
'''

try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import numpy as np
import unittest
import warnings

import re

import os
import pdb
import shutil

import ccdproc

import MOSFIRE
import CSU
import Options
from MosfireDrpLog import debug, info, warning, error

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
    try:
        edges = np.load(fn)
    except:
        error("Cannot load slit edges file")
        raise Exception("Cannot load slit edges file")
    edges,meta = edges[0:-1], edges[-1]

    if meta['maskname'] != maskname:
        warning("The maskname for the edge file '%s' does not match "
                "that in the edge file '%s'" % (maskname, meta['maskname']))
        warning("Continuing")

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


def load_lambdaslit(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    if False:
        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_solution_{0}.fits".format(fnum))

    fn = "lambda_solution_{0}.fits".format(fnum)

    print fn

    ret = readfits(fn, options)
    if ret[0]['filter'] != band:
        error ("Band name mismatch")
        raise Exception("band name mismatch")

    if ret[0]['maskname'] != maskname:
        warning("The maskname for the edge file '%s' does not match "
                "that in the edge file '%s'" % (maskname, ret[0]['maskname']))
        warning("Continuing")

    
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
            if k in hdu.header: continue

            if k == 'COMMENT': continue
            if k == '': continue

            k = k.rstrip()
            hdu.header[k] = (value,comment)

    warnings.filterwarnings('always')
    if overwrite:
        try: 
            os.remove(fn)
            debug("Removed old file '{0}'".format(fn))
        except: pass

    info("Wrote to '%s'" % (fn))
    warnings.filterwarnings('ignore','Card is too long, comment will be truncated.')
    hdu.writeto(fn)
    warnings.filterwarnings('always')
    if lossy_compress: os.system("gzip --force {0}".format(fn))



def readfits(path, use_bpm=False):
    '''Read a fits file from path and return a tuple of (header, data, 
    Target List, Science Slit List (SSL), Mechanical Slit List (MSL),
    Alignment Slit List (ASL)).'''

    if os.path.exists(path + ".gz"):
        path = path + ".gz"

    if not os.path.exists(path):
        error("The file at path '%s' does not exist." % path)
        raise Exception("The file at path '%s' does not exist." % path)

    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data
    datasec = ""
    try:
        datasec = header["DATASEC"]
        debug("%s contains a DATASEC keyword not compatible with the pipeline" % path)
        debug("The content of the keyword will be erased on the reduced data")
        del header["DATASEC"]
    except:
        pass
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
        error("The file at path '%s' does not exist." % path)
        raise Exception("The file at path '%s' does not exist." % path)

    hdulist = pf.open(path)
    output = []

    for hdu in hdulist:
        output.append(hdu.header)

        if "DRPVER" in hdu.header:
            itsver = hdu.header["DRPVER"]
            if itsver != MOSFIRE.__version__:
                error("The file requested '%s' uses DRP version %f "
                    "but the current DRP version is %f. There might be an "
                    "incompatibility" % (path, itsver, MOSFIRE.__version__))
                raise Exception("The file requested '%s' uses DRP version %f "
                    "but the current DRP version is %f. There might be an "
                    "incompatibility" % (path, itsver, MOSFIRE.__version__))

        else:
            error("The file requested '%s' does not seem to be "
                    "the result of this DRP. This should never be the "
                    " case.")
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
        warning("Could not parse date out of file name: %s" % (fname))
    
    return yr, month, dy

def fname_to_path(fname, options):
    '''Take a filename like m120507_0123, parse date, and return full path'''

    if os.path.isabs(fname): return fname

    yr, month, dy = fname_to_date_tuple(fname)
    path = os.path.join(options["indir"], yr + month + "%2.2i" % dy)
    if not os.path.exists(os.path.join(path, fname)):
        path = os.path.join(options["indir"], yr + month + "%2.2i" % (dy-1))

        if not os.path.exists(path):
            error("Could not find file '%s' in '%s' out of parsed "
                "%s, %s, %s" % (fname,
                options["indir"], yr, month, dy))
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
        debug( "Loading: %s" % fname)
        inputs = np.loadtxt(fname, dtype= [("f", "S100")])
        path = ""
        start_index = 0
        if len(inputs):
            if os.path.isabs(inputs[0][0]):
                path = inputs[0][0]
                start_index = 1

            for i in xrange(start_index, len(inputs)):
                output.append(os.path.join(path, inputs[i][0]))

    return output

def fix_long2pos_headers(filelist):
    '''Fixes old long2pos observations which have a wrong set of keywords'''
    files = list_file_to_strings(filelist)
    # Print the filenames to Standard-out
    info("Fixing long2pos headers for files in "+str(filelist))

    # Iterate through files
    for fname in files:
        if os.path.isabs(fname): path = fname
        else: path = os.path.join(fname_to_path(fname, options), fname)

        hdulist = pf.open(path, mode='update')
        header = hdulist[0].header

        # determine if this file really needs to be updated (for example,
        # prevents a second update of an already updated file
        if 'long2pos' in header['MASKNAME'] and header['FRAMEID']=='object'\
            and (header['PATTERN']=='long2pos' or header['PATTERN']=='Stare'):
            info( "File "+str(fname)+" will be updated")

            # make a copy of the original file
            newname = path+".original"
            info("copying ... "+str(path))
            info("into ...... "+str(newname))
            shutil.copyfile(path,newname)
            if not os.path.exists(newname):
                errstr = "Error in generating original file:  '%s' does not exist"\
                         "(could not be created)." % newname
                error(errstr)
                raise Exception(errstr)

            #updating header
            # assign FRAMEID to narrow slits
            if header['YOFFSET']==21 or header['YOFFSET']==-7:
                header['FRAMEID']="B"
            if header['YOFFSET']==-21 or header['YOFFSET']==7:
                header['FRAMEID']="A"

            # assign FRAMEID to wide slits
            if header['YOFFSET']==14 or header['YOFFSET']==-14:
                header['FRAMEID']="A"
                
            #reverse sign of offsets for narrow slits
            if header['YOFFSET']==-21:
                header['YOFFSET']=7
            if header['YOFFSET']==21:
                header['YOFFSET']=-7

            #transform Xoffset from pixels to arcseconds
            header['XOFFSET'] = header['XOFFSET']*0.18
        else:
            info("File "+str(fname)+" does not need to be updated")
        hdulist.flush()
        hdulist.close()


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
        datasec = ""
        try:
            datasec = header["DATASEC"]
            debug("%s contains a DATASEC keyword not compatible with the pipeline" % path)
            debug("The content of the keyword will be erased on the reduced data")
            del header["DATASEC"]
        except:
            pass
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        error("Improper MOSFIRE FITS File: %s" % path)
        raise Exception("Improper MOSFIRE FITS File: %s" % path)

#     if np.abs(header["REGTMP1"] - 77) > 0.1:
#         warning("**************************************")
#         warning("The temperature of the detector is %3.3f where it "
#                 "should be 77.000 deg. Please notify Keck support staff." %
#                 header["REGTMP1"])

    ssl = ssl[ssl.field("Slit_Number") != ' ']
    msl = msl[msl.field("Slit_Number") != ' ']
    asl = asl[asl.field("Slit_Number") != ' ']
        

    # ELIMINATE POSITION B of the long2pos slit
    ssl = ssl[ssl.field("Target_Name") != 'posB']
    msl = msl[msl.field("Target_in_Slit") != 'posB']
    asl = asl[asl.field("Target_in_Slit") != 'posBalign']
    targs = targs[targs.field("Target_Name") !='posB']
    targs = targs[targs.field("Target_Name") != "posBalign"]

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
        warning("Improper MOSFIRE FITS File: %s" % path)

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
        error("Found %i bars instead of %i" % 
                (lens(poss), CSU.numbars))
        raise CSU.MismatchError("Found %i bars instead of %i" % 
                (lens(poss), CSU.numbars))
        

    return np.array(poss)


def floatcompress(data, ndig=14):
    '''Adapted from Finkbeiner IDL routine floatcompress'''

    t = data.dtype
    if not ((t == 'float32') or (t == 'float64')):
         error("Only works on floating point numbers")
         raise Exception("Only works on floating point numbers")

    wzer = np.where(data == 0)
    data[wzer] = 1.0

    log2 = np.ceil(np.log(np.abs(data)) / np.log(2.0))
    mant = np.round(data/2.0**(log2 - ndig))/2.0**ndig
    out = mant*2.0**log2

    out[wzer] = 0.0
    return out

def imarith(operand1, op, operand2, result):
    info( "%s %s %s -> %s" % (operand1, op, operand2, result))
    assert type(operand1) == str
    assert type(operand2) == str
    assert os.path.exists(operand1)
    assert os.path.exists(operand2)
    assert op in ['+', '-']

    ## Strip off the [0] part of the operand as we are assuming that we are
    ## operating on the 0th FITS HDU.
    if re.match('(\w+\.fits)\[0\]', operand1):
        operand1 = operand1[:-3]
    if re.match('(\w+\.fits)\[0\]', operand2):
        operand2 = operand2[:-3]

    import operator
    operation = { "+": operator.add, "-": operator.sub,\
                  "*": operator.mul, "/": operator.truediv}

    hdulist1 = pf.open(operand1, 'readonly')
    hdulist2 = pf.open(operand2, 'readonly')
    data = operation[op](hdulist1[0].data, hdulist2[0].data)
    header = hdulist1[0].header
    header['history'] = 'imarith {} {} {}'.format(operand1, op, operand2)
    header['history'] = 'Header values copied from {}'.format(operand1)
    if 'EXPTIME' in hdulist1[0].header and\
       'EXPTIME' in hdulist2[0].header and\
       operation in ['+', '-']:
        exptime = operation[op](float(hdulist1[0].header['EXPTIME']),\
                                float(hdulist2[0].header['EXPTIME']))
        header['history'] = 'Other than exposure time which was edited'
        header['EXPTIME'] = exptime

    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(result)


def imcombine(filelist, out, options, method="average", reject="none",\
              lsigma=3, hsigma=3, mclip=False,\
              nlow=None, nhigh=None):
    '''Combines images in input list with optional rejection algorithms.

    Args:
        filelist: The list of files to imcombine
        out: The full path to the output file
        method: either "average" or "median" combine
        options: Options dictionary
        bpmask: The full path to the bad pixel mask
        reject: none, minmax, sigclip
        nlow,nhigh: Parameters for minmax rejection, see iraf docs
        mclip: use median as the function to calculate the baseline values for
               sigclip rejection?
        lsigma, hsigma: low and high sigma rejection thresholds.
    
    Returns:
        None

    Side effects:
        Creates the imcombined file at location `out'
    '''
    assert method in ['average', 'median']
    if os.path.exists(out):
        os.remove(out)

    if reject == 'none':
        info('Combining files using ccdproc.combine task')
        info('  reject=none')
        for file in filelist:
            debug('  Combining: {}'.format(file))
        ccdproc.combine(filelist, out, method=method,\
                        minmax_clip=False,\
                        iraf_minmax_clip=True,\
                        sigma_clip=False,\
                        unit="adu")
        info('  Done.')
    elif reject == 'minmax':
        ## The IRAF imcombine minmax rejection behavior is different than the
        ## ccdproc minmax rejection behavior.  We are using the IRAF like
        ## behavior here.  To support this a pull request for the ccdproc
        ## package has been made:
        ##    https://github.com/astropy/ccdproc/pull/358
        ##
        ## Note that the ccdproc behavior still differs slightly from the
        ## nominal IRAF behavior in that the rejection does not consider whether
        ## any of the rejected pixels have been rejected for other reasons, so
        ## if nhigh=1 and that pixel was masked for some other reason, the
        ## new ccdproc algorithm, will not mask the next highest pixel, it will
        ## still just mask the highest pixel even if it is already masked.
        ##
        ## From IRAF (help imcombine):
        ##  nlow = 1,  nhigh = 1 (minmax)
        ##      The number of  low  and  high  pixels  to  be  rejected  by  the
        ##      "minmax"  algorithm.   These  numbers are converted to fractions
        ##      of the total number of input images so  that  if  no  rejections
        ##      have  taken  place  the  specified number of pixels are rejected
        ##      while if pixels have been rejected by masking, thresholding,  or
        ##      non-overlap,   then   the  fraction  of  the  remaining  pixels,
        ##      truncated to an integer, is used.
        ##

        ## Check that minmax rejection is possible given the number of images
        if nlow is None:
            nlow = 0
        if nhigh is None:
            nhigh = 0
        if nlow + nhigh >= len(filelist):
            warning('nlow + nhigh >= number of input images.  Combining without rejection')
            nlow = 0
            nhigh = 0
        
        if ccdproc.version.major >= 1 and ccdproc.version.minor >= 1\
           and ccdproc.version.release:
            info('Combining files using ccdproc.combine task')
            info('  reject=clip_extrema')
            info('  nlow={}'.format(nlow))
            info('  nhigh={}'.format(nhigh))
            for file in filelist:
                info('  {}'.format(file))
            ccdproc.combine(filelist, out, method=method,\
                            minmax_clip=False,\
                            clip_extrema=True,\
                            nlow=nlow, nhigh=nhigh,\
                            sigma_clip=False,\
                            unit="adu")
            info('  Done.')
        else:
            ## If ccdproc does not have new rejection algorithm in:
            ## https://github.com/astropy/ccdproc/pull/358
            ## Manually perform rejection using ccdproc.combiner.Combiner object
            info('Combining files using local clip_extrema rejection algorithm')
            info('and the ccdproc.combiner.Combiner object.')
            info('  reject=clip_extrema')
            info('  nlow={}'.format(nlow))
            info('  nhigh={}'.format(nhigh))
            for file in filelist:
                info('  {}'.format(file))
            ccdlist = []
            for file in filelist:
                ccdlist.append(ccdproc.CCDData.read(file, unit='adu', hdu=0))
            c = ccdproc.combiner.Combiner(ccdlist)
            nimages, nx, ny = c.data_arr.mask.shape
            argsorted = np.argsort(c.data_arr.data, axis=0)
            mg = np.mgrid[0:nx,0:ny]
            for i in range(-1*nhigh, nlow):
                where = (argsorted[i,:,:].ravel(),
                         mg[0].ravel(),
                         mg[1].ravel())
                c.data_arr.mask[where] = True
            if method == 'average':
                result = c.average_combine()
            elif method == 'median':
                result = c.median_combine()
            for key in ccdlist[0].header.keys():
                header_entry = ccdlist[0].header[key]
                if key != 'COMMENT':
                    result.header[key] = (header_entry,
                                          ccdlist[0].header.comments[key])
            hdul = result.to_hdu()
#             print(hdul)
#             for hdu in hdul:
#                 print(type(hdu.data))
            hdul[0].writeto(out)
#             result.write(out)
            info('  Done.')
    elif reject == 'sigclip':
        info('Combining files using ccdproc.combine task')
        info('  reject=sigclip')
        info('  mclip={}'.format(mclip))
        info('  lsigma={}'.format(lsigma))
        info('  hsigma={}'.format(hsigma))
        baseline_func = {False: np.mean, True: np.median}
        ccdproc.combine(filelist, out, method=method,\
                        minmax_clip=False,\
                        clip_extrema=False,\
                        sigma_clip=True,\
                        sigma_clip_low_thresh=lsigma,\
                        sigma_clip_high_thresh=hsigma,\
                        sigma_clip_func=baseline_func[mclip],\
                        sigma_clip_dev_func=np.std,\
                        )
        info('  Done.')
    else:
        raise NotImplementedError('{} rejection unrecognized by MOSFIRE DRP'.format(reject))




class TestIOFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_readfits(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()

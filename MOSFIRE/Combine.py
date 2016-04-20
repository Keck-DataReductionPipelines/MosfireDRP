import os
import re

from astropy.io import fits as pf
import ccdproc

from MosfireDrpLog import debug, info, warning, error

theBPM = None # the Bad pixel mask


def imarith(operand1, op, operand2, result):
    imarith_noiraf(operand1, op, operand2, result)

def imarith_iraf(operand1, op, operand2, result):
    from pyraf import iraf
    iraf.images()

    pars = iraf.imarith.getParList()
    iraf.imcombine.unlearn()

    info( "%s %s %s -> %s" % (operand1, op, operand2, result))
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)

    iraf.imarith.setParList(pars)

def imarith_noiraf(operand1, op, operand2, result):
    '''Python replacement for the IRAF imarith function.
    
    Currently this only supports + and - operators.
    '''
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
    operation = { "+": operator.add, "-": operator.sub }

    hdulist1 = fits.open(operand1, 'readonly')
    hdulist2 = fits.open(operand2, 'readonly')
    data = operation[op](hdulist1[0].data, hdulist2[0].data)
    header = hdulist1[0].header
    header['history'] = 'imarith {} {} {}'.format(operand1, op, operand2)
    header['history'] = 'Header values copied from {}'.format(operand1)
    if 'EXPTIME' in hdulist1[0].header and 'EXPTIME' in hdulist2[0].header:
        exptime = operation[op](float(hdulist1[0].header['EXPTIME']),\
                                float(hdulist2[0].header['EXPTIME']))
        header['EXPTIME'] = exptime
    header['history'] = 'Other than exposure time which was edited'

    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(result)

def imcombine(filelist, out, options, bpmask=None, reject="none", nlow=None,
        nhigh=None):
    imcombine_noiraf(filelist, out, options, bpmask=bpmask, reject="none", nlow=nlow, nhigh=nlow)


def imcombine_iraf(filelist, out, options, bpmask=None, reject="none", nlow=None,
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

    f = open("flatcombinelog.txt", "w")
    if reject == 'minmax':
        t = iraf.imcombine("@%s" % path, out, Stdout=f,
            reject=reject, nlow=nlow, nhigh=nhigh)
    else:
        t = iraf.imcombine(s, out, Stdin=filelist, Stdout=f,
            reject=reject)
    f.close()
    f=open("flatcombinelog.txt")
    for line in f:
        info(line.rstrip("\n"))
    f.close()

    iraf.imcombine.setParList(pars)


def imcombine_noiraf(filelist, out, options, bpmask=None, reject="none", nlow=None,
        nhigh=None):
    '''Replacement for iraf.imcombine which uses the ccdproc.combine method.

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
    if reject == 'none':
        info('Combining files using ccdproc.combine task')
        info('  reject=none')
        for file in filelist:
            info('  {}'.format(file))
        ccdproc.combine(filelist, out, method='average',\
                        minmax_clip=False,\
                        sigma_clip=False,\
                        unit="adu")
        info('  Done.')
    elif reject == 'minmax':
        ## The IRAF imcombine parameter for minmax rejection specifies the
        ## number of pixels to clip, while the analogous parameters for
        ## ccdproc.combine specify the pixel values above and below which to
        ## clip, so a conversion will have to be made.
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
        ## From ccdproc.combine doc string:
        ##  minmax_clip : Boolean (default False)
        ##      Set to True if you want to mask all pixels that are below
        ##      minmax_clip_min or above minmax_clip_max before combining.
        ##  
        ##      Parameters below are valid only when minmax_clip is set to True.
        ##  
        ##      minmax_clip_min: None, float
        ##           All pixels with values below minmax_clip_min will be masked.
        ##      minmax_clip_max: None or float
        ##           All pixels with values above minmax_clip_max will be masked.
        raise NotImplementedError('minmax rejection is not yet implemented')
#         clip_min = 
#         clip_max = 
#         ccdproc.combine(filelist, out, method='average',\
#                         minmax_clip=True,\
#                         minmax_clip_min=clip_min, minmax_clip_max=clip_max,\
#                         sigma_clip=False)
    elif reject == 'sigclip':
        raise NotImplementedError('sigclip rejection is not yet implemented')
#         ccdproc.combine(filelist, out, method='average',\
#                         minmax_clip=False,\
#                         sigma_clip=True,\
#                         sigma_clip_low_thresh=
#                         sigma_clip_high_thresh=
#                         sigma_clip_func=np.mean,\
#                         sigma_clip_dev_func=np.std,\
#                         )
    elif reject == 'avsigclip':
        raise NotImplementedError('avsigclip rejection is not yet implemented')
    elif reject == 'pclip':
        raise NotImplementedError('pclip rejection is not yet implemented')
    else:
        raise NotImplementedError('{} rejection unrecognized by MOSFIRE DRP'.format(reject))



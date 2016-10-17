#!/usr/env/python

from __future__ import division, print_function

## Import General Tools
import os
import sys
import textwrap

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import table
from astropy import wcs
from astropy.modeling import models, fitting
import scipy.signal as signal

from matplotlib import pyplot as pl
pl.rcParams['keymap.fullscreen'] = ''
pl.rcParams['keymap.grid'] = ''
pl.rcParams['keymap.save'] = ''

from MOSFIRE import Detector
from MOSFIRE.MosfireDrpLog import info, debug, warning, error

RN = Detector.RN * u.electron
Gain = Detector.gain * u.electron/u.adu
py3 = sys.version_info[0] > 2 #creates boolean value for test that Python > 2

##-------------------------------------------------------------------------
## Find Stellar Traces
##-------------------------------------------------------------------------
def print_instructions():
    '''Print instructions for interactive fit to screen.
    '''
    text = [
    ('The figure shows the raw spectrum collapsed in the wavelength'
    ' direction as black points.  If there are one or more object traces '
    'with reasonable '
    'signal to noise, you should see one or more signals (e.g. a positive '
    'gaussian profile).'
    ),

    ('You can use this tool to add, modify, and delete apertures.  '
    'The apertures are indicated by a yellow shaded region and their center '
    'position and half width in pixels is annotated near the top of each '
    'shaded region.  There may be automatically generated regions already '
    'shown if that option was selected when the software was run.'
    ),
    
    ('The apertures define the pixels which will be used as input to the '
    'optimal spectral extraction (Horne 1986) algorithm.  Having wide a '
    'wide aperture should not add additional noise as that will be '
    'optimized during the spectral extraction step.  The apertures are shown '
    'here in order for the user to verify 1) that there is no overlap between '
    'adjacent objects, 2) that the apertures are wide enough to reasonably '
    'encompass all flux from the object, and 3) that all objects have '
    'properly defined apertures.'
    ),

    ('To delete an existing aperture: place the mouse near the center of the '
    'aperture and press the "d" key.'
    ),
    
    ('To add an aperture by fitting a gaussian to the profile: place the mouse '
    'near the peak of the profile and press the "g" key.  The half width of '
    'the aperture will be set at 5 times the sigma of the fitted gaussian.  '
    'If a gaussian fit has been used to definte an aperture, the fit will be '
    'shown as a blue line.'
    ),

    ('To add an aperture manually: place the mouse in the X position where the '
    'new aperture should be centered and press the "a" key.  Then type the half'
    ' width (in pixels) for that aperture in response to the query in the '
    'terminal.'
    ),

    ('To modify the half width of an existing aperture: place the mouse near '
    'the center of the aperture and press the "w" key.  Then type the half'
    ' width (in pixels) for that aperture in response to the query in the '
    'terminal.'
    ),

    ('To modify the center position of an existing aperture: place the mouse '
    'near the center of the aperture and press the "p" key.  Then type the '
    'position (in pixels) for that aperture in response to the query in the '
    'terminal.'
    ),

    ('When you are done adding or removing apertures, close the interactive '
    'plot window by clicking the close button in the upper right corner '
    '(or by whatever method is typical for your OS or windowing system) '
    'or press the "q" or "n" keys (for "quit" or "next" respectively).'
    ),
    ]

    print('#'*80)
    print('                     Aperture Definition Tool Instructions')
    for paragraph in text:
        print()
        print(textwrap.fill(textwrap.dedent(paragraph).strip('\n'), width=80))
    print('#'*80)
    print()


class ApertureEditor(object):
    def __init__(self, hdu, title=None):
        self.header = hdu.header
        self.data = np.ma.MaskedArray(data=hdu.data, mask=np.isnan(hdu.data))
        self.ydata = np.mean(self.data, axis=1)
        self.xdata = range(0,len(self.ydata), 1)
        self.fig = pl.figure()
        self.ax = self.fig.gca()
        self.title = title
        if title:
            info('Editing apertures for {}'.format(title))
        else:
            info('Editing apertures')
        self.apertures = table.Table(names=('id', 'position', 'width',\
                                              'amplitude', 'sigma'),\
                                       dtype=('i4', 'f4', 'f4', 'f4', 'f4'))


    def add_aperture(self, pos, width):
        info('  Adding aperture at position {} of width {}'.format(pos, width))
        id = len(self.apertures)
        data = {'id': id,
                'position': pos,
                'width': width,
                'amplitude': None,
                'sigma': None,
               }
        self.apertures.add_row(data)
        self.plot_apertures()


    def delete_aperture(self, id):
        pos = self.apertures[id]['position']
        info('  Removing aperture at position {:.0f}'.format(pos))
        self.apertures.remove_row(id)
        self.plot_apertures()


    def set_position(self, id=None, pos=None):
        assert id is not None
        info('  Moving aperture {} from {:.0f} to {:.0f}'.format(id,
             float(self.apertures[id]['position']), float(pos)))
        if pos:
            self.apertures[id]['position'] = pos
        else:
            if py3:
                response = input("Please enter new position in pixels: ")
            else:
                response = raw_input("Please enter new position in pixels: ")
            pos = int(response)
            self.apertures[id]['position'] = pos
        self.plot_apertures()


    def set_width(self, id=None, width=None):
        assert id is not None
        info('  Changing width for aperture {} at position {:.0f}'.format(id,
             self.apertures[id]['position']))
        if width:
            self.apertures[id]['width'] = width
        else:
            if py3:
                response = input("Please enter new half width in pixels: ")
            else:
                response = raw_input("Please enter new half width in pixels: ")
            width = int(response)
            self.apertures[id]['width'] = width
        self.plot_apertures()


    def guess(self):
        ## Try to guess at aperture positions if no information given
        ## Start with maximum pixel value
        valid_profile = list(self.ydata[~self.ydata.mask])
        maxval = max(valid_profile)
        maxind = valid_profile.index(maxval)
        info('  Guessing at aperture near position {}'.format(maxind))
        return (maxind, maxval)
#         self.fit_trace(maxind, maxval)


    def fit_trace(self, pos, amp):
        info('  Adding fitted aperture near position {:.0f}'.format(pos))
        id = len(self.apertures)
        g0 = models.Gaussian1D(mean=pos, amplitude=amp,
                              bounds={'amplitude': [0, float('+Inf')]})
        fitter = fitting.LevMarLSQFitter()
        g = fitter(g0, self.xdata, self.ydata)
        # Set maximum width of aperture to the YOFFSET parameter
        try:
            maxap = np.floor(float(self.header['YOFFSET'])/0.1799)
        except:
            maxap = self.data.shape[0]
        # Set minimum width of aperture to an estimate of 3x the seeing
        seeing = 0.5 # arcsec
        minap = np.ceil(3*seeing/0.1799)
        width = np.ceil(5.*g.param_sets[2][0])
        width = max(min(maxap, width), minap)
        data = {'id': id,
                'position': g.param_sets[1][0],
                'width': width,
                'amplitude': g.param_sets[0][0],
                'sigma': g.param_sets[2][0],
               }
        self.apertures.add_row(data)
        self.plot_apertures()


    def plot_data(self):
        '''Plot the raw data without apertures.
        '''
        debug('  Plotting profile data')
        pl.plot(self.xdata, self.ydata, 'k-',
                label='Spatial Profile', drawstyle='steps-mid')
        pl.xlim(min(self.xdata), max(self.xdata))
        yspan = self.ydata.max() - self.ydata.min()
        pl.ylim(self.ydata.min()-0.02*yspan, self.ydata.max()+0.18*yspan)
        pl.xlabel('Pixel Position')
        pl.ylabel('Flux (e-/sec)')
        if self.title is None:
            pl.title('Spatial Profile')
        else:
            pl.title(self.title)


    def plot_apertures(self):
        pl.cla()
        self.plot_data()
        yspan = self.ydata.max() - self.ydata.min()
        if len(self.apertures) == 0:
            pl.draw()
        for ap in self.apertures:
            debug('  Plotting aperture at position {}'.format(ap['position']))
            if ap['sigma'] is not None and ap['amplitude'] is not None:
                g = models.Gaussian1D(mean=ap['position'],
                                      amplitude=ap['amplitude'],
                                      stddev=ap['sigma'])
                fit = [g(x) for x in self.xdata]
                pl.plot(self.xdata, fit, 'b-', label='Fit', alpha=0.5)
            shadeymin = np.floor(self.ydata.min())
            shadeymax = np.ceil(self.ydata.max())
            pl.axvspan(ap['position']-ap['width'],
                            ap['position']+ap['width'],
                            ymin=shadeymin,
                            ymax=shadeymax,
                            facecolor='y', alpha=0.3,
                            )
            pl.text(ap['position']-ap['width']+1,
                         self.ydata.max() + 0.05*yspan,
                         'position={:.0f}\nwidth={:.0f}'.format(ap['position'],
                                                                ap['width']),
                         )
            pl.draw()
        pl.show()

    def connect(self):
        '''Connect keypresses to matplotlib for interactivity.
        '''
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event',
                                                   self.keypress)


    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_key)


    def determine_id(self, event):
        x = event.xdata
        closest = None
        for i,ap in enumerate(self.apertures):
            id = ap['id']
            d = abs(ap['position'] - x)
            if closest is None:
                closest = (id, d)
            elif d < closest[1]:
                closest = (id, d)
        debug('  ID of aerture nearest to keypress is {}'.format(closest[0]))
        return closest[0]


    def keypress(self, event):
        '''Based on which key is presses on a key press event, call the
        appropriate method.
        '''
        if event.key == 'a':
            if py3:
                response = input("Please enter new half width in pixels: ")
            else:
                response = raw_input("Please enter new half width in pixels: ")
            width = int(response)
            self.add_aperture(event.xdata, width)
            print('Adding aperture at position {}, width {}'.format(event.xdata,
                  width))
        elif event.key == 'w':
            id = self.determine_id(event)
            self.set_width(id=id)
        elif event.key == 'p':
            id = self.determine_id(event)
            self.set_position(id=id)
        elif event.key == 'g':
            self.fit_trace(event.xdata, event.ydata/abs(event.ydata))
        elif event.key == 'd':
            id = self.determine_id(event)
            self.delete_aperture(id)
        elif event.key == 'n':
            self.quit(event)
        elif event.key == 'q':
            self.quit(event)
        elif event.key == 's':
            print(self.apertures)
        elif event.key == 'r':
            info('  Plotting apertures')
            self.plot_apertures()


    def savefig(self, plotfile):
        '''Save the figure to a png file.
        '''
        self.fig.savefig(plotfile, bbox_inches='tight')


    def quit(self, event):
        info('  Done with aperture edits for this slit')
        self.disconnect()
        pl.close(self.fig)


##-------------------------------------------------------------------------
## Find Stellar Traces
##-------------------------------------------------------------------------
def find_apertures(hdu, guesses=[], title=None, interactive=True):
    '''Finds targets in spectra by simply collapsing the 2D spectra in the
    wavelength direction and fitting Gaussian profiles to the positional profile
    '''
    pl.ioff()
    ap = ApertureEditor(hdu, title=title)
    if interactive:
        ap.connect()

    if (guesses == []) or (guesses is None):
        # Guess at object position where specified in header by CRVAL2
        pos = int(-hdu.header['CRVAL2'])
        amp = ap.ydata[int(pos)]
        # Estimate signal to noise of profile at that spot
        # First, clip top and bottom 10 percent of pixels to roughly remove
        # contribution by a single bright source.
        pct = 10
        pctile = (np.percentile(ap.ydata, pct), np.percentile(ap.ydata, 100-pct))
        w = np.where((ap.ydata > pctile[0]) & (ap.ydata < pctile[1]))
        # Use the unclipped pixels to estimate signal to noise.
        std = np.std(ap.ydata[w])
        snr = amp/std
        # If SNR is strong, fit a gaussian, if not, just blindly add an aperture
        if snr > 5:
            ap.fit_trace(pos, amp)
        else:
            width = 10
            ap.add_aperture(pos, width)
    else:
        for guess in guesses:
            if len(guess) == 1:
                ap.fit_trace(guess, 1)
            else:
                ap.add_aperture(guess[0], guess[1])
    return ap.apertures


##-------------------------------------------------------------------------
## Standard Spectral Extraction
##-------------------------------------------------------------------------
def standard_extraction(data, variance):
    spect_1D = np.sum(data, axis=0)
    variance_1D = np.sum(variance, axis=0)
    return spect_1D, variance_1D


##-------------------------------------------------------------------------
## Iterate Spatial Profile
##-------------------------------------------------------------------------
def iterate_spatial_profile(P, DmS, V, f,\
                            smoothing=5, order=3, minpixels=50,\
                            sigma=4, nclip=2, verbose=True):
    poly0 = models.Polynomial1D(degree=order)
    fit_poly = fitting.LinearLSQFitter()    

    Pnew = np.zeros(P.shape)
    for i,row in enumerate(P):
        weights = f**2/V[i]
        weights.mask = weights.mask | np.isnan(weights)
        srow = np.ma.MaskedArray(data=signal.medfilt(row, smoothing),\
                     mask=(np.isnan(signal.medfilt(row, smoothing)) | weights.mask))
        xcoord = np.ma.MaskedArray(data=np.arange(0,len(row),1),\
                                   mask=srow.mask)

        for iter in range(nclip+1):
            nrej_before = np.sum(srow.mask)
            fitted_poly = fit_poly(poly0,\
                                   xcoord[~xcoord.mask], srow[~srow.mask],\
                                   weights=weights[~weights.mask])
            fit = np.array([fitted_poly(x) for x in xcoord])
            resid = (DmS[i]-f*srow)**2 / V[i]
            newmask = (resid > sigma)
            
            weights.mask = weights.mask | newmask
            srow.mask = srow.mask | newmask
            xcoord.mask = xcoord.mask | newmask

            nrej_after = np.sum(srow.mask)
            if nrej_after > nrej_before:
                if verbose:\
                    info('Row {:3d}: Rejected {:d} pixels on clipping '+\
                          'iteration {:d}'.format(\
                           i, nrej_after-nrej_before, iter))
        
        ## Reject row if too few pixels are availabel for the fit
        if (srow.shape[0] - nrej_after) < minpixels:
            if verbose:
                warning('Row {:3d}: WARNING! Only {:d} pixels remain after '+\
                        'clipping and masking'.format(\
                      i, srow.shape[0] - nrej_after))
            fit = np.zeros(fit.shape)
        ## Set negative values to zero
        if np.sum((fit<0)) > 0 and verbose:
            info('Row {:3d}: Reset {:d} negative pixels in fit to 0'.format(\
                  i, np.sum((fit<0))))
        fit[(fit < 0)] = 0
        Pnew[i] = fit
    
    return Pnew


##-------------------------------------------------------------------------
## Optimal Spectral Extraction
##-------------------------------------------------------------------------
def optimal_extraction(image, variance_image, aperture_table,
                       fitsfileout=None,
                       plotfileout=None,
                       plot=None):
    '''Given a 2D spectrum image, a 2D variance image, and a table of apertures
    (e.g. as output by find_apertures() above), this function will optimally
    extract a 1D spectrum for each entry in the table of apertures.
    
    
    '''
    if type(image) == fits.HDUList:
        hdu = image[0]
    elif type(image) == fits.PrimaryHDU:
        hdu = image
    else:
        error('Input to standard_extraction should be an HDUList or an HDU')
        raise TypeError

    if type(variance_image) == fits.HDUList:
        vhdu = variance_image[0]
    elif type(image) == fits.PrimaryHDU:
        vhdu = variance_image
    else:
        error('Input to standard_extraction should be an HDUList or an HDU')
        raise TypeError

    spectra2D = hdu.data
    variance2D = vhdu.data
    header = hdu.header
    w = wcs.WCS(hdu.header)
    
    ## State assumptions
    assert header['DISPAXIS'] == 1
    assert w.to_header()['CTYPE1'] == 'AWAV'
    assert header['CD1_2'] == 0
    assert header['CD2_1'] == 0
    flattened_wcs = w.dropaxis(1)
    assert flattened_wcs.to_header()['CTYPE1'] == 'AWAV'

    ## Replace old WCS in header with the collapsed WCS
    for key in w.to_header().keys():
        if key in header.keys():
            header.remove(key)
    for key in flattened_wcs.to_header().keys():
        header[key] = (flattened_wcs.to_header()[key],\
                       flattened_wcs.to_header().comments[key])

    spectra = []
    variances = []
    for i,row in enumerate(aperture_table):
        pos = row['position']
        width = int(row['width'])
        info('Extracting aperture {:d} at position {:.0f}'.format(i, pos))

        ymin = max([int(np.floor(pos-width)), 0])
        ymax = min([int(np.ceil(pos+width)), spectra2D.shape[0]])

        DmS = np.ma.MaskedArray(data=spectra2D[ymin:ymax,:],\
                                mask=np.isnan(spectra2D[ymin:ymax,:]))
        V = np.ma.MaskedArray(data=variance2D[ymin:ymax,:],\
                              mask=np.isnan(spectra2D[ymin:ymax,:]))
        info('  Performing standard extraction')
        f_std, V_std = standard_extraction(DmS, V)
        info('  Forming initial spatial profile')
        P_init_data = np.array([row/f_std for row in DmS])
        P_init = np.ma.MaskedArray(data=P_init_data,\
                                   mask=np.isnan(P_init_data))
        info('  Fitting spatial profile')
        P = iterate_spatial_profile(P_init, DmS, V, f_std, verbose=False)
        info('  Calculating optimally extracted spectrum')
        f_new_denom = np.ma.MaskedArray(data=np.sum(P**2/V, axis=0),\
                                        mask=(np.sum(P**2/V, axis=0)==0))
        f_opt = np.sum(P*DmS/V, axis=0)/f_new_denom
        var_fopt = np.sum(P, axis=0)/f_new_denom
        sig_fopt = np.sqrt(var_fopt)
        typical_sigma = np.mean(sig_fopt[~np.isnan(sig_fopt)])
        spectra.append(f_opt)
        variances.append(var_fopt)
        info('  Typical level = {:.1f}'.format(np.mean(f_opt)))
        info('  Typical sigma = {:.1f}'.format(typical_sigma))

    mask = np.isnan(np.array(spectra)) | np.isnan(np.array(variances))
    spectra = np.ma.MaskedArray(data=np.array(spectra), mask=mask)
    variances = np.ma.MaskedArray(data=np.array(variances), mask=mask)

    w = wcs.WCS(header).dropaxis(1)
    wh = w.to_header()
    for key in wh.keys():
        header[key] = wh[key]

    for i,row in enumerate(aperture_table):
        hdulist = fits.HDUList([])
        if plotfileout:
            fig = pl.figure(figsize=(16,6))
            wunit = getattr(u, w.to_header()['CUNIT1'])

        sp = spectra[i]
        hdulist.append(fits.PrimaryHDU(data=sp.filled(0), header=header))
        hdulist[0].header['APPOS'] = row['position']
        if plotfileout:
            sigma = np.sqrt(variances[i])
            pix = np.arange(0,sp.shape[0],1)
            wavelengths = w.wcs_pix2world(pix,1)[0]*wunit.to(u.micron)*u.micron
            fillmin = sp-sigma
            fillmax = sp+sigma
            mask = np.isnan(fillmin) | np.isnan(fillmax) | np.isnan(sp)
            pl.fill_between(wavelengths[~mask], fillmin[~mask], fillmax[~mask],\
                             label='uncertainty',\
                             facecolor='black', alpha=0.2,\
                             linewidth=0,\
                             interpolate=True)
            pl.plot(wavelengths, sp, 'k-',
                    label='Spectrum for Aperture {} at {:.0f}'.format(i,
                          row['position']))
            pl.xlabel('Wavelength (microns)')
            pl.ylabel('Flux (e-/sec)')
            pl.xlim(wavelengths.value.min(),wavelengths.value.max())
            pl.ylim(0,1.05*sp.max())
            pl.legend(loc='best')

            bn, ext = os.path.splitext(plotfileout)
            plotfilename = '{}_{:02d}{}'.format(bn, i, ext)
            pl.savefig(plotfilename, bbox_inches='tight')
            pl.close(fig)

        var = variances[i]
        hdulist.append(fits.ImageHDU(data=var.filled(0), header=header))
        hdulist[1].header['APPOS'] = row['position']
        hdulist[1].header['COMMENT'] = 'VARIANCE DATA'
        if fitsfileout:
            bn, ext = os.path.splitext(fitsfileout)
            fitsfilename = '{}_{:02d}{}'.format(bn, i, ext)
            hdulist.writeto(fitsfilename, clobber=True)


##-------------------------------------------------------------------------
## Extract Spectra Function
##-------------------------------------------------------------------------
def extract_spectra(maskname, band, interactive=True, target='default'):

    if target == 'default':
        ## Get objectnames from slit edges
        edges = np.load('slit-edges_{}.npy'.format(band))
        objectnames = [edge['Target_Name'] for edge in edges[:-1]]
        eps_files = ['{}_{}_{}_eps.fits'.format(maskname, band, objectname)
                     for objectname in objectnames]
        sig_files = ['{}_{}_{}_sig.fits'.format(maskname, band, objectname)
                     for objectname in objectnames]
        spectrum_plot_files = ['{}_{}_{}.png'.format(maskname, band, objectname)
                               for objectname in objectnames]
        fits_files = ['{}_{}_{}_1D.fits'.format(maskname, band, objectname)
                      for objectname in objectnames]
    else:
        objectnames = [target]
        eps_files = ['{}_{}_eps.fits'.format(target, band)]
        sig_files = ['{}_{}_sig.fits'.format(target, band)]
        spectrum_plot_files = ['{}_{}.png'.format(target, band)]
        fits_files = ['{}_{}_1D.fits'.format(target, band)]

    aperture_tables = {}
    if interactive:
        print_instructions()
    # First, iterate through all slits and define the apertures to extract
    for i,eps_file in enumerate(eps_files):
        objectname = objectnames[i]
        eps = fits.open(eps_file, 'readonly')[0]
        aperture_tables[objectname] = find_apertures(eps, title=objectname,
                                                     interactive=interactive,
                                                     )
    # Second, iterate through all slits again and perform spectral extraction
    # using the apertures defined above
    for i,eps_file in enumerate(eps_files):
        objectname = objectnames[i]
        sig_file = sig_files[i]
        spectrum_plot_file = spectrum_plot_files[i]
        fits_file = fits_files[i]
        if len(aperture_tables[objectname]) == 0:
            info('No apertures defined for {}. Skipping extraction.'.format(
                 objectname))
        else:
            info('Extracting spectra for {}'.format(objectname))
            eps = fits.open(eps_file, 'readonly')[0]
            sig = fits.open(sig_file, 'readonly')[0]
            try:
                optimal_extraction(eps, sig, aperture_tables[objectname],
                                   fitsfileout=fits_file,
                                   plotfileout=spectrum_plot_file,
                                   )
            except Exception as e:
                warning('Failed to extract spectra for {}'.format(objectname))
                warning(e)


if __name__ == "__main__":
    cwd = os.path.abspath('.')
    upone, band = os.path.split(cwd)
    uptwo, date = os.path.split(upone)
    upthree, maskname = os.path.split(uptwo)
    extract_spectra(maskname, band, interactive=True)

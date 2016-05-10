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

import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import connect

from MOSFIRE import Detector

RN = Detector.RN * u.electron
Gain = Detector.gain * u.electron/u.adu


##-------------------------------------------------------------------------
## Find Stellar Traces
##-------------------------------------------------------------------------
class TraceFitter(object):
    '''A callback object for matplotlib to make an interactive trace fitting
    window.  Contains funcationality for fitting gaussians to a 2D spectrum
    which has been collapsed in the wavelength direction.
    '''
    def __init__(self, xdata, ydata, traces=None):
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        self.xdata = xdata
        self.ydata = ydata
        self.traces = traces
        self.fitter = fitting.LevMarLSQFitter()
        if self.traces:
            self.fit = [self.traces(x) for x in self.xdata]
        else:
            self.fit = [0. for x in self.xdata]
        self.fit_line = None
        self.fit_markers = None
        self.text = None
        self.trace_table = table.Table(names=('position', 'amplitude', 'sign',\
                                              'sigma', 'FWHM'),\
                                       dtype=('f4', 'f4', 'i1', 'f4', 'f4'))

    def print_instructions(self):
        '''Print instructions for interactive fit to screen.
        '''


        text = [
        ('The figure shows the raw spectrum collapsed in the wavelength'
        ' direction as black points.  If there are one or more object traces '
        'with reasonable '
        'signal to noise, you should see one or more signals (e.g. a positive '
        'gaussian profile).'
        ),

        ('The model trace profile is a sum of '
        'gaussians and (if it exists) it is shown as a blue line.  The peak of '
        'each gaussian is marked with a red + sign.'
        ),

        ('To add a gaussian profile, to the model which will be fit to the '
        'data, put the mouse near the apex (+ sign) of that profile '
        'and press the "a" key.'
        ),

        ('To delete a gaussian profile, put the mouse near the apex (+ sign) '
        'of the profile you wish'
        ' to delete and press the "d" key.'
        ),

        ('When you are done adding or removing traces, close the interactive '
        'plot window by clicking the red close button in the upper right corner'
        ' (or by whatever method is typical for your OS or windowing system).'
        ),
        ]

        print('#'*80)
        print('                Interactive Trace Finder Instructions')
        for paragraph in text:
            print()
            print(textwrap.fill(textwrap.dedent(paragraph).strip('\n'), width=80))
        print('#'*80)
        print()

    def plot_data(self):
        '''Plot the raw data without the fit to the traces.
        '''
        self.ax.plot(self.xdata, self.ydata, 'ko-')
        plt.xlim(min(self.xdata), max(self.xdata))

    def plot_traces(self):
        '''Plot the trace fit on the raw data.
        '''
        if self.traces:
            self.fit_traces()
            self.fit = [self.traces(x) for x in self.xdata]
            xmarks = [self.traces.param_sets[3*i+1][0]\
                      for i in range(int(len(self.traces.param_sets)/3))]
            ymarks = [self.traces.param_sets[3*i][0]\
                      for i in range(int(len(self.traces.param_sets)/3))]
        else:
            self.fit = [0. for x in self.xdata]
            xmarks = []
            ymarks = []

        if not self.fit_line:
            self.fit_line, = self.ax.plot(self.xdata, self.fit, 'b-', alpha=0.7)
            self.fit_markers, = self.ax.plot(xmarks, ymarks, 'r+',\
                                             markersize=20.0,\
                                             markeredgewidth=3,\
                                             alpha=0.7)
        else:
            self.fit_line.set_ydata(self.fit)
            self.fit_markers.set_xdata(xmarks)
            self.fit_markers.set_ydata(ymarks)
        self.fig.canvas.draw()

    def fit_traces(self):
        '''Fit the trace model to the data (or update the fit based on new
        model parameters).
        '''
        if self.traces:
            self.traces = self.fitter(self.traces, self.xdata, self.ydata)
            params = self.traces.param_sets
            ntraces = int(len(params)/3)
            table_data = {'position': [params[3*i+1][0]\
                                       for i in range(ntraces)],\
                          'amplitude': [params[3*i][0]\
                                        for i in range(ntraces)],\
                          'sign': [int(params[3*i][0]/abs(params[3*i][0]))\
                                   for i in range(ntraces)],\
                          'sigma': [params[3*i+2][0]\
                                    for i in range(ntraces)],\
                          'FWHM': [2.355*params[3*i+2][0]\
                                   for i in range(ntraces)],\
                         }
            self.trace_table = table.Table(table_data)

    def connect(self):
        '''Connect keypresses to matplotlib for interactivity.
        '''
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event',\
                                                   self.keypress)
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event',\
                                                     self.click)

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_click)
        self.fig.canvas.mpl_disconnect(self.cid_key)

    def click(self, event):
        '''Print out coordinates of mouse click.  Primarily for testing.
        '''
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            print('Clicked at: {:.1f}, {:.2f}'.format(event.xdata, event.ydata))

    def keypress(self, event):
        '''Based on which key is presses on a key press event, call the
        appropriate method.
        '''
        if event.key == 'a':
            self.add_trace(event.ydata/abs(event.ydata), event.xdata)
        elif event.key == 'd':
            self.delete_trace(event)
        elif event.key == 'q':
            self.quit(event)

    def add_trace(self, amp, pos):
        '''Add an additional gaussian to the tracel model.
        '''
        if not self.traces:
            self.traces = models.Gaussian1D(mean=pos, amplitude=amp)
        else:
            self.traces += models.Gaussian1D(mean=pos, amplitude=amp)
        self.fit_traces()
        self.plot_traces()

    def delete_trace(self, event):
        '''Delete one of the gaussians from the trace model.
        '''
        x = event.xdata
        y = event.ydata
        if len(self.trace_table) > 0:
            distances = [np.sqrt( (row['amplitude']-y)**2 + (row['position']-x)**2 )
                         for row in self.trace_table]
            i = distances.index(min(distances))
            assert type(i) == int
            self.trace_table.remove_row(i)
            ## Rebuild model with only remaining traces
            self.traces = None
            for row in self.trace_table:
                self.add_trace(row['amplitude'], row['position'])
                self.fit_traces()
                self.plot_traces()

    def savefig(self, plotfile):
        '''Save the figure to a png file.
        '''
        self.fig.savefig(plotfile, bbox_inches='tight')

    def quit(self, event):
        print('quitting')
        self.disconnect()


##-------------------------------------------------------------------------
## Find Stellar Traces
##-------------------------------------------------------------------------
def find_traces(data, guesses=[], interactive=True, plotfile=None):
    '''Finds targets in spectra by simply collapsing the 2D spectra in the
    wavelength direction and fitting Gaussian profiles to positions provided.
    '''
    mdata = np.ma.MaskedArray(data=data, mask=np.isnan(data))
    vert_profile = np.mean(mdata, axis=1)
    xpoints = range(0,len(vert_profile), 1)

    if guesses == []:
        traces0 = None
    elif guesses == None:
        traces0 = None
    else:
        traces0 = models.Gaussian1D(mean=guesses[0][0], amplitude=guesses[0][1])
        if len(guesses) > 1:
            for i in range(1,len(guesses),1):
                traces0 += models.Gaussian1D(mean=guesses[i][0],\
                                             amplitude=guesses[i][1])

    tf = TraceFitter(xpoints, vert_profile, traces=traces0)
    tf.plot_data()
    tf.plot_traces()
    if interactive:
        tf.print_instructions()
        tf.connect()
        plt.show()
    if plotfile:
        tf.savefig(plotfile)

    tf.trace_table.sort('position')
    
    return tf.trace_table


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
                if verbose: print('Row {:3d}: Rejected {:d} pixels on clipping iteration {:d}'.format(\
                                  i, nrej_after-nrej_before, iter))
        
        ## Reject row if too few pixels are availabel for the fit
        if (srow.shape[0] - nrej_after) < minpixels:
            if verbose: print('Row {:3d}: WARNING! Only {:d} pixels remain after clipping and masking'.format(\
                              i, srow.shape[0] - nrej_after))
            fit = np.zeros(fit.shape)
        ## Set negative values to zero
        if np.sum((fit<0)) > 0 and verbose:
            print('Row {:3d}: Reset {:d} negative pixels in fit to 0'.format(\
                  i, np.sum((fit<0))))
        fit[(fit < 0)] = 0
        Pnew[i] = fit
    
    return Pnew


##-------------------------------------------------------------------------
## Optimal Spectral Extraction
##-------------------------------------------------------------------------
def optimal_extraction(image, variance_image, trace_table,\
                       fitsfileout=None,\
                       combine=False):
    '''Given a 2D spectrum image, a 2D variance image, and a table of trace
    lines (e.g. as output by find_traces() above), this function will optimally
    extract a 1D spectrum for each entry in the table of traces.
    
    
    '''
    if type(image) == fits.HDUList:
        hdu = image[0]
    elif type(image) == fits.PrimaryHDU:
        hdu = image
    else:
        print('Input to standard_extraction should be an HDUList or an HDU')
        raise TypeError

    if type(variance_image) == fits.HDUList:
        vhdu = variance_image[0]
    elif type(image) == fits.PrimaryHDU:
        vhdu = variance_image
    else:
        print('Input to standard_extraction should be an HDUList or an HDU')
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
        header[key] = (flattened_wcs.to_header()[key], flattened_wcs.to_header().comments[key])

    spectra = []
    variances = []
    for i,row in enumerate(trace_table):
        pos = row['position']
        width = 5.*row['sigma']
        sign = row['sign']
        print('Extracting data for trace {:d} at position {:.1f}'.format(i, pos))
        DmS = np.ma.MaskedArray(data=sign*spectra2D[pos-width:pos+width,:],\
                                mask=np.isnan(eps.data[pos-width:pos+width,:]))
        V = np.ma.MaskedArray(data=variance2D[pos-width:pos+width,:],\
                              mask=np.isnan(eps.data[pos-width:pos+width,:]))
        print('  Performing standard extraction')
        f_std, V_std = standard_extraction(DmS, V)
        print('  Forming initial spatial profile')
        P_init_data = np.array([row/f_std for row in DmS])
        P_init = np.ma.MaskedArray(data=P_init_data,\
                                   mask=np.isnan(P_init_data))
        print('  Fitting spatial profile')
        P = iterate_spatial_profile(P_init, DmS, V, f_std, verbose=False)
        print('  Calculating optimally extracted spectrum')
        f_new_denom = np.ma.MaskedArray(data=np.sum(P**2/V, axis=0),\
                                        mask=(np.sum(P**2/V, axis=0)==0))
        f_opt = np.sum(P*DmS/V, axis=0)/f_new_denom
        var_fopt = np.sum(P, axis=0)/f_new_denom
        sig_fopt = np.sqrt(var_fopt)
        spectra.append(f_opt)
        variances.append(var_fopt)
        print('  Typical level = {:.1f}'.format(np.mean(f_opt)))
        print('  Typical sigma = {:.1f}'.format(np.mean(sig_fopt[~np.isnan(sig_fopt)])))

    mask = np.isnan(np.array(spectra)) | np.isnan(np.array(variances))
    spectra = np.ma.MaskedArray(data=np.array(spectra), mask=mask)
    variances = np.ma.MaskedArray(data=np.array(variances), mask=mask)

    if not combine:
        hdulist = []
        for i,sp in enumerate(spectra):
            if i == 0:
                hdulist.append(fits.PrimaryHDU(data=sp, header=header))
            else:
                hdulist.append(fits.ImageHDU(data=sp, header=header))
        for i,var in enumerate(variances):
            hdulist.append(fits.ImageHDU(data=var, header=header))
        if fitsfileout: hdulist.writeto(fitsfileout, clobber=True)
        return hdulist
    else:
        print('Combining individual trace spectra in to final spectrum')
        spectrum = np.average(spectra, axis=0, weights=1./variances)
        sigma = np.average(1./variances, axis=0)
        variance = sigma**2
        hdulist = fits.HDUList([fits.PrimaryHDU(data=spectrum, header=header),\
                                fits.ImageHDU(data=variance, header=header)])
        if fitsfileout: hdulist.writeto(fitsfileout, clobber=True)
        return hdulist


##-------------------------------------------------------------------------
## Process some test data if called directly
##-------------------------------------------------------------------------
if __name__ == '__main__':
    ## Load MOSFIRE Data
    ## This assume you are running from within the folder which contrains the
    ## MOSFIRE_DRP_MASK/ directory which was generated by running the DRP on the
    ## "test case" data.
    path = 'MOSFIRE_DRP_MASK/2012sep10/H'
    eps_file = 'MOSFIRE_DRP_MASK_H_TARG2_eps.fits'
    sig_file = 'MOSFIRE_DRP_MASK_H_TARG2_sig.fits'
    snr_file = 'MOSFIRE_DRP_MASK_H_TARG2_snrs.fits'
    itime_file = 'MOSFIRE_DRP_MASK_H_TARG2_itime.fits'

    eps = fits.open(os.path.join(path, eps_file), 'readonly')[0]
    sig = fits.open(os.path.join(path, sig_file), 'readonly')[0]
#     snr = fits.open(os.path.join(path, snr_file), 'readonly')[0]
#     itime = fits.open(os.path.join(path, itime_file), 'readonly')[0]

    guesses = [(96, -1),\
               (124, +1),\
               (152, -1)
              ]

    trace_table = find_traces(eps.data, guesses=guesses,\
                              interactive=True, plotfile=None)
    hdulist = optimal_extraction(eps, sig, trace_table,\
                                 fitsfileout='1Dspectrum.fits',\
                                 combine=True)

    spectrum = hdulist[0].data
    w = wcs.WCS(hdulist[0].header).dropaxis(1)
    sigma = np.sqrt(hdulist[1].data)
    pix = np.arange(0,spectrum.shape[0],1)
    wavelengths = w.wcs_pix2world(pix,1)[0]
    wavelengths *= 1e6

    plt.figure(figsize=(16,6))

    plt.fill_between(wavelengths, spectrum-sigma, spectrum+sigma, label='uncertainty',\
                     facecolor='black', alpha=0.2,\
                     linewidth=0,\
                     interpolate=False)

    plt.plot(wavelengths, spectrum, 'k-', label='Final Spectrum')

#     plt.plot(pix, sigma, 'y-', alpha=0.7, label='Sigma')

    plt.title('Final Averaged Spectrum')
    plt.xlabel('Wavelength (microns)')
    plt.xlim(wavelengths.min(),wavelengths.max())
    plt.ylim(0,1.05*spectrum.max())
    plt.legend(loc='best')
    plt.savefig('final_result.png')

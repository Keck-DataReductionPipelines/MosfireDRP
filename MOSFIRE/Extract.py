#!/usr/env/python

from __future__ import division, print_function

## Import General Tools
import os
import sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import table
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
        self.trace_table = table.Table(names=('position', 'amplitude', 'sign', 'sigma', 'FWHM'),\
                                       dtype=('f4', 'f4', 'i1', 'f4', 'f4'))

    def plot_data(self):
        self.ax.plot(self.xdata, self.ydata, 'ko-')

    def plot_traces(self):
        if self.traces:
            self.fit_traces()
            self.fit = [self.traces(x) for x in self.xdata]
            xmarks = [self.traces.param_sets[3*i+1][0] for i in range(int(len(self.traces.param_sets)/3))]
            ymarks = [self.traces.param_sets[3*i][0] for i in range(int(len(self.traces.param_sets)/3))]
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

#             self.ax.text(self.trace_table[0]['position']+3.*self.trace_table[0]['sigma'],\
#                          self.trace_table[0]['amplitude'],\
#                          'Position={:.1f}'.format(self.trace_table[0]['position']))
        self.fig.canvas.draw()

    def fit_traces(self):
        if self.traces:
            self.traces = self.fitter(self.traces, self.xdata, self.ydata)
            table_data = {'position': [self.traces.param_sets[3*i+1][0]\
                                       for i in range(int(len(self.traces.param_sets)/3))],\
                          'amplitude': [self.traces.param_sets[3*i][0]\
                                        for i in range(int(len(self.traces.param_sets)/3))],\
                          'sign': [int(self.traces.param_sets[3*i][0]/abs(self.traces.param_sets[3*i][0]))\
                                   for i in range(int(len(self.traces.param_sets)/3))],\
                          'sigma': [self.traces.param_sets[3*i+2][0]\
                                    for i in range(int(len(self.traces.param_sets)/3))],\
                          'FWHM': [2.355*self.traces.param_sets[3*i+2][0]\
                                   for i in range(int(len(self.traces.param_sets)/3))],\
                         }
            self.trace_table = table.Table(table_data)

    def connect(self):
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self.keypress)
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self.click)

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_click)
        self.fig.canvas.mpl_disconnect(self.cid_key)

    def click(self, event):
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            print('Click position: {:.1f}, {:.1f}'.format(event.xdata, event.ydata))

    def keypress(self, event):
        if event.key == 'a':
            self.add_trace(event.ydata/abs(event.ydata), event.xdata)
        elif event.key == 'd':
            self.delete_trace(event)
        elif event.key == 'q':
            self.quit(event)

    def add_trace(self, amp, pos):
        if not self.traces:
            self.traces = models.Gaussian1D(mean=pos, amplitude=amp)
        else:
            self.traces += models.Gaussian1D(mean=pos, amplitude=amp)
        self.fit_traces()
        self.plot_traces()

    def delete_trace(self, event):
        x = event.xdata
        y = event.ydata
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
                traces0 += models.Gaussian1D(mean=guesses[i][0], amplitude=guesses[i][1])

    tf = TraceFitter(xpoints, vert_profile, traces=traces0)
    tf.plot_data()
    tf.plot_traces()
    if interactive:
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
    return spect_1D


##-------------------------------------------------------------------------
## Iterate Spatial Profile
##-------------------------------------------------------------------------
def iterate_spatial_profile(P, DmS, V, f, smoothing=5, order=3, minpixels=50, sigma=4, nclip=2, verbose=True):
    poly0 = models.Polynomial1D(degree=order)
    fit_poly = fitting.LinearLSQFitter()    

    Pnew = np.zeros(P.shape)
    for i,row in enumerate(P):
        weights = f**2/V[i]
        weights.mask = weights.mask | np.isnan(weights)
        srow = np.ma.MaskedArray(data=signal.medfilt(row, smoothing),\
                                 mask=(np.isnan(signal.medfilt(row, smoothing)) | weights.mask ))
        xcoord = np.ma.MaskedArray(data=np.arange(0,len(row),1),\
                                   mask=srow.mask)

        for iter in range(nclip+1):
            nrej_before = np.sum(srow.mask)
            fitted_poly = fit_poly(poly0, xcoord[~xcoord.mask], srow[~srow.mask], weights=weights[~weights.mask])
            fit = np.array([fitted_poly(x) for x in xcoord])
            resid = (DmS[i]-f*srow)**2 / V[i]
            newmask = (resid > sigma)
            
            weights.mask = weights.mask | newmask
            srow.mask = srow.mask | newmask
            xcoord.mask = xcoord.mask | newmask

            nrej_after = np.sum(srow.mask)
            if nrej_after > nrej_before:
                if verbose: print('Row {:3d}: Rejected {:d} pixels from fit on clipping iteration {:d}'.format(\
                                  i, nrej_after-nrej_before, iter))
        
        ## Reject row if too few pixels are availabel for the fit
        if (srow.shape[0] - nrej_after) < minpixels:
            if verbose: print('Row {:3d}: WARNING! Only {:d} pixels remain after sigma clipping and masking'.format(\
                              i, srow.shape[0] - nrej_after))
            fit = np.zeros(fit.shape)
        ## Set negative values to zero
        if np.sum((fit<0)) > 0 and verbose:
            print('Row {:3d}: Reset {:d} negative pixels in fit to 0'.format(i, np.sum((fit<0))))
        fit[(fit < 0)] = 0
        Pnew[i] = fit
    
    return Pnew


##-------------------------------------------------------------------------
## Optimal Spectral Extraction
##-------------------------------------------------------------------------
def optimal_extraction(spectra2D, variance2D, trace_table):
    '''Given a 2D spectrum image, a 2D variance image, and a table of trace
    lines (e.g. as output by find_traces() above), this function will optimally
    extract a 1D spectrum for each entry in the table of traces.
    
    
    '''
    spectra = []
    sigmas = []
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
        f_std = standard_extraction(DmS, V)
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
        sigmas.append(sig_fopt)
        print('  Typical level = {:.1f}'.format(np.mean(f_opt)))
        print('  Typical sigma = {:.1f}'.format(np.mean(sig_fopt[~np.isnan(sig_fopt)])))

    mask = np.isnan(np.array(spectra)) | np.isnan(np.array(sigmas))
    mspectra = np.ma.MaskedArray(data=np.array(spectra), mask=mask)
    msigmas = np.ma.MaskedArray(data=np.array(sigmas), mask=mask)
    return mspectra, msigmas


##-------------------------------------------------------------------------
## Process some test data if called directly
##-------------------------------------------------------------------------
if __name__ == '__main__':
    ## Load MOSFIRE Data
    path = '/Volumes/Internal_1TB/MOSFIRE_Data/TestCase/Reduced/MOSFIRE_DRP_MASK/2012sep10/H'
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

    trace_table = find_traces(eps.data, guesses=guesses, interactive=False, plotfile=None)
    spectra, sigmas = optimal_extraction(eps.data, sig.data, trace_table)


    print('Combining individual trace spectra in to final spectrum')
    spectrum = np.average(spectra, axis=0, weights=1./sigmas)

    plt.figure(figsize=(16,6))
    col = ['r', 'g', 'b']
    for i,row in enumerate(trace_table):
        plt.plot(spectra[i]-1*i, '{}-'.format(col[i]), alpha=0.5, label='Trace{:d} (offset by {:d})'.format(i, -1*i))
    plt.plot(spectrum, 'k-', label='Final Spectrum')
    plt.title('Final Averaged Spectrum')
    plt.xlim(0,spectrum.shape[0])
    plt.ylim(0,1.05*spectra.max())
    plt.legend(loc='best')
    plt.savefig('final_result.png')

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
from matplotlib.pylab import connect

from MOSFIRE import Detector
from MOSFIRE.MosfireDrpLog import info, debug, warning, error

RN = Detector.RN * u.electron
Gain = Detector.gain * u.electron/u.adu
py3 = sys.version_info[0] > 2 #creates boolean value for test that Python major version > 2

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

    ('The model trace profile is a sum of '
    'gaussians and (if it exists) it is shown as a blue line.  The peak of '
    'each gaussian is marked with a red + sign.'
    ),

    ('To add a gaussian profile, to the model which will be fit to the '
    'data, put the mouse near the apex of that profile '
    'and press the "a" key (for "add").'
    ),

    ('To delete a gaussian profile, put the mouse near the apex (+ sign) '
    'of the profile you wish'
    ' to delete and press the "d" key.'
    ),

    ('When you are done adding or removing traces, close the interactive '
    'plot window by clicking the red close button in the upper right corner'
    ' (or by whatever method is typical for your OS or windowing system) '
    'or press the "q" or "n" keys (for "quit" or "next" respectively).'
    ),
    ]

    print('#'*80)
    print('                Interactive Trace Finder Instructions')
    for paragraph in text:
        print()
        print(textwrap.fill(textwrap.dedent(paragraph).strip('\n'), width=80))
    print('#'*80)
    print()


class ApertureEditor(object):
    def __init__(self, xdata, ydata, title=None, traces=None):
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        self.title = title
        self.xdata = xdata
        self.ydata = ydata
        self.apertures = table.Table(names=('id', 'position', 'width',\
                                              'amplitude', 'sigma'),\
                                       dtype=('i4', 'f4', 'f4', 'f4', 'f4'))


    def add_aperture(self, pos, width):
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
        self.apertures.remove_row(id)
        self.plot_apertures()

    def set_position(self, id=None, pos=None):
        assert id is not None
        print('Changing position for aperture {} at position {:.1f}'.format(id,
              self.apertures[id]['position']))
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
        print('Changing width for aperture {} at position {:.1f}'.format(id,
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


    def fit_trace(self, pos, amp):
        id = len(self.apertures)
        g0 = models.Gaussian1D(mean=pos, amplitude=amp,
                              bounds={'amplitude': [0, float('+Inf')]})
        fitter = fitting.LevMarLSQFitter()
        g = fitter(g0, self.xdata, self.ydata)
        data = {'id': id,
                'position': g.param_sets[1][0],
                'width': np.ceil(5.*g.param_sets[2][0]),
                'amplitude': g.param_sets[0][0],
                'sigma': g.param_sets[2][0],
               }
        self.apertures.add_row(data)
        self.plot_apertures()


    def plot_data(self):
        '''Plot the raw data without apertures.
        '''
        self.ax.plot(self.xdata, self.ydata, 'ko-', label='Spatial Profile')
        plt.xlim(min(self.xdata), max(self.xdata))
        yspan = self.ydata.max() - self.ydata.min()
        plt.ylim(self.ydata.min()-0.02*yspan, self.ydata.max()+0.18*yspan)
        plt.xlabel('Pixel Position')
        plt.ylabel('Flux (e-/sec)')
        if self.title is None:
            plt.title('Spatial Profile')
        else:
            plt.title(self.title)


    def plot_apertures(self):
        plt.cla()
        self.plot_data()
        yspan = self.ydata.max() - self.ydata.min()
        for ap in self.apertures:
            if ap['sigma'] is not None and ap['amplitude'] is not None:
                g = models.Gaussian1D(mean=ap['position'],
                                      amplitude=ap['amplitude'],
                                      stddev=ap['sigma'])
                fit = [g(x) for x in self.xdata]
                self.ax.plot(self.xdata, fit, 'b-', label='Fit', alpha=0.7)
            self.ax.axvspan(ap['position']-ap['width'],
                            ap['position']+ap['width'],
                            ymin=self.ydata.min(),
                            ymax=self.ydata.max(),
                            facecolor='y', alpha=0.3,
                            )
            self.ax.text(ap['position']-ap['width']+1,
                         self.ydata.max() + 0.05*yspan,
                         'position={:.0f}\nwidth={:.0f}'.format(ap['position'],
                                                                ap['width']),
                         )
        self.fig.canvas.draw()



    def connect(self):
        '''Connect keypresses to matplotlib for interactivity.
        '''
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event',
                                                   self.keypress)

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_click)
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
            print('Adding aperture at position {}, width {}'.format(event.xdata, width))
        if event.key == 'w':
            id = self.determine_id(event)
            self.set_width(id=id)
        if event.key == 'p':
            id = self.determine_id(event)
            self.set_position(id=id)
        if event.key == 'g':
            self.fit_trace(event.xdata, event.ydata/abs(event.ydata))
        elif event.key == 'd':
            id = self.determine_id(event)
            self.delete_aperture(id)
        elif event.key == 'n':
            self.quit(event)
        elif event.key == 'q':
            self.quit(event)

    def savefig(self, plotfile):
        '''Save the figure to a png file.
        '''
        self.fig.savefig(plotfile, bbox_inches='tight')

    def quit(self, event):
        plt.close(self.fig)



class TraceFitter_orig(object):
    '''A callback object for matplotlib to make an interactive trace fitting
    window.  Contains funcationality for fitting gaussians to a 2D spectrum
    which has been collapsed in the wavelength direction.
    
    Most users will be interested in the self.traces propery which contains the
    fitted trace profiles or in the self.trace_table which is an
    astropy.tables.Table object containing properties of the Gaussians which
    have been fitted.
    
    Inputs:
        xdata : (list) The pixel numbers of the collapsed vertical profile of
                the 2D spectra.  This is usually just a simple range(#pix)
                for the ydata.
        ydata : (list) The collapsed vertical profile of the 2D spectra.
        traces : (optional, astropy.models.Gaussian1D or a model set) The
                 starting model (i.e. initial guess).  If fitting interactively,
                 then additional gaussians can be added with the "a" key.
    
    Example Use (non-interactive):
    
    Note that the non-interactive case assumes that the input traces0 contains
    an existing model (or model set) with at least 1 Gaussian profile.
    
    tf = TraceFitter(xdata, ydata, traces=traces0) # traces0 is sum of Gaussians
    tf.plot_data()   # Generate a plot of the collapsed spatial profile
    tf.fit_traces()  # Fit the model to the data
    tf.plot_traces() # Add the model points to the plot of the data
    
    The user can then examine the properties of the Gaussian model components by
    examining tf.trace_table.
    
    Example Use (interactive):
    
    tf = TraceFitter(xdata, ydata, traces=traces0) # traces0 is sum of Gaussians
    tf.plot_data()   # Generate a plot of the collapsed spatial profile
    tf.fit_traces()  # Fit the model to the data
    tf.plot_traces() # Add the model points to the plot of the data
    tf.print_instructions() # Print instructions on how to use the interactivity
    tf.connect()     # Connect the matplotlib plot to the keyboard input
    plt.show()       # Show the plot
    
    Use the "a" and "d" keys to add and delete Gaussians from the model and
    simply close the plot window when done.  The trace_table will be populated
    as in the non-interactive case.
    '''
    def __init__(self, xdata, ydata, title=None, traces=None):
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        self.xdata = xdata
        self.ydata = ydata
        self.traces = traces
        self.star_index = []
        self.title = title
        self.fitter = fitting.LevMarLSQFitter()
        if self.traces:
            self.fit = [self.traces(x) for x in self.xdata]
        else:
            self.fit = [0. for x in self.xdata]
        self.fit_line = None
        self.fit_markers = None
        self.apertures = None
        self.text = None
        self.trace_table = table.Table(names=('position', 'amplitude',\
                                              'sigma', 'width', 'FWHM'),\
                                       dtype=('f4', 'f4', 'f4', 'f4', 'f4'))


    def plot_data(self):
        '''Plot the raw data without the fit to the traces.
        '''
        self.ax.plot(self.xdata, self.ydata, 'ko-', label='Spatial Profile')
        plt.xlim(min(self.xdata), max(self.xdata))
        plt.xlabel('Pixel Position')
        plt.ylabel('Flux (e-/sec)')
        if self.title:
            plt.title(self.title)
        else:
            plt.title('Simple Spatial Profile')

    def plot_traces(self):
        '''Plot the trace fit on the raw data.
        '''
        if self.traces:
            self.fit = [self.traces(x) for x in self.xdata]
            xmarks = [self.traces.param_sets[3*i+1][0]\
                      for i in range(int(len(self.traces.param_sets)/3))]
            ymarks = [self.traces.param_sets[3*i][0]\
                      for i in range(int(len(self.traces.param_sets)/3))]
            xwidths = [3.*self.traces.param_sets[3*i+2][0]\
                       for i in range(int(len(self.traces.param_sets)/3))]
            xapertures = [[xmarks[i]-xwidths[i], xmarks[i]+xwidths[i] ]\
                          for i in range(int(len(self.traces.param_sets)/3))]
        else:
            self.fit = [0. for x in self.xdata]
            xmarks = []
            ymarks = []

        if not self.fit_line:
            self.fit_line, = self.ax.plot(self.xdata, self.fit, 'b-',\
                                          label='Fit',\
                                          alpha=0.7)
            self.fit_markers, = self.ax.plot(xmarks, ymarks, 'r+',\
                                             label='Apertures',\
                                             markersize=20.0,\
                                             markeredgewidth=3,\
                                             alpha=0.7)
            self.apertures = self.ax.axvspan(xapertures[0][0],
                                              xapertures[0][1],
                                              ymin=-200,
                                              ymax=200,
                                              facecolor='y', alpha=0.3
                                              )
        else:
            self.fit_line.set_ydata(self.fit)
            self.fit_markers.set_xdata(xmarks)
            self.fit_markers.set_ydata(ymarks)
        plt.legend(loc='best')
        self.fig.canvas.draw()

    def fit_traces(self):
        '''Fit the trace model to the data (or update the fit based on new
        model parameters).
        '''
        if self.traces:
            self.traces = self.fitter(self.traces, self.xdata, self.ydata)
            params = self.traces.param_sets
            ntraces = int(len(params)/3)
            table_data = {'position': [params[3*i+1][0]
                                       for i in range(ntraces)],
                          'amplitude': [params[3*i][0]
                                        for i in range(ntraces)],
                          'sigma': [params[3*i+2][0]
                                    for i in range(ntraces)],
                          'width': [3.*params[3*i+2][0]
                                    for i in range(ntraces)],
                          'FWHM': [2.355*params[3*i+2][0]
                                   for i in range(ntraces)],
                         }
            self.trace_table = table.Table(table_data)

    def connect(self):
        '''Connect keypresses to matplotlib for interactivity.
        '''
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event',
                                                   self.keypress)

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_click)
        self.fig.canvas.mpl_disconnect(self.cid_key)

    def keypress(self, event):
        '''Based on which key is presses on a key press event, call the
        appropriate method.
        '''
        if event.key == 'a':
            self.add_trace(event.ydata/abs(event.ydata), event.xdata)
        elif event.key == 'd':
            self.delete_trace(event)
        elif event.key == 'n':
            self.quit(event)
        elif event.key == 'q':
            self.quit(event)

    def add_trace(self, amp, pos, id=1):
        '''Add an additional gaussian to the tracel model.
        '''
        if not self.traces:
            self.traces = models.Gaussian1D(mean=pos, amplitude=amp,
                                 bounds={'amplitude': [0, float('+Inf')]})
        else:
            self.traces += models.Gaussian1D(mean=pos, amplitude=amp,
                                  bounds={'amplitude': [0, float('+Inf')]})
        self.star_index.append(id)
        self.fit_traces()
        self.plot_traces()


    def edit_width(self, amp, pos, id=1):
        '''Edit the width of the aperture.
        '''
#         if not self.traces:
#             self.traces = models.Gaussian1D(mean=pos, amplitude=amp,
#                                  bounds={'amplitude': [0, float('+Inf')]})
#         else:
#             self.traces += models.Gaussian1D(mean=pos, amplitude=amp,
#                                   bounds={'amplitude': [0, float('+Inf')]})
#         self.star_index.append(id)
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
            self.star_index = []
            for row in self.trace_table:
                self.add_trace(row['amplitude'], row['position'])
            self.plot_traces()

    def savefig(self, plotfile):
        '''Save the figure to a png file.
        '''
        self.fig.savefig(plotfile, bbox_inches='tight')

    def quit(self, event):
        plt.close(self.fig)


##-------------------------------------------------------------------------
## Find Stellar Traces
##-------------------------------------------------------------------------
def find_traces(data, guesses=[], title=None, interactive=True, plotfile=None):
    '''Finds targets in spectra by simply collapsing the 2D spectra in the
    wavelength direction and fitting Gaussian profiles to positions provided.
    '''
    mdata = np.ma.MaskedArray(data=data, mask=np.isnan(data))
    vert_profile = np.mean(mdata, axis=1)
    xpoints = range(0,len(vert_profile), 1)
    ap = ApertureEditor(xpoints, vert_profile, title=title)

    if (guesses == []) or (guesses is None):
        ## Try to guess at trace positions if no information given
        ## Start with maximum pixel value
        valid_profile = list(vert_profile[~vert_profile.mask])
        maxval = max(valid_profile)
        maxind = valid_profile.index(maxval)
        ap.fit_trace(maxind, maxval)
    else:
        for guess in guesses:
            if len(guess) == 1:
                ap.fit_trace(guess, 1)
            else:
                ap.add_aperture(guess[0], guess[1])
    ap.plot_apertures()
    if interactive:
        ap.connect()
        plt.show()
    return ap.apertures




def find_traces_orig(data, guesses=[], title=None, interactive=True, plotfile=None):
    '''Finds targets in spectra by simply collapsing the 2D spectra in the
    wavelength direction and fitting Gaussian profiles to positions provided.
    '''
    mdata = np.ma.MaskedArray(data=data, mask=np.isnan(data))
    vert_profile = np.mean(mdata, axis=1)
    xpoints = range(0,len(vert_profile), 1)

    if (guesses == []) or (guesses is None):
        ## Try to guess at trace positions if no information given
        ## Start with maximum pixel value, but only accept it if the
        ## neighboring pixels are also high
        valid_profile = list(vert_profile[~vert_profile.mask])
        maxval = max(valid_profile)
        maxind = valid_profile.index(maxval)
        traces0 = models.Gaussian1D(mean=maxind, amplitude=maxval)
    else:
        traces0 = models.Gaussian1D(mean=guesses[0][0], amplitude=guesses[0][1])
        if len(guesses) > 1:
            for i in range(1,len(guesses),1):
                traces0 += models.Gaussian1D(mean=guesses[i][0],\
                                             amplitude=guesses[i][1])

    tf = TraceFitter(xpoints, vert_profile, traces=traces0, title=title)
    tf.plot_data()
    tf.fit_traces()
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
def optimal_extraction(image, variance_image, trace_table,
                       fitsfileout=None,
                       plotfileout=None,
                       plot=None):
    '''Given a 2D spectrum image, a 2D variance image, and a table of trace
    lines (e.g. as output by find_traces() above), this function will optimally
    extract a 1D spectrum for each entry in the table of traces.
    
    
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
    for i,row in enumerate(trace_table):
        pos = row['position']
        width = row['width']
        info('Extracting data for trace {:d} at position {:.1f}'.format(i, pos))
        DmS = np.ma.MaskedArray(data=spectra2D[int(pos-width):int(pos+width),:],\
                                mask=np.isnan(spectra2D[int(pos-width):int(pos+width),:]))
        V = np.ma.MaskedArray(data=variance2D[int(pos-width):int(pos+width),:],\
                              mask=np.isnan(spectra2D[int(pos-width):int(pos+width),:]))
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
        spectra.append(f_opt)
        variances.append(var_fopt)
        info('  Typical level = {:.1f}'.format(np.mean(f_opt)))
        info('  Typical sigma = {:.1f}'.format(np.mean(sig_fopt[~np.isnan(sig_fopt)])))

    mask = np.isnan(np.array(spectra)) | np.isnan(np.array(variances))
    spectra = np.ma.MaskedArray(data=np.array(spectra), mask=mask)
    variances = np.ma.MaskedArray(data=np.array(variances), mask=mask)

    w = wcs.WCS(header).dropaxis(1)
    wh = w.to_header()
    for key in wh.keys():
        header[key] = wh[key]

    for i,row in enumerate(trace_table):
        hdulist = fits.HDUList([])
        if plotfileout:
            fig = plt.figure(figsize=(16,6))
            wavelength_units = getattr(u, w.to_header()['CUNIT1'])

        sp = spectra[i]
        hdulist.append(fits.PrimaryHDU(data=sp.filled(0), header=header))
        hdulist[0].header['TRACEPOS'] = row['position']
        if plotfileout:
            sigma = 1./variances[i]
            pix = np.arange(0,sp.shape[0],1)
            wavelengths = w.wcs_pix2world(pix,1)[0] * wavelength_units.to(u.micron)*u.micron
            plt.subplot(len(trace_table), 1, i+1)
            plt.fill_between(wavelengths, sp-sigma, sp+sigma,\
                             label='uncertainty',\
                             facecolor='black', alpha=0.2,\
                             linewidth=0,\
                             interpolate=True)
            plt.plot(wavelengths, sp, 'k-',
                     label='Spectrum for Trace {} at {}'.format(i, row['position']))
            plt.xlabel('Wavelength (microns)')
            plt.ylabel('Flux (e-/sec)')
            plt.xlim(wavelengths.value.min(),wavelengths.value.max())
            plt.ylim(0,1.05*sp.max())
            plt.legend(loc='best')

            bn, ext = os.path.splitext(plotfileout)
            plotfilename = '{}_{:02d}{}'.format(bn, i, ext)
            plt.savefig(plotfilename, bbox_inches='tight')
            plt.close(fig)

        var = variances[i]
        hdulist.append(fits.ImageHDU(data=var.filled(0), header=header))
        hdulist[1].header['TRACEPOS'] = row['position']
        hdulist[1].header['COMMENT'] = 'VARIANCE DATA'
        if fitsfileout:
            bn, ext = os.path.splitext(fitsfileout)
            fitsfilename = '{}_{:02d}{}'.format(bn, i, ext)
            hdulist.writeto(fitsfilename, clobber=True)

    return hdulist


##-------------------------------------------------------------------------
## Extract Spectra Function
##-------------------------------------------------------------------------
def extract_spectra(maskname, band, interactive=True):

    ## Get objectnames from slit edges
    edges = np.load('slit-edges_{}.npy'.format(band))
    objectnames = [edge['Target_Name'] for edge in edges[:-1]]

    trace_tables = {}
    if interactive:
        print_instructions()
    for objectname in objectnames:
        eps_file = '{}_{}_{}_eps.fits'.format(maskname, band, objectname)
        eps = fits.open(eps_file, 'readonly')[0]
#         trace_plot_file = '{}_{}_{}_trace.png'.format(maskname, band, objectname)
        trace_plot_file = None
        print('Finding traces for {}'.format(objectname))
        trace_tables[objectname] = find_traces(eps.data, title=objectname,
                                               interactive=interactive,
                                               plotfile=trace_plot_file)

    for objectname in objectnames:
        eps_file = '{}_{}_{}_eps.fits'.format(maskname, band, objectname)
        sig_file = '{}_{}_{}_sig.fits'.format(maskname, band, objectname)
        eps = fits.open(eps_file, 'readonly')[0]
        sig = fits.open(sig_file, 'readonly')[0]
        spectrum_plot_file = '{}_{}_{}.png'.format(maskname, band, objectname)
        fits_file = '{}_{}_{}_1D.fits'.format(maskname, band, objectname)
        info('Extracting traces for {}'.format(objectname))
#         try:
        hdulist = optimal_extraction(eps, sig, trace_tables[objectname],
                                     fitsfileout=fits_file,
                                     plotfileout=spectrum_plot_file,
                                     )
#         except:
#             warning('Failed to extract trace for {}'.format(objectname))

##-------------------------------------------------------------------------
## Process some test data if called directly
##-------------------------------------------------------------------------
if __name__ == '__main__':
    extract_spectra('MOSFIRE_DRP_MASK', 'H', interactive=True)

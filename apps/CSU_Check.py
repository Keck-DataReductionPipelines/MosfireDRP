'''

Written March 3rd 2011 by npk
'''

import numpy as np
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
from matplotlib import pyplot as pl
from pyraf import iraf
from MOSFIRE import CSU, Detector, IO, Fit


def fit_line_with_sigclip(xs, data, i = 0):

        ps = np.polyfit(xs, data, 1)
        pf = np.poly1d(ps)

        residual = np.abs(pf(xs) - data)
        sd = np.std(residual)
        
        ok = np.where(residual < 2*sd)[0]

        ps = np.polyfit(xs[ok], data[ok], 1)
        pf = np.poly1d(ps)
        return [pf, ok]

        if len(ok) == len(residual):
                return [pf, ok]
        elif i > 2:
                return [pf, ok]
        else:
                return fit_line_with_sigclip(xs[ok], data[ok], i+1)

def median_tails(v):
        a = np.median(v[0:2])
        b = np.median(v[-3:-1])

        t = v - np.float(a+b)/2.

        return t

def make_slice(pos, w, h):
        '''Returns [ [xslice, ylsice], [x0, x1, y0, x1] ] where
        xslice is used as Array[x0:x1]
        yslice is used as Array[y0:y1]'''
        x0 = pos[0]-w
        if x0 < 0: x0 = 0
        x1 = pos[0] + w
        if x1 > 2047: x1 = 2047
        if x0 > x1: x0 = x1
        xs = slice(x0, x1)

        y0 = pos[1]-h
        if y0 < 0: y0 = 0
        y1 = pos[1]+h
        if y1 > 2047: y1 = 2047
        if y0 > y1: y0 = y1
        ys = slice(y0,y1)

        return [[xs,ys],[x0,x1,y0,y1]]

(header, data) = IO.readfits("/users/npk/desktop/c8/m101029_0233.ref.fits")
(header2, data2) = IO.readfits("/users/npk/desktop/c8/m101029_0425.ref.fits")
(header3, data3) = IO.readfits("/users/npk/desktop/c8/m101029_0427.ref.fits")

data = data3


deg = np.pi/180.
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

reload(CSU)
reload(IO)
reload(Fit)
reload(Detector)

pl.ion()

bs = CSU.Barset()
bs.set_header(header)

pl.figure(1)
pl.clf()

means = []
sds = []
deltas = []
poss = []
qs = []

cnt = 1

cntfit = 1
pl.figure(5)
pl.clf()
pl.subplot(7, 7, cntfit)
pl.figure(6)
pl.clf()
pl.subplot(7, 7, cntfit)

for bar in range(4, 92, 2):
        print bar
        pos = bs.get_bar_pix(bar)
        [[xslice, yslice],extent] = make_slice(pos, 6,25)
        if extent[0] == extent[1]: 
                cnt += 1
                continue
        if extent[2] == extent[3]: 
                cnt += 1
                continue
        cnt+=1

        fits = []
        xs = np.arange(-10,10)
        for i in xs:
                tofit = data[pos[1]-i, xslice]
                y = median_tails(tofit)

                ps = Fit.do_fit(y, Fit.residual_pair)
                fits.append(ps[0])

        fits = np.array(fits)

        m = [np.mean(fits[:,i]) for i in range(5)]
        s = [np.std(fits[:,i]) for i in range(5)]
        means.append(m)
        sds.append(s)

        [ff, ok] = fit_line_with_sigclip(xs, fits[:,1])

        pl.figure(5)
        pl.subplot(7,7,cntfit)
        pl.plot(xs, fits[:,1] - ff(xs), '*-')
        pl.plot(xs[ok], fits[ok,1] - ff(xs[ok]), 'or-')
        pl.ylim([-.1,.1])
        pl.title("%2i" % bar)

        pl.figure(6)
        pl.subplot(7,7,cntfit)
        pl.plot(xs, fits[:,4],'*-')
        pl.plot(xs[ok], fits[ok,4],'or-')
        pl.ylim([2.9,4])

        cntfit += 1

        delta = (extent[0] + ff[0]) - pos[0]
        poss.append(extent[0] + ff[0])
        deltas.append(delta)
        q = np.degrees(ff[1]) - np.degrees(CSU.rotation)
        qs.append(q)

        pl.figure(1)
        pl.text(pos[0], pos[1], 'b%2.0i: w=%3.2f p=%5.2f q=%3.2f d=%1.3f' % (bar, np.mean(fits[:,4]), extent[0]+ff[0], q, delta), fontsize=11, family='monospace', horizontalalignment='center')

pl.xlim([0,2048])
pl.ylim([0,2048])


means = np.array(means)
sds = np.array(sds)
pl.draw()

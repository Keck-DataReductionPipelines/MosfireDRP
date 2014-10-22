'''

Written March 18th 2011 by npk


'''

import sys, datetime, getpass, os
import numpy as np, pyfits as pf
import pylab as pl
import scipy.ndimage.filters, scipy.io
from pyraf import iraf
from MOSFIRE import CSU, Detector, IO, Fit

from IPython.Shell import IPShellEmbed

start_shell = IPShellEmbed()

NaN = np.nan


def save_region(bs, fname):
        s = bs.to_ds9_region()
        try:
                f = open(fname, "w")
                f.writelines(s)
                f.close()
        except:
                return Exception("Could not write barset region file to: %s" % fname)

def fit_line_with_sigclip(xs, data, i=0):


        ps = Fit.do_fit_edge(xs, data)
        pf = lambda x: Fit.fit_bar_edge(ps, x)


        residual = np.abs(pf(xs) - data)
        sd = np.std(residual)

        ok = np.where(residual < 2.5 * sd)[0]
        if len(ok) == 0:
                return [lambda x: NaN, []]

        ps = Fit.do_fit_edge(xs[ok], data[ok])
        pf = lambda x: Fit.fit_bar_edge(ps, x)

        return [pf, ok]



if False:
        def fit_line_with_sigclip(xs, data, i = 0):

                ps = np.polyfit(xs, data, 1)
                pf = np.poly1d(ps)


                residual = np.abs(pf(xs) - data)
                sd = np.std(residual)
                
                ok = np.where(residual < 2.5*sd)[0]

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

def make_slice(pos, offset, w, h):
        '''Returns [ [xslice, ylsice], [x0, x1, y0, x1] ] where
        xslice is used as Array[x0:x1]
        yslice is used as Array[y0:y1]'''
        x0 = round(pos[0]- w + offset)
        if x0 < 0: x0 = 0
        x1 = round(pos[0] + w + offset)
        if x1 > 2047: x1 = 2047
        if x0 > x1: x0 = x1
        xs = slice(x0, x1)

        y0 = round(pos[1]-h)
        if y0 < 0: y0 = 0
        y1 = round(pos[1]+h)
        if y1 > 2047: y1 = 2047
        if y0 > y1: y0 = y1
        ys = slice(y0,y1)


        return [[xs,ys],[x0,x1,y0,y1],round(pos[1])]

deg = np.pi/180.
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

reload(CSU)
reload(IO)
reload(Fit)
reload(Detector)


def sigclip(data, low=4, high=4):
        c = data.ravel()
        delta = 1
        while delta:
                s = c.std()
                m = c.mean()
                size = c.size

                c = c[(c > (m - s*low)) & (c < (m + s*high))]

                delta = size-c.size

        return c.mean()


def is_in_bounds(extent):
        if extent[1] > 2047: return False
        if extent[0] < 0: return False
        if extent[2] > 2047: return False
        if extent[0] < 0: return False
        if extent[0] == extent[1]: 
                return False
        if extent[2] == extent[3]: 
                return False

        for i in [0,1]:
                for j in [2,3]:
                        if not CSU.in_field(extent[i],extent[j]): return False
        return True

def is_odd(n):
        return (n % 2) == 1


def go(fname):
        global plot_arr

        print fname
        
        (header, data) = IO.readfits(fname)
        bs = CSU.Barset()
        bs.set_header(header)

        bars = []
        means = []
        sds = []
        deltas = []
        deltas_mm = []
        poss = []
        poss_mm = []
        poss_y = []
        request = []
        qs = []
        bars = range(1, 93)
        fit_fun = Fit.residual_single

        for bar in bars:
                pos = bs.get_bar_pix(bar)
                if bar % 8 == 0:
                        print "%2.0i: (%7.2f, %7.2f)" % (bar, pos[0], pos[1])


                if is_odd(bar):
                        if (bs.pos[bar] - bs.pos[bar-1]) < 2.7:
                                fit_fun = Fit.residual_pair
                        else:
                                fit_fun = Fit.residual_single

                width = 19 
                [[xslice, yslice],extent,ystart] = make_slice(pos,0,width,30)

                

                if not is_in_bounds(extent):
                        fits = [0,0,0,0,0]
                        [ff,ok] = [np.poly1d(0,0), []]
                        means.append(fits)
                        sds.append(fits)
                        drop_this = True
                else:
                        drop_this = False
                        fits = []
                        ys = np.arange(-10,10, dtype=np.float32)
                        for i in ys:
                                tofit = data[ystart-i, xslice]
                                y = median_tails(tofit)

                                ps = Fit.do_fit(y, fit_fun)
                                fits.append(ps[0])
                        
                        fits = np.array(fits)
                        fits[:,1] += 1

                        # fit to the ridgeline
                        [ff, ok] = fit_line_with_sigclip(ys, fits[:,1])
                        m = [np.mean(fits[:,i]) for i in range(5)]
                        s = [np.std(fits[:,i]) for i in range(5)]
                        means.append(m)
                        sds.append(s)


                slit_center_offset = pos[1] - ystart
                fc = ff(slit_center_offset)
                slit_center_pos = np.float(extent[0] + fc )

                if drop_this: 
                        poss.append(NaN)
                        poss_y.append(NaN)
                        poss_mm.append(NaN)
                else: 
                        poss.append(slit_center_pos)
                        poss_y.append(ystart)
                        poss_mm.append(CSU.csu_pix_to_mm_poly(slit_center_pos, ystart)[0])

                delta = np.float(slit_center_pos - pos[0])
                if drop_this: 
                        deltas.append(NaN)
                        deltas_mm.append(NaN)
                else: 
                        deltas.append(delta)
                        b = CSU.csu_pix_to_mm_poly(slit_center_pos + delta, ystart)[0]
                        deltas_mm.append(b - poss_mm[-1])

                q = np.float(np.degrees(np.tan(ff(1)-ff(0))))
                if drop_this: qs.append(NaN)
                qs.append(q)


        means = np.array(means)
        f = lambda x: np.array(x).ravel()
        sds = f(sds)
        deltas = f(deltas)
        poss = f(poss)
        poss_y = f(poss_y)
        poss_mm = f(poss_mm)
        deltas_mm = f(deltas_mm)
        qs  = f(qs)
        bars = f(bars)



        fout = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".sav"
        print "saving"
        tosav = {"bars": bars, "request": bs.pos, "deltas_mm": deltas_mm, "poss": poss, "poss_mm": poss_mm, "deltas": deltas, "means": means, "qs": qs}
        scipy.io.savemat(fout, tosav)
        save_region(bs, "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".reg")
        print "saved"

        regout = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".meas.reg"
        pairs = np.array([poss,poss_y]).transpose()
        s = CSU.to_ds9_region(pairs, dash=0, color="blue", label=False)
        try:
                f = open(regout, "w")
                f.writelines(s)
                f.close()
        except:
                print "Couldn't write: %s" % regout


        return [tosav, bs]



path  = "/users/npk/desktop/c9/"
def generate_fname(num):
        global path
        files = os.listdir(path)
        for fn in files:
                if fn.find("12_%4.4i" % num) > 0:
                        return fn
                if fn.find("13_%4.4i" % num) > 0:
                        return fn
                if fn.find("16_%4.4i" % num) > 0:
                        return fn
                if fn.find("17_%4.4i" % num) > 0:
                        return fn
                if fn.find("18_%4.4i" % num) > 0:
                        return fn
                if fn.find("19_%4.4i" % num) > 0:
                        return fn
                if fn.find("20_%4.4i" % num) > 0:
                        return fn
                if fn.find("21_%4.4i" % num) > 0:
                        return fn
                if fn.find("24_%4.4i" % num) > 0:
                        return fn

[ts, bs] = go(path + generate_fname(2846))

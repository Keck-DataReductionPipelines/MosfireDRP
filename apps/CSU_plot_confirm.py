
import numpy as np
import pylab as pl
import scipy as sp
import scipy.io
import os
from MOSFIRE import CSU
from matplotlib.backends.backend_pdf import PdfPages

from IPython.Shell import IPShellEmbed

start_shell = IPShellEmbed()

def fit_line_with_sigclip(xs, data, i = 0):

        ps = np.polyfit(xs, data, 1)
        pf = np.poly1d(ps)


        residual = np.abs(pf(xs) - data)
        sd = np.std(residual)
        
        ok = np.where(residual < 3.0*sd)[0]

        ps = np.polyfit(xs[ok], data[ok], 1)
        pf = np.poly1d(ps)
        return [pf, ok]

path = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" 
proto = "m11031%1.1i_%4.4i.fits.sav.mat"

pl.ion()
rs = []
for i in range(1754, 1793):
        if os.path.exists(path + proto % (6, i)):
                ld = scipy.io.loadmat(path + proto % (6, i))
        elif os.path.exists(path + proto % (7, i)):
                ld = scipy.io.loadmat(path + proto % (7, i))
        elif os.path.exists(path + proto % (8, i)):
                ld = scipy.io.loadmat(path + proto % (8, i))
        elif os.path.exists(path + proto % (9, i)):
                ld = scipy.io.loadmat(path + proto % (9, i))
        else:
                print "No luck"
                continue

        ld["img_num"] = i
        rs.append(ld)



bar_posns = [[] for i in range(92)]

for r in rs:
        for i in range(92):
                assert(r["bars"][i] == i+1)
                p = r["poss_mm"][i]
                d = r["deltas_mm"][i]

                if not np.isfinite(p): continue
                if not np.isfinite(d): continue
                if p < 20: continue
                if p > 200: continue
                bar_posns[i].append([p, d, r["bars"][i]])


even = np.arange(1, 92, 2)
odd = np.arange(0, 92, 2)

ds = []
ff = np.poly1d([.08/45, -0.04])
pl.figure(1)
pl.clf()
fits = []
for r in rs:
        pm = r["poss_mm"]
        ok = np.where(np.isfinite(pm))[0]
        if len(ok) < 1: continue

        b = r["bars"][ok]
        #pl.plot(r["deltas_mm"][ok] + 0.0416 - 0.0010392*b,'o')
        pl.plot(r["deltas_mm"][ok],'o')
        fits.append(np.polyfit(r["bars"][ok].ravel(), r["deltas_mm"][ok].ravel(), 1))
        #ds.extend(r["deltas_mm"][ok] + 0.0416 - 0.0010392*b)
        ds.extend(r["deltas_mm"][ok])


ds = np.array(ds).ravel()

# Measured by hand
groups = [75, 100, 125, 150, 172, 196]
groups = [25, 50, 75, 100, 125, 150, 172, 196, 220, 244]

pl.figure(2)
pl.clf()
pl.axis("scaled")
pl.ylim(-140,140)
pl.xlim(-140,140)
pl.xlabel("Keck Field Position (mm)")
pl.ylabel("Keck Field Position (mm)")
pl.grid()

scale = 300
mn = ds.mean()
for bar in range(92):
        bv = np.array(bar_posns[bar])
        if len(bv) == 0: continue
        pos = bv[:,0]
        delts = bv[:,1]
        yshift = (bar % 2)*1.
        for group in groups:
                p = np.where(np.isfinite(pos) & (np.abs(pos-group)<5))[0]
                position = pos[p].mean() - mn - 137.4
                delta = delts[p].mean() - mn

                ypos = (np.round(bar/2.)-23.0) * 5.8 * CSU.tempscale + yshift
                if delta > 0: pl.arrow(position, ypos, delta*scale, 0, color='red', head_width=.8)
                else: pl.arrow(position, ypos, delta*scale, 0, color='blue',head_width=.8)

pl.text(-132,-117,"0.1 arcsecond")
pl.arrow(-120,-120,scale*.725*.1,0, head_width=.9)

pl.figure(3)
pl.clf()
arc = (ds-mn)/.725
pl.hist(arc,40,color='w')
pl.title("(Achieved - Requested) Bar Positions")
pl.xlabel("Bar Offsets [arcsecond]")
pl.ylabel("N (total = %i)" % len(arc))
pl.text(-0.1, 50, "SD: %3.2f arcsecond\nP-V: %3.2f arcsecond" % (arc.std(), arc.max()-arc.min()))

print 
print "Mean offset of: %4.3f" % mn

pdf = PdfPages("CSU_confirmation_mar_16_2011.pdf")
pl.figure(2).savefig(pdf, format="pdf")
pl.figure(3).savefig(pdf, format="pdf")
pdf.close()

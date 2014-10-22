
from MOSFIRE import CSU
import pylab as pl, numpy as np
from numpy import sin, cos

reload(CSU)

pl.ion()

pl.figure(1)
pl.clf()
pl.xlim([-100,2200])
pl.ylim([-200,2200])

dxs = []
for s in range(4):
        for x in range(0,281,280):
                p1 = CSU.csu_mm_to_pix_poly(x, s)
                p2 = CSU.csu_mm_to_pix(x,s,Please_Use=True)

                dxs.append(p1[0]-p2[0])

                dx = p1[0]-p2[0]
                dy = p1[1]-p2[1]
                scale=50
                pl.arrow(p1[0], p1[1], dx*scale, dy*scale)
                if s == 0:
                        pl.axvline(p2[0],color='black',ls='-.',lw=.5)



pl.figure(2)
pl.clf()


try:
        f = open("../platescale/finalshift.xyXrotYrot.4.972.120k.dat")
        lines = f.readlines()
        f.close()
except:
        print "Failed to read file"


pl.figure(3)
pl.clf()
pl.xlim([-100,2200])
pl.ylim([-100,2200])
for line in lines:
        [px,py,xmm,ymm]= map(np.float,line.split())

        px_lin = 7.670726 * xmm
        py_lin = 7.670726 * ymm

        d = np.radians(360-359.762)
        px_lin = cos(d) * px_lin - sin(d) * py_lin
        py_lin = sin(d) * px_lin + cos(d) * py_lin

        px_lin += 1042.986
        py_lin += 1035.879
        
        dx = px-px_lin
        dy = py-py_lin

        scale =200
        print "%3.2f %3.2f" % (dx, dy)
        pl.arrow(px_lin, py_lin, scale*dx, scale*dy)
        pl.plot([px_lin],[py_lin],'ok')

pl.arrow(100,1800,scale/2.,0)
pl.title("Distortion Map: ccs solution - linear solution")

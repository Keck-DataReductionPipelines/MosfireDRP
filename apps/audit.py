

import time
import traceback
import getpass
import os
import pdb
import pprint
import sets
import sqlite3
import sys
import textwrap

import numpy as np
from matplotlib import pyplot as pl

from operator import itemgetter
from itertools import groupby

import MOSFIRE 

from MOSFIRE import Options, IO, Wavelength


def audit(filename):
    
    header, data = IO.readfits(filename)


    ll0 = header['crval1']
    dlam = header['cd1_1']

    ls = ll0 + dlam * np.arange(data.shape[1])

    linelist = Wavelength.pick_linelist(header)


    deltas = []
    sigs = []
    xpos = []
    ys = []
    for y in np.linspace(750, 1100, 30):
    #for y in np.linspace(5, 640, 50):
        sp = np.ma.array(data[y,:])

        xs, sxs, sigmas = Wavelength.find_known_lines(linelist, ls, sp,
                Options.wavelength)

        xpos.append(xs)
        ys.append([y] * len(xs))

        deltas.append(xs - (linelist - ll0)/dlam)
        sigs.append(sxs)


    xpos, ys, deltas, sigs = map(np.array, [xpos, ys, deltas, sigs])

    deltas[np.abs(deltas) > .75] = np.nan
    sigs[np.abs(sigs) > .001] = np.nan

    pl.clf()
    size = 0.003/sigs 
    size[size > 30] = 30
    size[size < 1] = 1
    pl.scatter( xpos, ys, c=deltas, s=size)
    pl.xlim([0, data.shape[1]])
    pl.ylim([0, data.shape[0]])
    pl.xlabel("Spectral pixel")
    pl.ylabel("Spatial pixel")
    pl.title("Night sky line deviation from solution [pixel]")
    pl.colorbar()

    pl.savefig("audit.pdf")
    
    pdb.set_trace()


def usage():
    print """
    audit [filename]
Commands: """

    print("\n")

if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
        sys.exit()

    audit(sys.argv[-1])

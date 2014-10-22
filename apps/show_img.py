#!/usr/bin/python

'''
Written  March 15th 2011 by npk
'''


import os, numpy, scipy, sys
path  = "/users/npk/desktop/c9/"
def generate_fname(num):
        global path
        files = os.listdir(path)
        for fn in files:
                if fn.find("12_%4.4i" % num) > 0:
                        return fn
                if fn.find("13_%4.4i" % num) > 0:
                        return fn
                if fn.find("14_%4.4i" % num) > 0:
                        return fn
                if fn.find("15_%4.4i" % num) > 0:
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



fn = generate_fname(int(sys.argv[1]))
print fn

regpath = "/users/npk/dropbox/mosfire/cooldown\ 9/csu_meas/"
os.system("ds9 %s -regions %s -regions %s" %(path+fn, regpath+fn+".reg", regpath+fn+".meas.reg"))
#os.system("ds9 %s -regions %s" %(path+fn, regpath+fn+".reg"))


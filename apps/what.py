#!/usr/local/bin/python

'''
MOSFIRE 'what' command:
    Spits out an informative summary of files
    in the current directory. Or files selected
    via a glob.

    i.e. what *0311.fits

npk - March 23 2011
'''

import MOSFIRE
import MOSFIRE.IO
import glob
import sys



files = []
if len(sys.argv) == 1:
    files = glob.iglob('*')
else:
    for i in range(1, len(sys.argv)):
        files.extend(glob.iglob(sys.argv[i]))



#print("filename          object  exptime        maskname lamp  filt   Turret")
for fname in files:

    try:
        header = MOSFIRE.IO.readheader(fname)
    except IOError, err:
        print("Couldn't IO %s" % fname)
        continue
    except:
        print("%s is unreadable" % fname)
        continue

    lamps = ""
    try:
        if header["pwstata7"] == 1:
            lamps += header["pwloca7"][0:2]
        if header["pwstata8"] == 1:
            lamps += header["pwloca8"][0:2]
    except KeyError:
        lamps = "???"
        
    header.update("lamps", lamps)

    try:
        if header["aborted"]:
            header.update("object", "ABORTED")
    except:
        print("Missing header file in: %s" % fname)

    try:
        print("%(datafile)12s %(object)40s %(truitime)6.1f s %(maskname)35s %(lamps)3s %(filter)4s %(mgtname)7s" % (header))
    except:
        try:
            print("%(datafile)12s %(object)25s %(truitime)6.1f s %(lamps)3s %(filter)6s %(mgtname)7s" % (header))
        except:
            print("%s Skipped" % fname)


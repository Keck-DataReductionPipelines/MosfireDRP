

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

from operator import itemgetter
from itertools import groupby

import MOSFIRE 

from MOSFIRE import Options, IO




def load_db():
    indir = Options.indir
    outname = os.path.join(Options.outdir, "mosfire_files.db")

    print("Database: {0}".format(outname))
    conn = sqlite3.connect(outname)

    return conn

def create(cursor):
    cursor.execute('''
        CREATE TABLE if not exists files
        (id integer primary key, path text, fdate text, number integer)
    ''')

keys = []

def append_column(cursor, name, typename):
    qry = "alter table files\nadd {0} {1}".format(name, typename)
    try:
        cursor.execute(qry)
        print("Added {0} as {1}".format(name, typename))
    except sqlite3.OperationalError:
        pass


def make():
    """Make the database"""

    db = load_db()
    c = db.cursor()
    create(c)
    dirs = os.walk(Options.indir)

    Options.indir = Options.indir.rstrip("/")

    for root, dirs, files in dirs:
        if root == Options.indir: continue
        ignore, path = root.split(Options.indir)

        if len(path.split("/")) != 2: continue

        try: date = int(path.split("/")[1][0:4])
        except: continue

        if (date < 2012) or (date > 2030): continue

        for file in files:
            if len(file) != 17: continue
            p = os.path.join(root, file)

            num = db.execute('select count(*) from files where path = "%s"' %
                    p).fetchall()
            if num[0][0] > 0: 
                print("Skipping: " + p + " [already in db]")
                continue
            print(p)

            hdr = IO.readheader(p)

            try:
                fdate = file.split("_")[0][1:]
                number = file.split("_")[1][:-5]
            except:
                print("Skipping: " + p)
                continue

            insert_sql = "insert into files(path,fdate,number,"
            vals = "?,?,?,"
            values = [p, fdate, number]

            for key in hdr.keys():

                if key == 'COMMENT': continue

                value = hdr[key]
                T = type(value)
                key = key.replace("-","_")

                insert_sql += key + ","
                vals += "?,"
                values.append(value)

                if key in keys: continue
                keys.append(key)
            


                if T == int: typename = 'integer'
                if T == float: typename = 'real'
                else: typename = 'text'
                append_column(c, key, typename)


            insert_sql = insert_sql[:-1] + ") values (" + vals[:-1] + ")"
            try:
                c.execute(insert_sql, tuple(values))
            except:
                print "Query failed on:"
                print insert_sql
                traceback.print_exc()
                #sys.exit()
                 

    db.commit()

def find_continuous(data):
    '''Find all continuous numbers in a list'''
    # http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    ranges = []
    for k, g in groupby(enumerate(data), lambda (i,x):i-x):
        group = map(itemgetter(1), g)
        ranges.append((group[0], group[-1]))
    return ranges

def underline_ansi(str):
    return chr(27) + '[4m' + str + chr(27) + '[0m'

def bold_ansi(str):
    return chr(27) + '[1m' + str + chr(27) + '[0m'
      
def boldunderline_ansi(str):
    return chr(27) + '[1m' + chr(27) + '[4m' + str + chr(27) + '[0m'


def sql_for_mask_group_filter(db, maskname):
    cur = db.execute(
    '''
    select count(filter), filter, itime/1000.0, yoffset 
    from files 

    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy"

    group by filter'''.
                format(maskname))

    return cur.fetchall()

def sql_for_mask_filter_flats(db, maskname, filter):

    query = '''
    select path, fdate, number
    from files
    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy" and
    filter = "{1}" and (el-45) < .1 and flatspec = 1
    order by fdate, number
            '''.format(maskname, filter)

    print "Flat Query is:", query
    cursor = db.execute(query)
    return cursor.fetchall()

def sql_for_mask_filter_spectra(db, maskname, filter):
    query = '''
    select fdate
    from files
    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy" and
    filter = "{1}" and (itime/1000.0) > 30 and flatspec = 0 and (domestat =
    "tracking" or domestat = 0) and aborted = 0

    group by fdate

            '''.format(maskname, filter)
    
    #print "DB Query is: ", query
    cur = db.execute(query)
    return cur.fetchall()

def sql_for_mask_filter_date(db, maskname, filter, date):
    query = '''
    select path, fdate, number, yoffset, itime/1000.0
    from files
    where maskname = "{0}" and filter = "{1}" and (itime/1000.0) > 30 and 
            fdate = {2} and flatspec = 0  and (domestat = "tracking" or
            domestat = 0) and aborted = 0
    order by fdate, number
    '''.format(maskname, filter, date)

    print "DB Query is: ", query
    cur = db.execute(query)

    return cur.fetchall()
 
def plan_to_fname(plan):
    return "%s_%s.py" % (plan["maskname"], plan["filter"])

longslit_plan_file ='''
# This file was automatically generated by the mospy db application
#  The application was executed by {uid} on {createdate}
#
# Help, bugs to: http://mosfire.googlecode.com
#
# Instructions
#   1. edit band = 'fixme' to band = 'Y' or 'J' or 'H' or 'K'
#       e.g. band = 'J'
#   2. edit range(a,b) to be a list of flat names
#   3. edit range(c,d) to be a list of long names
#       Note for steps 2&3 most likely these will be a continuous sequence
#   4. edit [709, 1350] to be the pixel values at the beginning and end 
#       of the long slit. Look at the raw data.

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify
from MOSFIRE import Wavelength, Longslit

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = '{maskname}'
band = '{band}'

flatnames = {flatnames}
longnames = {longnames}

flatops = Options.flat
waveops = Options.wavelength

{lslitoptions}

Flats.handle_flats(flatnames, maskname, band, flatops)
Wavelength.imcombine(longnames, maskname, band, waveops)
Wavelength.fit_lambda_interactively(maskname, band, longnames, waveops)
Wavelength.fit_lambda(maskname, band, longnames, longnames, waveops,
    longslit=longslit)
Wavelength.apply_lambda_simple(maskname, band, longnames, waveops,
    longslit=longslit, smooth=True)

Longslit.go(maskname, band, longnames, waveops, longslit)



'''

plan_file ='''
# This file was automatically generated by the mospy db application
#  The application was executed by {uid} on {createdate}
#
# Help, bugs to: http://mosfire.googlecode.com

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify
from MOSFIRE import Wavelength

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = '{maskname}'
band = '{band}'
num_dates = {num_dates}

flatnames = {flatnames}

sciframes = {sciframes}

wavenames = {wavenames}

flatops = Options.flat
waveops = Options.wavelength

Flats.handle_flats(flatnames, maskname, band, flatops)

{wavecombine}

Combine.handle_combine(wavenames, maskname, band, waveops)


'''

def plan_to_python(plans):
    '''Convert the python list/dictionary created by masks() into a python
    script'''


    '''
    A plan is a structure that looks something like

    Plan {
        filter -> Filter name (string)
        maskname -> maskname (string)
        flatlist -> 
                    ["YYmmDD_####"...] list of flats
        dates -> [{
            date -> YYmmDD (string)
            observations -> [{
                observation -> (a,b) (tuple of file number range)
                offsets [{
                    "name" -> ["YYmmDD_####" ...] list of sci frames at offset
                                                    "name"
                }]
            }]
        }]
    }


    This function unpacks the above structure into a python program that will
    produce a data reduction plan file.
    '''

    
    for plan in plans:
        fname = plan_to_fname(plan)
        if os.path.exists(fname):
            #print("Plan '%s' already exists, remove the plan file " 
                    #"to overwrite" % fname)
            os.remove(fname)
            #REMOVE COMMENT BELOW:
            #continue

        outf = open(fname, "w")
        num_dates = len(plan["dates"])

        waves = []
        scis = []
        for date in plan["dates"]:  
            for observation in date["observations"]:
                obs_wave = []
                obs_sci = {}
                offsets = observation["offsets"].keys()


                if (len(offsets) == 1) and offsets[0] is 'Unknown':
                    fnames = observation["offsets"]['Unknown']['fname']
                    obs_sci["A"] = fnames[0:-1:2]
                    obs_sci["B"] = fnames[1:-1:2]
                    obs_wave.extend(fnames)
                else:
                    for offset in offsets:
                        fnames = observation["offsets"][offset]['fname']
                        obs_sci[offset] = fnames
                        obs_wave.extend(fnames)

                scis.append(obs_sci)
                waves.append(obs_wave)


        wavecombine = ""
        for i in xrange(len(waves)):
            wavecombine += "Wavelength.imcombine(wavenames[%i], maskname, " \
                "band, waveops)\n" % (i)
            if i == 0:
                wavecombine += "Wavelength.fit_lambda_interactively(" \
                    "maskname, band, wavenames[0], waveops)\n"

            wavecombine += "Wavelength.fit_lambda(" \
                    "maskname, band, wavenames[%i], wavenames[0], " \
                    " waveops)\n" % i

            wavecombine += "Wavelength.apply_lambda_simple(maskname, band, " \
                    " wavenames[%i], waveops)\n" % i

            pos = scis[i].keys()
            if len(pos) != 2:
                print "Only handling A/B subtraction currently"
                continue

            wavecombine += \
                    "Background.handle_background(sciframes[%i]['%s'], " \
                    "sciframes[%i]['%s'], wavenames[%i], maskname, band, " \
                    "waveops)\n" % (i, pos[0], i, pos[1], i)

            wavecombine += \
                    "Rectify.handle_rectification(maskname, ['A', 'B'], " \
                    "wavenames[%i], band, waveops)" % (i)

            wavecombine += "\n"


        res = { "uid": getpass.getuser(), 
                "createdate": time.asctime(),
                "maskname": plan["maskname"], 
                "band": plan["filter"],
                "flatnames": plan["flatlist"], 
                "sciframes": scis,
                "wavenames": waves, 
                "wavecombine": wavecombine, 
                "num_dates": num_dates}

        outf.write(plan_file.format(**res))
        outf.close()
        

def longslits():
    """List all longslits"""

    if len(sys.argv) == 4:
        db = load_db()
        fdate = int(sys.argv[3])

        query = """
            select object, path, fdate, number, filter, yoffset, maskname, 
                gratmode, itime, el
            from files
            where substr(maskname,0,9) == 'LONGSLIT' and fdate = "{0}"
            order by number
            """.format(fdate, fdate)
        
        cur = db.execute(query)

        ress = cur.fetchall()

        if len(ress) == 0:
            raise Exception("No such objects in database. Query:\n%s" % query)
            return

        print("{0}".format(ress[0][-1]))
        print("{0}".format(object))
        print("{0:6s} {1:6s} {2:3s} {3:6s} {4:4s} {5:15s}".format("type", "date", "num", 
            "band", "offset", "object"))
        objs = {}
        for res in ress:
            obj, path, fdate, number, filter, yoffset, maskname, gratmode, exptime, el = res
            guess = '?'
            if gratmode == 'imaging':
                guess = "align"
            elif filter == 'Dark':
                guess = 'dark'
            elif filter == 'Moving':
                guess = 'bad'
            elif len(obj) > 4 and obj[0:4] == 'Flat':
                guess = 'flat'
                key = "flat_{0}".format(filter)
            else:
                guess = "sci"
                key = "{0}_{1}".format(obj,filter)

            if guess == 'flat' or guess == 'sci':
                if objs.has_key(key):
                    objs[key].append(path)
                else:
                    objs[key] = [path]

            if res[5] is None: offset = -999
            else: offset = float(res[5])
            print("{0:6s} {1:6s} {2:3g} {3:6s} {4:5.1f} {5:15s}".format(guess, res[2], 
                res[3], res[4], offset, obj))

        print("")
        print("--- SUMMARY ---")
        for key, value in objs.iteritems():
            print("{0:10s}: {1:5g} frames".format(key, len(value)))


    else:
        print "Not enough arguments"
        sys.exit()

    res = {
        "uid": getpass.getuser(), 
        "createdate": time.asctime(),
        'maskname': "longslit_%s" % (fdate),
        'band' : 'fixme',
        'flatnames': "['m%s_%%4.4i.fits' %% i for i in range(a,b)]" % (fdate),
        'longnames': "['m%s_%%4.4i.fits' %% i for i in range(c,d)]" % (fdate),
        'lslitoptions': "longslit = {'yrange': [709, 1350]}"
    }
    
    fout = "%s_longslit.py" % fdate
    try:
        f = open(fout, "w")
        f.write(longslit_plan_file.format(**res))
        f.close()
    except:
        print("Could not open and write to {0}".format(fout))
    

def masks():
    """List all slit masks"""

    db = load_db()

    if len(sys.argv) == 3:
        cur = db.execute("select maskname, count(maskname) from files group by maskname")
        ress = cur.fetchall()
        print("{0:74s} {1:5s}".format("Mask Name", "Count"))
        print("-"*80)
        bold_on = False
        for res in ress:
            output = "{0:74s} {1:5g}".format(res[0], res[1])
            if bold_on: print bold_ansi(output)
            else: print output
            bold_on = not bold_on

        print
        print('''Execute:
        mospy db masks [maskname]
        to generate a mask plan''')

    if len(sys.argv) == 4:
        maskname = sys.argv[3]
       
        FILTERS = sql_for_mask_group_filter(db, maskname)

        plans = []
        for res in FILTERS:

            num_frames, filter, itime, yoffset = res
            if yoffset is None: yoffset='Unknown'
            this_plan = {"maskname": maskname, "filter": filter}

            print
            print(boldunderline_ansi("{0:45s} {1:4s}".format(maskname, filter)))

            if filter == 'Dark':
                print "  Dark frames not fully supported yet"
                continue

            FL = sql_for_mask_filter_flats(db, maskname, filter)

            print "%i flats on %i nights " % (len(FL), len(set([str(S[1]) for
                S in FL])))

            this_plan["flatlist"] = [str("m%s_%4.4i.fits" % (f[1],f[2])) for f
                    in FL]
            

            DATES = sql_for_mask_filter_spectra(db, maskname, filter)

            this_plan["dates"] = []
            for date in DATES:
                date = date[0]
                this_date = {"date": date}
                              
                FRAMES = sql_for_mask_filter_date(db, maskname, filter, date)

                print(underline_ansi("{0}: {1} frames:".format(date,
                    len(FRAMES))))

                nums = [int(S[2]) for S in FRAMES]
                observations = find_continuous(nums)

                this_date["observations"] = []
                for observation in observations:
                    this_observation = {"observation": observation}
                    offsets = {}
                    for frame in FRAMES:
                        path, fdate, number, yoffset, itime = frame
                        if yoffset is None: yoffset = "Unknown"

                        if (number < observation[0]) or (number > 
                                observation[1]):
                            continue

                        if float(yoffset) == 0: pdb.set_trace()
                        if offsets.has_key(yoffset):

                            offsets[yoffset]["fname"].append(
                                    str("m%s_%4.4i.fits" % (fdate,number)))

                            offsets[yoffset]["itime"] += itime
                        else:
                            offsets[yoffset] = {}

                            offsets[yoffset]["fname"] = [
                                    str("m%s_%4.4i.fits" %
                                    (fdate, number))]

                            offsets[yoffset]["itime"] = itime
                            offsets[yoffset]["start/stop"] = observation

                        this_observation["offsets"] = offsets
                    this_date["observations"].append(this_observation)

                for observation in this_date["observations"]:
                    for k,v in observation["offsets"].iteritems():
                        print("\tOffset {0:5s} has {1:3g} frames ({2}-{3}) "
                            "total exptime is {4:5g} s".format(str(k),
                                len(v["fname"]), v["start/stop"][0],
                                v["start/stop"][1], v["itime"]))

                this_plan["dates"].append(this_date)

            plans.append(this_plan)
            plan_to_python(plans) 




                

commands = [make, masks, longslits]
def usage():
    print """
Commands: """

    for command in commands:
        print("\t" + command.__name__ + ": " + command.__doc__)

    print("\n")

if __name__ == '__main__':

    if len(sys.argv) < 3:
        usage()
        sys.exit()

    if sys.argv[2] == 'make':
        print "Making database"
        make()
    if sys.argv[2] == 'masks':
        masks()
    if sys.argv[2] == 'longslits':
        longslits()

    else:
        usage()
        sys.exit()




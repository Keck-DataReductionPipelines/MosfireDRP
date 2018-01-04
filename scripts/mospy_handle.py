#!/usr/local/bin/python

'''
MOSFIRE 'handle' command:

(c) npk - Dec 2013
'''
import MOSFIRE.IO as IO
import os
import numpy as np
import sys
import glob

from MOSFIRE.MosfireDrpLog import debug, info, warning, error

if len(sys.argv) < 3:
    print('''Usage: mospy handle [target]''')
    sys.exit()

## Output the file list to a text file for later examination
if os.path.exists('filelist.txt'):
    debug('Removing old filelist.txt')
    os.remove('filelist.txt')
fl = open('filelist.txt', 'w')


files = []
for i in range(1, len(sys.argv)):
    files.extend(glob.iglob(os.path.abspath(sys.argv[i])))

files = [file for file in files if os.path.splitext(file)[1] != '.original']

masks = {}


info('Examining {} files'.format(len(files)))
for fname in files:

    try:
        header = IO.readheader(fname)
    except IOError:#, err:
        fl.write("Couldn't IO %s\n" % fname)
        continue
    except:
        fl.write("%s is unreadable\n" % fname)
        continue

    lamps = ""
    try:
        if header["pwstata7"] == 1:
            lamps += header["pwloca7"][0:2]
        if header["pwstata8"] == 1:
            lamps += header["pwloca8"][0:2]
    except KeyError:
        lamps = "???"
        
    header['lamps'] = lamps

    try:
        if header["aborted"]:
            header['object' ] = 'ABORTED'
    except:
        fl.write("Missing header file in: %s\n" % fname)

    try:
        fl.write("%(datafile)12s %(object)35s %(truitime)6.1fs %(maskname)35s %(lamps)3s %(filter)4s %(mgtname)7s\n" % (header))
    except:
        try:
            fl.write("%(datafile)12s %(object)25s %(truitime)6.1fs %(lamps)3s %(filter)6s %(mgtname)7s\n" % (header))
        except:
            fl.write("%s Skipped\n" % fname)
            continue


    datafile = header['datafile'] + '.fits'
    maskname = str(header['maskname'])
    target = str(header['targname'])
    filter = header['filter']
    yr,mn,dy = IO.fname_to_date_tuple(datafile)
    date = str(yr)+mn+str(dy)
    object = header['object']
    frameid = header['FRAMEID'].strip()

    itime = header['truitime']
    grating_turret = header['mgtname']

    if object.find("MIRA") == -1: 
        mira = False
    else: 
        mira = True

    if header['MGTNAME'] is not 'mirror':
        mira = False
        
    if maskname.find(" (align)") == -1:
        align = False
    else:
        maskname = maskname.replace(" (align)", "")
        align = True

    if maskname.find('LONGSLIT') != -1:
#         print("longslit file")
        align = False

    if maskname.find('long2pos') != -1:
        if grating_turret != 'mirror':
            align = False

    empty_files = {'Align': [], 'Ne': [], 'Ar': [], 'Flat': [], 'FlatThermal': [],
            'Dark': [], 'Aborted': [], 'Image': [], 'MIRA': [], 'Unknown': []}

    if maskname not in masks:
        masks[maskname] = {date: {filter: empty_files}}

    if date not in masks[maskname]:
        masks[maskname][date] = {filter: empty_files}

    if filter not in masks[maskname][date]:
        masks[maskname][date][filter] = empty_files

    # convert numbers such as 1.0 to 1, but leaves 1.5 as 1.5
    # - added to match AutoDriver.py code
    offset_hdr = float(header['YOFFSET'])
    if offset_hdr % 1 == 0:
        offsetvalue = int(offset_hdr)
    else:
        offsetvalue = offset_hdr

    offset = 'Offset_' + str(offsetvalue)
    if (maskname.find('long2pos') != -1 and align is False) or maskname.find('LONGSLIT') != -1:
        # if the target name contains a /, replace it with _
        target_name = target.replace("/","_")
        # if the target name contains a space, remove it
        target_name = target_name.replace(" ","")
        # add a posC and posA to the offset names
        position = ''
        if header['XOFFSET']>0:
            position = 'PosC'
        if header['XOFFSET']<0:
            position = 'PosA'
        offset = offset+'_'+str(target_name)
        if position is not '':
            offset = offset+'_'+position

    
    if mira:
        masks[maskname][date][filter]['MIRA'].append(fname)
    elif align:
        masks[maskname][date][filter]['Align'].append(fname)
    elif 'Ne' in header['lamps']:
        masks[maskname][date][filter]['Ne'].append(fname)
    elif 'Ar' in header['lamps']:
        masks[maskname][date][filter]['Ar'].append(fname)
    elif header['ABORTED']:
        masks[maskname][date][filter]['Aborted'].append(fname)
    elif header['FILTER'] == 'Dark':
        masks[maskname][date][filter]['Dark'].append(fname)
    elif header['FLATSPEC'] == 1:
        masks[maskname][date][filter]['Flat'].append(fname)
    elif object.find("Flat:") != -1 and ( object.find("lamps off") != -1 or object.find("Flat:Off")) != -1 :
        masks[maskname][date][filter]['FlatThermal'].append(fname)
    elif header['mgtname'] == 'mirror':
        masks[maskname][date][filter]['Image'].append(fname)
    elif offset != 0:
#         print "offset is now:"+str(offset)
        if frameid in ["A", "B", "A'", "B'","D","C", "E"]:
            if offset in masks[maskname][date][filter]: 
                masks[maskname][date][filter][offset].append((fname, itime))
#                 print("adding file to existing offset file")
            else: 
                masks[maskname][date][filter][offset] = [(fname, itime)]
#                 print("creating new offset file")
        else:
            fl.write('{} has unexpected FRAMEID: {}\n'.format(fname, frameid))
    else:
        masks[maskname][date][filter]['Unknown'].append(fname)



##### Now handle mask dictionary

def descriptive_blurb():
    import getpass, time

    uid = getpass.getuser()
    date = time.asctime()

    return "# Created by '%s' on %s\n" % (uid, date)


# Write out the list of files in filepath
#   list = ['/path/to/mYYmmDD_####.fits' ...]
#   filepath is absolute path to the file name to write to
#
#   Result, is a file called filepath is written with
# fits files in the list.
def handle_file_list(output_file, files):
    '''Write a list of paths to MOSFIRE file to output_file.'''

    if os.path.isfile(output_file):
        print("%s: already exists, skipping" % output_file )
#         pass

    if len(files) > 0:
        with open(output_file, "w") as f:
            f = open(output_file, "w")
            f.write(descriptive_blurb())

            picker = lambda x: x
            if len(files[0]) == 2: picker = lambda x: x[0]

            # Identify unique path to files:
            paths = [os.path.dirname(picker(file)) for file in files]
            paths = list(set(paths))

            if len(paths) == 1:
                path_to_all = paths[0]
                converter = os.path.basename
                f.write("%s # Abs. path to files [optional]\n" % path_to_all)
            else:
                converter = lambda x: x

            info('Writing {} files to {}'.format(len(files), output_file))
            for path in files:
                if len(path) == 2:
                    to_write = "%s # %s s\n" % (converter(path[0]), path[1])
                else:
                    to_write = "%s\n" % converter(path)
                f.write("%s" % to_write)


def handle_date_and_filter(mask, date, filter, mask_info):

    path = os.path.join(mask,date,filter)
    try:
        os.makedirs(path)
    except OSError:
        pass

    for type in list(mask_info.keys()):
        handle_file_list(os.path.join(path, type + ".txt"), mask_info[type])


for mask in list(masks.keys()):
    for date in list(masks[mask].keys()):
        for filter in list(masks[mask][date].keys()):
            handle_date_and_filter(mask, date, filter, masks[mask][date][filter])



#!/usr/local/bin/python                                                                                                                                                    

import MOSFIRE
from MOSFIRE import IO, Wavelength
from MOSFIRE.IO import fname_to_path
import os
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import time
import sys
import glob

class Driver:

    def __init__(self,outputFile,type):
        self.outputFile = outputFile
        self.type = type
        self.offsetFiles = []
        allowedTypes = ['slitmask', 'longslit', 'long2pos', 'long2pos_specphot']
        if self.type not in allowedTypes:
            print "Unknown driver type"
        else:
            print "Generating automatic driver file "+outputFile
            self.target = open(outputFile,'w')
            self.import_section()

    def addLine(self, line):
        self.target.write(line+"\n")

    def import_section(self):
        self.addLine("import os, time, logging")
        self.addLine("import MOSFIRE")
        self.addLine("from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength, Extract")
        self.addLine("from MOSFIRE.MosfireDrpLog import info, debug, warning, error")
        self.addLine("logger = logging.getLogger(__name__)")
        self.addLine("import numpy as np")
        self.addLine("from matplotlib import pyplot as pl")
        self.addLine("from astropy.io import fits as pf")
        self.addLine("np.seterr(all='ignore')")
        self.addLine("flatops = Options.flat")
        self.addLine("waveops = Options.wavelength")
        self.addLine("")

    def addOffsetFiles(self,offsetFiles, resetList=False):
        # might not be needed
        if resetList:
            self.offsetFiles = []
        for offsetfile in offsetFiles:
            self.offsetFiles.append(offsetfile)

    def printObsfiles(self,obsfiles):
        for obsfile in obsfiles:
            self.addLine(obsfile)
        self.addLine("")

    def printnoninteractive(self,noninteractive=False):
        self.addLine("#Set noninteractive to True to autofit wavelenth solution instead of manually fitting.")
        self.addLine("noninteractiveflag="+str(noninteractive))

    def printMaskAndBand(self):
        offsetfile = self.offsetFiles[0]
        fname = IO.list_file_to_strings(offsetfile)

        if os.path.isabs(fname[0]): path = fname[0]
        else: path = os.path.join(fname_to_path(fname[0]), fname[0])
        hdulist = pf.open(path)
        header = hdulist[0].header

        self.maskName = header['maskname']
        self.band = header['filter']
        self.addLine("maskname = '"+str(self.maskName)+"'")
        self.addLine("band = '"+str(self.band)+"'")
        self.addLine("")
                     
    def isEmpty(self,file):
        fname = IO.list_file_to_strings(file)
        if len(fname):
            return False
        else:
            return True

    def printFlat(self):
        longslit=""
        if self.type is 'long2pos' or self.type is 'long2pos_specphot' or self.type is 'longslit':
            longslit=",longslit=longslit"
    
        # using only Flat.txt
        if os.path.isfile('Flat.txt'):
            if self.isEmpty('Flat.txt') is True:
                self.addLine("### WARNING: Flat.txt is empty! ###")
            flatLine = "Flats.handle_flats('Flat.txt', maskname, band, flatops"+longslit+")"

        # using both Flat.txt and FlatThermal.txt
        if os.path.isfile('FlatThermal.txt') and self.band is 'K':
            if self.isEmpty('FlatThermal.txt') is True:
                self.addLine("### WARNING: FlatThermal.txt is empty! ###")
            flatLine = "Flats.handle_flats('Flat.txt', maskname, band, flatops,lampOffList='FlatThermal.txt'"+longslit+")"

        # write the resulting line
        self.addLine(flatLine)
        self.addLine("")        

    def addLongslit(self):
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
            self.addLine("# Note: for long2pos, the row position is ignored, and the middle point of the slit is used")
            self.addLine("longslit = {'yrange': [[1062,1188],[887,1010]], 'row_position': 0, 'mode':'long2pos'}")
        if self.type is 'longslit':
            # use the slitname to determine the range (such as LONGSLIT-3x0.7)
            numberOfSlits = int(self.maskName.lstrip("LONGSLIT-").split("x")[0])
            verticalOffset = 10 # this is the vertical offset to apply to each measurement to shift the position up in the detector. It seems to be around 10
            slitSizePixels = int(numberOfSlits*(2048/46))
            slitTop = 1024+slitSizePixels/2+verticalOffset
            slitBottom = 1024-slitSizePixels/2+verticalOffset
            RowPosition = 1024+verticalOffset
            self.addLine("longslit = {'yrange':["+str(slitBottom)+","+str(slitTop)+"],'row_position':"+str(RowPosition)+",'mode':'longslit'}")
            
        
    def printWavelengthFit(self):
        if self.type is 'longslit' or self.type is 'long2pos':
            addLongSlit = ",longslit=longslit"
        else:
            addLongSlit = ""
            
        
        if self.type is 'slitmask' or self.type is 'longslit':

            self.useNeon = False
            self.useArgon = False
            # determine is Argon and Neon files contain data for K bands
            if self.isEmpty('Ar.txt') is False and self.band is 'K':
                self.useArgon = True
            if self.isEmpty('Ne.txt') is False and self.band is 'K':
                self.useNeon = True

            self.addLine("Wavelength.imcombine(obsfiles, maskname, band, waveops)")
            if self.useArgon:
                self.addLine("Wavelength.imcombine('Ar.txt', maskname, band, waveops)")
            if self.useNeon:
                self.addLine("Wavelength.imcombine('Ne.txt', maskname, band, waveops)")

            self.addLine("Wavelength.fit_lambda_interactively(maskname, band, obsfiles,waveops"+addLongSlit+", noninteractive=noninteractiveflag)")

            if self.useArgon:
                self.addLine("Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ar.txt', argon=True)")
            if self.useNeon:
                self.addLine("Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ne.txt', neon=True)")

            self.addLine("Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,waveops"+addLongSlit+")")

            if self.useArgon and self.useNeon:
                self.addLine("Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt',waveops, wavenames2='Ar.txt'"+addLongSlit+")")
            if self.useArgon and not self.useNeon:
                self.addLine("Wavelength.fit_lambda(maskname, band, 'Ar.txt', 'Ar.txt',waveops"+addLongSlit+")")
            if self.useNeon and not self.useArgon:
                self.addLine("Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt',waveops"+addLongSlit+")")

            if self.useNeon or self.useArgon:
                self.addLine("LROI = [[21000,22800]]*1")
            if self.useNeon:
                self.addLine("LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles, 'Ne.txt', LROI, waveops)")
            if self.useArgon and not self.useNeon:
                self.addLine("LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles, 'Ar.txt', LROI, waveops)")

            self.addLine("Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops"+addLongSlit+")")

            if self.useArgon and self.useNeon:
                self.addLine("Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ne.txt', LROIs, waveops)")
            if self.useArgon and not self.useNeon:
                self.addLine("Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ar.txt', LROIs, waveops)")
            if self.useNeon and not self.useArgon:
                self.addLine("Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ne.txt', LROIs, waveops)")

            # determine waveleng name
            files = IO.list_file_to_strings(self.offsetFiles)
            if self.useNeon:
                neon_files = IO.list_file_to_strings('Ne.txt')
                self.waveName = "merged_lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,"")).rstrip(".fits")+"_and_"+str(Wavelength.filelist_to_wavename(neon_files, self.band, self.maskName,""))
            elif self.useArgon and not self.useNeon:
                argon_files = IO.list_file_to_strings('Ar.txt')
                self.waveName = "merged_lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,"")).rstrip(".fits")+"_and_"+str(Wavelength.filelist_to_wavename(argon_files, self.band, self.maskName,""))           
            else:
                self.waveName = "lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,""))            
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
            calibWith = ""
            if self.isEmpty('Ar.txt') is False: 
                self.addLine("argon = ['Ar.txt']")
                calibWith = "argon"
                waveFiles = IO.list_file_to_strings('Ar.txt')
            if self.isEmpty('Ne.txt') is False:
                self.addLine("neon = ['Ne.txt']")
                calibWith = "neon"
                waveFiles = IO.list_file_to_strings('Ne.txt')            
            if calibWith:
                # we have either Argon, or Neon, or both, so we can use arcs for the reduction
                self.addLine("Wavelength.imcombine("+str(calibWith)+", maskname, band, waveops)")
                self.addLine("Wavelength.fit_lambda_interactively(maskname, band, "+str(calibWith)+",waveops,longslit=longslit, "+str(calibWith)+"=True, noninteractive=noninteractiveflag)")
                self.addLine("Wavelength.fit_lambda(maskname, band, "+str(calibWith)+","+str(calibWith)+",waveops,longslit=longslit)")
                self.addLine("Wavelength.apply_lambda_simple(maskname, band, "+str(calibWith)+", waveops, longslit=longslit, smooth=True)")            
                self.waveName = "lambda_solution_"+str(Wavelength.filelist_to_wavename(waveFiles, self.band, self.maskName,""))
            else:
                # we have no arcs. For the time being, we can try with sky lines but this only works with long2pos specphot
                print "#####################################################################################################"
                print "WARNING: There are no arc calibration files"
                print "         The pipeline will try to use sky lines but this only works if the observation is long enough"
                print "         and if you are only using long2pos. It will NOT work on long2pos_specphot"
                print "         Please contact the MosfireDRP team to obtain a standard wavelength solution"
                print "#####################################################################################################" 
                self.addLine("obsfiles = obsfiles_posAnarrow + obsfiles_posCnarrow")
                self.addLine("Wavelength.imcombine(obsfiles, maskname, band, waveops)")
                self.addLine("Wavelength.fit_lambda_interactively(maskname, band, obsfiles ,waveops,longslit=longslit, noninteractive=noninteractiveflag)")
                self.addLine("Wavelength.fit_lambda(maskname, band, obsfiles,obsfiles ,waveops,longslit=longslit)")
                self.addLine("Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops, longslit=longslit, smooth=True)")            
                files = IO.list_file_to_strings(self.offsetFiles)
                self.waveName = "lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,""))                
                
        self.addLine("")
        self.addLine("Wavelength_file = '"+str(self.waveName)+"'")
        self.addLine("")
            
    def printBackground(self):
        if self.type is 'long2pos_specphot':
            for slit in ['posAnarrow','posCnarrow','posAwide','posCwide']:
                self.addLine("Background.handle_background(obsfiles_"+str(slit)+",Wavelength_file,maskname,band,waveops, target=target_"+str(slit)+")")
        if self.type is 'long2pos':
            for slit in ['posAnarrow','posCnarrow']:
                self.addLine("Background.handle_background(obsfiles_"+str(slit)+",Wavelength_file,maskname,band,waveops, target=target_"+str(slit)+")")
        if self.type is 'slitmask':
            self.addLine("Background.handle_background(obsfiles,Wavelength_file,maskname,band,waveops)")
        if self.type is 'longslit':
            self.addLine("Background.handle_background(obsfiles,Wavelength_file,maskname,band,waveops,target=target)")

        self.addLine("")

    def printRectification(self):
        if self.type is 'slitmask': 
            self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles]')
            self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles,waveops)')
        if self.type is 'longslit':
            self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles]')
            self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles,waveops, target=target)')
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
            for slit in ['posAnarrow','posCnarrow']:
                self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles_'+str(slit)+']')            
                self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_'+str(slit)+',waveops, target=target_'+str(slit)+')')
        if self.type is 'long2pos_specphot': 
            for slit in ['posAwide','posCwide']:
                self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles_'+str(slit)+']')
                self.addLine('redfiles = [redfiles[0]]')
                self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_'+str(slit)+',waveops, target=target_'+str(slit)+')')
        self.addLine("")

    def printExtraction(self):
        self.addLine('Extract.extract_spectra(maskname, band, interactive=True, combine=True)')


    def printHeader(self):
        now = time.strftime("%c")
        self.addLine("#Driver file automatically generated on "+str(now))
        self.addLine("#For questions and comments, email mosfiredrp@gmail.com, submit a ticket on the ticketing system, or contact Luca Rizzi @ WMKO")
        self.addLine("")
                     
        
        
    def CloseFile(self):
        self.target.close()



def OffsetPairs():

    offsetFiles = glob.glob("Offset_*.txt")
    # remove Offset_ and remote .txt
    tmpOffsets = [off.replace("Offset_","").replace(".txt","") for off in offsetFiles]
    # for each name, separate using _ as a separator
    slitmaskOffset = []
    processedTargets= []
    targets_and_offsets= {}
    for off in tmpOffsets:
        # separate using _
        off_array = off.split('_')
        # the first value of the array is the offset value

        # if the array has only one element, this is a slitmask (Offset_1.5.txt), add this value to the slitmask offsets
        if len(off_array)==1:
            type = 'slitmask'
            if "slitmask" in targets_and_offsets:
                tmp = targets_and_offsets["slitmask"]
                tmp.append(float(off_array[0]))
                targets_and_offsets["slitmask"]=tmp
            else:
                targets_and_offsets["slitmask"]=[float(off_array[0]),]
        else:
            # if the array has more than one element, we are in a long2pos or longslit mode
            # if the last element is a PosC or PosA, then it's long2pos
            #print off_array
            if off_array[-1] in ['PosA','PosC']:
                type = 'long2pos'
                # check if we have already seen this target
                tname = "_".join(off_array[1:-1])
                # we are doing this for cases in which the file is Offset_-7_HIP87_7.25_PosA.txt (the _ in the file name is a problem)
            else:
                type = 'longslit'
                tname = "_".join(off_array[1:])

            if tname not in processedTargets:
                # add the new target to the list
                processedTargets.append(tname)
                # add the current offset to the list of offsets files for this target
            if tname in targets_and_offsets:
                tmp=targets_and_offsets[tname]
                tmp.append(float(off_array[0]))
                #print "adding new offset to target "+str(tname)
                targets_and_offsets[tname]=tmp
            else:
                #print "creating new offset set for target "+str(tname)
                targets_and_offsets[tname]=[float(off_array[0]),]

    return targets_and_offsets,type

def SetupFiles(target=None, offsets=None, type=None):
    # convert numbers such as 1.0 to 1, but leaves 1.5 as 1.5
    offsets = [int(f) if f % 1 ==0 else f for f in offsets]
    setupLines = []
    obsFiles = []
    specphot = False
    type=type

    # slitmask
    if type is 'slitmask':
        offsets = [f for f in offsets if f>0]
        for off in offsets:
            obsFiles.append("Offset_"+str(off)+".txt")
            obsFiles.append("Offset_"+str(off*-1)+".txt")
        setupLines.append("obsfiles=['"+str("','".join(obsFiles))+"']")
    elif type is 'longslit':
        # files are assumed to be in pairs, and we drop the "0" value is present.
        # remove negative and 0 offsets
        offsets = [f for f in offsets if f>0]
        for off in offsets:
            obsFiles.append("Offset_"+str(off)+"_"+str(target)+".txt")
            obsFiles.append("Offset_"+str(off*-1)+"_"+str(target)+".txt")

        setupLines.append("obsfiles=['"+str("','".join(obsFiles))+"']")
        setupLines.append('target="'+str(target)+'"')
    elif type is 'long2pos' or type is 'long2pos_specphot':
        # old long 2 pos (-7,-14,-21, 7,14,21)
        # narrow slits
        if set([-7,-21,7,21]).issubset(offsets):
            setupLines.append("obsfiles_posCnarrow = ['Offset_-21_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            obsFiles.append("Offset_7_"+str(target)+"_PosA.txt")  # we are using this to determine maskname and band
            obsFiles.append("Offset_-7_"+str(target)+"_PosC.txt")  # we are using this to determine maskname and band            
            setupLines.append('target_posCnarrow = "'+str(target)+'_POSC_NARROW"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posCnarrow)")
            setupLines.append("obsfiles_posAnarrow = ['Offset_7_"+str(target)+"_PosA.txt', 'Offset_21_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAnarrow = "'+str(target)+'_POSA_NARROW"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posAnarrow)")
        # wide slits
        if set([-7,-14,7,14]).issubset(offsets):
            setupLines.append("obsfiles_posCwide = ['Offset_-14_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            setupLines.append('target_posCwide = "'+str(target)+'_POSC_WIDE"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posCwide)") 
            setupLines.append("obsfiles_posAwide = ['Offset_14_"+str(target)+"_PosA.txt', 'Offset_21_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAwide = "'+str(target)+'_POSA_WIDE"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posAwide)") 
            specphot = True
        # new long 2 pos (-7,0,7)
        # narrow slits
        if set([-7,7]).issubset(offsets) and not(set([21,21]).issubset(offsets)):
            setupLines.append("obsfiles_posCnarrow = ['Offset_7_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            obsFiles.append("Offset_7_"+str(target)+"_PosA.txt")
            obsFiles.append("Offset_-7_"+str(target)+"_PosC.txt")  # we are using this to determine maskname and band
            setupLines.append('target_posCnarrow = "'+str(target)+'_POSC_NARROW"')
            setupLines.append("obsfiles_posAnarrow = ['Offset_7_"+str(target)+"_PosA.txt', 'Offset_-7_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAnarrow = "'+str(target)+'_POSA_NARROW"')
        # wide slits
        if set([-7,0,7]).issubset(offsets):
            setupLines.append("obsfiles_posCwide = ['Offset_0_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            setupLines.append('target_posCwide = "'+str(target)+'_POSC_WIDE"')
            setupLines.append("obsfiles_posAwide = ['Offset_0_"+str(target)+"_PosA.txt', 'Offset_-7_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAwide = "'+str(target)+'_POSA_WIDE"')
            specphot=True

    return setupLines, obsFiles, specphot



#set noninteractive variable
if len(sys.argv) > 3:
    print "Usage: mospy AutoDriver [True|False]"
    sys.exit()

noninteractiveval=False
if len(sys.argv) == 3: 
    if str(sys.argv[2]) in ("t", "T" "true", "True"):
        noninteractiveval=True
    elif str(sys.argv[2]) in ("f", "F" "false", "False"):
        noninteractiveval=False
    else:
        print "Usage: mospy AutoDriver [True|False]"
        sys.exit()

targets_and_offsets,type = OffsetPairs()
print targets_and_offsets
print type


if 'slitmask' in targets_and_offsets:
    print "slitmask mode"
    mydriver=Driver("Driver.py","slitmask")
    mydriver.printHeader()    
    obsLines,obsFiles,specphot = SetupFiles('slitmask',targets_and_offsets['slitmask'],type)   
    mydriver.addOffsetFiles(obsFiles)
    mydriver.printMaskAndBand()
    mydriver.printnoninteractive(noninteractive=noninteractiveval)
    mydriver.printObsfiles(obsLines)
    mydriver.printFlat()
    mydriver.printWavelengthFit()
    mydriver.printBackground()
    mydriver.printRectification()
    mydriver.printExtraction()
    mydriver.CloseFile()

elif type is 'long2pos' or type is 'longslit':
    Targets = targets_and_offsets.keys()
    for target in Targets:
        print str(type)+" mode"
        obsLines,obsFiles,specphot = SetupFiles(target,targets_and_offsets[target],type)
        if type is 'longslit':
            mydriver=Driver("Longslit_"+str(target)+".py","longslit")
        elif specphot:
            mydriver=Driver("Long2pos_"+str(target)+".py","long2pos_specphot")
        else:
            mydriver=Driver("Long2pos_"+str(target)+".py","long2pos")
        mydriver.printHeader()
        mydriver.addOffsetFiles(obsFiles)
        mydriver.printMaskAndBand()
        mydriver.printnoninteractive(noninteractive=noninteractiveval)
        mydriver.printObsfiles(obsLines)
        mydriver.addLongslit()        
        mydriver.printFlat()
        mydriver.printWavelengthFit()
        mydriver.printBackground()
        mydriver.printRectification()
        mydriver.CloseFile()

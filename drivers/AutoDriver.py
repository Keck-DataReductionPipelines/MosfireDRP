

import MOSFIRE
from MOSFIRE import IO, Wavelength
import os, glob
import pyfits as pf
import time

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
        self.addLine("import os, time")
        self.addLine("import MOSFIRE")
        self.addLine("from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength")
        self.addLine("import numpy as np, pylab as pl, pyfits as pf")
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
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
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
        self.addLine("# Note: for long2pos, the row position is ignored, and the middle point of the slit is used")
        self.addLine("longslit = {'yrange': [[1062,1188],[887,1010]], 'row_position': 0, 'mode':'long2pos'}")
        
    def printWavelengthFit(self):
        if self.type is 'slitmask':
            # regular mask, with sky lines
            self.addLine("Wavelength.imcombine(obsfiles, maskname, band, waveops)")
            self.addLine("Wavelength.fit_lambda_interactively(maskname, band, obsfiles,waveops)")
            self.addLine("Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,waveops)")
            self.addLine("Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)")
            # determine waveleng name
            files = IO.list_file_to_strings(self.offsetFiles)
            self.waveName = "lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,""))
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
            self.addLine("argon = ['Ar.txt']")
            self.addLine("Wavelength.imcombine(argon, maskname, band, waveops)")
            self.addLine("Wavelength.fit_lambda_interactively(maskname, band, argon,waveops,longslit=longslit, argon=True)")
            self.addLine("Wavelength.fit_lambda(maskname, band, argon, argon,waveops,longslit=True)")
            self.addLine("Wavelength.apply_lambda_simple(maskname, band, argon, waveops, longslit=longslit, smooth=True)")            
            # determine waveleng name
            files = IO.list_file_to_strings('Ar.txt')
            self.waveName = "lambda_solution_"+str(Wavelength.filelist_to_wavename(files, self.band, self.maskName,""))
            self.addLine("Wavelength_file = '"+str(self.waveName)+"'")
        self.addLine("")
            
    def printBackground(self):
        if self.type is 'long2pos_specphot' or self.type is 'long2pos':
            for slit in ['posAnarrow','posCnarrow','posAwide','posCwide']:
                self.addLine("Background.handle_background(obsfiles_"+str(slit)+",Wavelength_file,maskname,band,waveops, target=target_"+str(slit)+")")
        if self.type is 'slitmask':
            self.addLine("Background.handle_background(obsfiles,'"+str(self.waveName)+"',maskname,band,waveops)")

        self.addLine("")

    def printRectification(self):
        if self.type is 'slitmask':
            self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles]')
            self.addLine('Rectify.handle_rectification(maskname, redfiles,"'+str(self.waveName)+'",band,obsfiles,waveops)')
        if self.type is 'long2pos' or self.type is 'long2pos_specphot':
            for slit in ['posAnarrow','posCnarrow']:
                self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles_'+str(slit)+']')            
                self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_'+str(slit)+',waveops, target=target_'+str(slit)+')')
        if self.type is 'long2pos': 
            for slit in ['posAwide','posCwide']:
                self.addLine('redfiles = ["eps_" + file + ".fits" for file in obsfiles_'+str(slit)+']')            
                self.addLine('Rectify.handle_rectification(maskname, redfiles,Wavelength_file,band,obsfiles_'+str(slit)+',waveops, target=target_'+str(slit)+')')
        self.addLine("")

    def printHeader(self):
        now = time.strftime("%c")
        self.addLine("Driver file automatically generated on "+str(now))
        self.addLine("For questions and comments, email mosfiredrp@gmail.com, submit a ticket on the ticketing system, or contact Luca Rizzi @ WMKO")
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
                # check if we have already seen this target
                tname = "_".join(off_array[1:-1])    # we are doing this for cases in which the file is Offset_-7_HIP87_7.25_PosA.txt (the _ in the file name is a problem)
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


    return targets_and_offsets

def SetupFiles(target=None, offsets=None):
    # convert numbers such as 1.0 to 1, but leaves 1.5 as 1.5
    offsets = [int(f) if f % 1 ==0 else f for f in offsets]
    setupLines = []
    obsFiles = []
    specphot = False

    # slitmask
    if target is 'slitmask':
        # files are assumed to be in pairs, and we drop the "0" value is present.
        for off in offsets:
            if off>0:
                setupLines.append("obsfiles=['Offset_"+str(off)+".txt','Offset_"+str(off*-1)+".txt']")
                obsFiles.append("Offset_"+str(off)+".txt")
                obsFiles.append("Offset_"+str(off*-1)+".txt")
    else:
        # old long 2 pos (-7,-14,-21, 7,14,21)
        # narrow slits
        if set([-7,-21,7,21]).issubset(offsets):
            setupLines.append("obsfiles_posCnarrow = ['Offset_-21_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            obsFiles.append("Offset_-21_"+str(target)+"_PosC.txt")  # we are using this to determine maskname and band
            setupLines.append('target_posCnarrow = "'+str(target)+'_posC_narrow"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posCnarrow)")
            setupLines.append("obsfiles_posAnarrow = ['Offset_7_"+str(target)+"_PosA.txt', 'Offset_21_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAnarrow = "'+str(target)+'_posA_narrow"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posAnarrow)")
        # wide slits
        if set([-7,-14,7,14]).issubset(offsets):
            setupLines.append("obsfiles_posCwide = ['Offset_-14_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            setupLines.append('target_posCwide = "'+str(target)+'_posC_wide"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posCwide)") 
            setupLines.append("obsfiles_posAwide = ['Offset_14_"+str(target)+"_PosA.txt', 'Offset_21_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAwide = "'+str(target)+'_posA_wide"')
            setupLines.append("IO.fix_long2pos_headers(obsfiles_posAwide)") 
            specphot = True
        # new long 2 pos (-7,0,7)
        # narrow slits
        if set([-7,7]).issubset(offsets) and not(set([21,21]).issubset(offsets)):
            setupLines.append("obsfiles_posCnarrow = ['Offset_7_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            obsFiles.append("Offset_7_"+str(target)+"_PosC.txt")  # we are using this to determine maskname and band
            setupLines.append('target_posCnarrow = "'+str(target)+'_posC_narrow"')
            setupLines.append("obsfiles_posAnarrow = ['Offset_7_"+str(target)+"_PosA.txt', 'Offset_-7_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAnarrow = "'+str(target)+'_posA_narrow"')
        # wide slits
        if set([-7,0,7]).issubset(offsets):
            setupLines.append("obsfiles_posCwide = ['Offset_0_"+str(target)+"_PosC.txt', 'Offset_-7_"+str(target)+"_PosC.txt']")
            setupLines.append('target_posCwide = "'+str(target)+'_posC_wide"')
            setupLines.append("obsfiles_posAwide = ['Offset_0_"+str(target)+"_PosA.txt', 'Offset_-7_"+str(target)+"_PosA.txt']")
            setupLines.append('target_posAwide = "'+str(target)+'_posA_wide"')
            specphot=True
    return setupLines, obsFiles, specphot



targets_and_offsets=OffsetPairs()

if 'slitmask' in targets_and_offsets:
    print "slitmask mode"
    mydriver=Driver("Driver_test.py","slitmask")
    mydriver.printHeader()    
    obsLines,obsFiles,specphot = SetupFiles('slitmask',targets_and_offsets['slitmask'])    
    mydriver.printObsfiles(obsLines)
    mydriver.addOffsetFiles(obsFiles)
    mydriver.printMaskAndBand()
    mydriver.printFlat()
    mydriver.printWavelengthFit()
    mydriver.printBackground()
    mydriver.printRectification()
    mydriver.CloseFile()

else:
    Targets = targets_and_offsets.keys()
    for target in Targets:
        print "long2pos mode"
        obsLines,obsFiles,specphot = SetupFiles(target,targets_and_offsets[target])
        if specphot:
            mydriver=Driver("Long2pos_"+str(target)+"_test.py","long2pos_specphot")
        else:
            mydriver=Driver("Long2pos_"+str(target)+"_test.py","long2pos")
        mydriver.printHeader()
        mydriver.addLongslit()
        mydriver.printObsfiles(obsLines)
        mydriver.addOffsetFiles(obsFiles)
        mydriver.printMaskAndBand()
        mydriver.printFlat()
        mydriver.printWavelengthFit()
        mydriver.printBackground()
        mydriver.printRectification()
        mydriver.CloseFile()



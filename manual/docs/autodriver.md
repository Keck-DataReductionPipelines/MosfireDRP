# AutoDriver

The pipeline is able to produce a driver file automatically for most cases, thus removing the need to copy one of the standard files and manually edit it.

To generate the driver file, go to the directory where your Offset files live, and where the reduction is going to happen and type:

    mospy AutoDriver

This will generate a file called Driver.py, which you should inspect before running it. Highly specialized cases such as particular combinations of sky lines and arcs might not be dealt with correctly. Note that the automatic generation of the driver file works for long2pos and longslit as well.

To handle special cases, the automatic generation of the driver file makes a number of assumptions, that might not be correct for your science case.

1. If either Ar.txt or Ne.txt or both are available, they are being used.
2. If the band is K, and FlatThermal.txtis available, it is used
3. For long2pos: if no arcs are available, only specphot in non spectrophotometric mode can be reduced and the pipeline will attempt to use the sky lines. Note that is likely to fail, as sky lines are too faint in short exposures for an accurate wavelength determination. The spectrophotometric mode contains wide slits that cannot be reduced using sky lines only.
4. In longslit, the pipeline will try to determine the size of the slit using the mask name. For example, if the maskname is LONGSLIT-3x0.7, the pipeline assumes that you have used 3 slits to generate the longslit and that they are centered around the middle line of the detector.
5. If any of the mandatory observations are missing (such as the flat fields), the pipeline will still generate a Driver.py file, but it will contain warnings about the missing files, and it will NOT run correctly.
6. If multiple observations of different stars are done in long2pos or in longslit mode, the pipeline will generate multiple driver files, one for each object. If the same object is observed multiple times during the same night, all the observations will end up in the same driver file. If you are observing a telluric standard at different times and you need to have separate spectra, you need to manually create Offset files and Driver files.


# The driver.py File

__Important Note:__ The information in this section is included for reference purposes only.  The supported method for use of the DRP is to use the `AutoDriver` step detailed above.  The example driver files discussed below and the methods for using them are not maintained.

The driver file controls all the pipeline steps, and in the drivers sub-directory, you will find a number of driver files: `Driver.py`, `K_Driver.py`, `Long2pos_driver.py`, and `Longslit_Driver.py`. The `Driver` and `K_Driver` will reduce your science data for bands Y,J, and H (this includes the sample data set). The K band requires a special approach because there are too few bright night-sky emission lines at the red end and so the `K_Driver` synthesizes arclamps and night sky lines. The `Long2pos_driver.py` handles `long2pos` and `long2pos_specphot` observations, while the `Longslit_driver.py` deals with observations of single objects using a longslit configuration.
 
The driver.py files included with the code download contains execution lines that are commented out. For this example, we will run the driver file one line at a time, but as you become familiar with the DRP process, you will develop your own driver file execution sequencing. Although in the future we hope to further automate the driver file, currently some steps require you to update the inputs with filenames created from previous steps. 

Below is a driver.py file:

    import os, time
    import MOSFIRE
    
    from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, \
         Rectify
    from MOSFIRE import Wavelength
    
    import numpy as np, pylab as pl, pyfits as pf
    
    np.seterr(all="ignore")
    
    #Update the insertmaskname with the name of the mask
    #Update S with the filter band Y,J,H,or K
    maskname = 'insertmaskname'
    band = 'S'    
    
    flatops = Options.flat
    waveops = Options.wavelength
    
    obsfiles = ['Offset_1.5.txt', 'Offset_-1.5.txt']
    
    #Flats.handle_flats('Flat.txt', maskname, band, flatops)
    #Wavelength.imcombine(obsfiles, maskname, band, waveops)
    #Wavelength.fit_lambda_interactively(maskname, band, obsfiles,
        #waveops)
    #Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,
        #waveops)

    #Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)
    #Background.handle_background(obsfiles,
        #'lambda_solution_wave_stack_H_m130429_0224-0249.fits',
        #maskname, band, waveops)

    redfiles = ["eps_" + file + ".fits" for file in obsfiles]
    #Rectify.handle_rectification(maskname, redfiles,
    #    "lambda_solution_wave_stack_H_m130429_0224-0249.fits",
    #    band, 
    #    "/scr2/npk/mosfire/2013apr29/m130429_0224.fits",
    #    waveops)
    #

To set up your driver file do the following:

1. Navigate to the desired output directory created by handle: `cd ~/Data/reducedMOSFIRE_DRP_MASK/2012sep10/H`
2. Copy the appropriate driver file: `cp ~/MosfireDRP-master/drivers/Driver.py .`  NOTE: If you are observing a K band mask you’ll want to copy the `K_driver.py` file over.
3. Edit driver.py (see bold text in driver file example)
    * Update maskname
    * Update band to be Y,J,H
    * Update the `Offset_#.txt` name. Handle creates offset files with names that are specific to the nod throw. The default driver file uses 1.5 arcsec offsets in the file name. 

In the sections that follow, we will describe the function and outputs of the commented lines found in the driver file starting with the creation of flats.

If you prefer to override the standard naming convention of the output files, you can specify

    target = “targetname” 

at the beginning of the driver file. If you do so, remember to also add target=target to both the Background and Rectify steps. Example:

    Background.handle_background(obsfiles,
        'lambda_solution_wave_stack_H_m150428_0091-0091.fits',
        maskname, band, waveops, target=target)










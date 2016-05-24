---
layout: page
title: Wavelength Calibration – Y, J, and H bands
permalink: /manual/wavelengthYJH
---

# Wavelength Calibration – Y, J, and H bands

In the shorter wavebands, when using the recommended exposure times, the wavelength calibration is performed on night sky lines. The mospy Wavelength module is responsbile for these operations. See the example driver file in section 7.

## Combine files

First step is to produce a file with which you will train your wavelength solution. Since we’re using night sky lines for training, the approach is to combine individual science exposures. This is performed by the python Wavelength.imcombine routine. For a lot of users, this will look something like in the Driver.py file:

    Wavelength.imcombine(obsfiles, maskname, band, waveops)

The first parameter is obsfiles which is a python string array indicating the list of files in the offset positions. Note that obsfiles has defaults of “Offset_1.5.txt” and “Offset_-1.5.txt” and may need to be updated as described in section 6. 

Suppose you want to exclude a file for reasons such as weather or telescope fault, simply remove the offending file from the appropriate Offset_*.txt. Likewise, you are welcome to add files in as you like, such as observations from the previous night.

Outputs of this step are:

| Filename                         | Contains                                                                    |
|----------------------------------|-----------------------------------------------------------------------------|
| `wave_stack_[band]_[range].fits` | A median-combined image of the files to be used for the wavelength solution.|

## Interactive wavelength fitting

The next step is to use the wave_stack_*.fits file and determine an initial wavelength solution for each slit. During this process, we interactively fit the lines using a gui that displays. To initiate this process, uncomment the line in the Driver.py file:

    #Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles, waveops)

And then re-execute the driver file: 

    mospy Driver.py 

when you run this step, a GUI window appears.


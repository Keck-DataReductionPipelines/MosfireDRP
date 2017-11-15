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

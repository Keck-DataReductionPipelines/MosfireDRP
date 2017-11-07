# Wavelength Calibration (K)

The night sky lines at the red end of the K-band are too faint to achieve small-fraction of a pixel RMS wavelength calibration. You will have to observe a Neon and Argon arc lamps during your afternoon calibrations. By default, the calibration script at the observatory is setup to acquire both the Ne and Argon arcs.

Because the beams emminating from the arclamp do not follow the same path as the beams coming from the sky, there will be a slight difference between the two solutions. For the afformentioned beam matching reason, the most accurate solution is the night sky lines. Thus, the code has to be clever about merging the two solutions.

The following subsections describe the additional steps that are necessary to process the arcline data and combine the arcs and night sky line wavelength solutions.

## Combine the arc line spectra

Just like the step in section 8.1 where you combined the science frames to create nightsky line spectra, we first need to combine the arcline data. The arcs are typically three files and you should see them listed in the Ne.txt and Ar.txt file lists in your K band sub directory. To combine the images simply uncomment and run:

    Wavelength.imcombine('Ne.txt', maskname, band, waveops)                                                                                                    
    Wavelength.imcombine('Ar.txt', maskname, band, waveops)                                                                                                    

## Identify arc lines using night sky solution 

Instead of having to interactively determine the wavelenth solution for the arcs like we did in section 8.2 for the night sky lines, we are going to use the solutions for the night sky lines as a first approximation for the arcs. This may usually be done because the arcs differ from the night sky lines by a fractions of pixels. You are welcome to interactively solve the neon lamp solution with the Wavelength.fit_lambda_interactively routine; however, the need to run the interactive solution method should be rare. 

To apply the solution from the night sky lines to the arcs center slit position, uncomment and run the following lines.

    Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ne.txt', neon=True)                                                              
    Wavelength.apply_interactive(maskname, band, waveops, apply=obsfiles, to='Ar.txt', argon=True)                                                             

This step, when run will produce output like:

    slitno  1 STD: 0.16 MAD: 0.06
    slitno  2 STD: 0.03 MAD: 0.02
    slitno  3 STD: 0.04 MAD: 0.04
    slitno  4 STD: 0.05 MAD: 0.01

For each slit, a new solution is generated for the neon line. The output mimics that described previously where STD is the standard deviation and MAD is the median absolute deviation in angstroms.

## Wavelength fitting for the entire slit using arcs

The next step in the wavelength fitting process is to propogate the arc solution spatially along each slit. Again this is the same process essentially as the fit for the night sky lines. This moves along each row for the slit to determine a wavelenth solution. The output files are comperable to those in step 8.3

    Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt', waveops, wavenames2='Ar.txt')                                                                    

You will note that the Ar.txt file is listed as an optional argument. If you choose not to use the Argon lamps, then you may simply remove the optional wavenames2 and execute this using only the Ne arcs.

Again, this process takes some time to complete.

## Merge the arc and sky lists

In this portion of the procedure, we merge the two lists. These commands may not be run individually. Instead any command containing the variable LROI needs to be run in one mospy driver file session in order to pass the LROI variable. In this section we determin the offsets between the region of overlap between the nightskylines and the arclines. A plot of that region is displayed. To move on you will have to close the plot. 

To execute this step you will need to uncomment the following lines in the driver file.

    LROI = [[21000, 22800]] * 1                                                                                                                                
    LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles, 'Ne.txt', LROI, waveops)        
    Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops)                                                                                          
    Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,  'Ne.txt', LROIs, waveops, neon=True)    

The merged output solution will have a filename that looks like:

    merged_lambda_coeffs_wave_stack_K_m130114_0451-0453_and_wave_stack_K_m140508_0197-0199.npy
    merged_lambda_solution_wave_stack_K_m130114_0451-0453_and_wave_stack_K_m140508_0197-0199.fits
    merged_rectified_wave_stack_K_m130114_0451-0453_and_wave_stack_K_m140508_0197-0199.fits.gz

The ouput files have the same format as those in section 8.4 and will need to be used as inputs to the Background and Rectify section below.


---
layout: page
title: Changes in Version 2016
permalink: /manual/changes
---

# 2016

## New features
* DRP is no longer dependent on IRAF/PyRAF
    * The DRP should now work with any python install which has the [required python packages](/manual/installing#Requirements)
* Updated instruction manual
* The DRP now does optimal spectral extraction (Horne 1986) and outputs a 1D spectrum

## Improvements and bug fixes

* 

# 2015A

## New features

* Reduction of long2pos and long2pos_specphot
* Reduction of longslit data
* Automatic generation of driver file
* Logging and diagnostic information on screen and on disk
* Package-style installation as a Ureka sub-package
* Support for Ureka 1.5.1

## Improvements and bug fixes

* Fix incorrect determination of the slit parameters which prevented the use of large slits
* Fix incorrect determination of the average wavelength dispersion for the long2pos mode
* Added ability of specifying the output name of the files
* Improved robustness of non-interactive wavelength solution, added possibilty of switching from interactive to non-interactive during the reduction, added k-sigma clipping of the sky or arc lines
* Fixed the problem with the interactive wavelength window not closing at the end of the fit
* Fixed the problem with the interactive fitting window showing up empty on the first fit (no need to use the x key to unzoom)
* Added procedure to fix the header of old images acquired with an outdated version of long2pos
* Disabled cosmic ray rejection for the case of less than 5 exposures
* There is no need to specify one of the observations in Rectify: Rectify will use the first of the files listed in the Offset files.
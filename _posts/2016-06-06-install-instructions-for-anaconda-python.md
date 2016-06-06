---
layout: post
title:  "Install Instructions for Anaconda Python"
date:   2016-06-06 12:00:00 -1000
---
The Ureka package has been [deprecated](http://ssb.stsci.edu/ureka/) as of April 26, 2016.  As a result, the MOSFIRE pipeline will migrate to a version which is not dependent on IRAF/PyRAF, but only on python packages.  It __should__ work with any python install which provides the required packages and versions (look for more info on the dependencies when the next DRP release comes out in late 2016).

This post is a quick instruction on how to install Anaconda and all required packages for the non-Ureka version of the MOSFIRE DRP.

* Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).
    * The pipeline should run via either python 2.7 or python 3.  We recommend python 3 as it will (in the very long term) be the only version supported by various community packages such as astropy (see this [discussion](https://github.com/astropy/astropy/issues/4454)).
* Once anaconda is installed, you can use the `conda` command line tool to get the other packages you will need.  Begin by updating conda itself by running:
    * `conda update conda`
    * If you like you can update all packages which are out of date by running: `conda update --all`.
* Install [ccdproc](http://ccdproc.readthedocs.io/en/latest/index.html) (an "astropy affiliated" package) using the astropy channel:
    * `conda install -c astropy ccdproc`

You should now have all the requirements to run the non-Ureka version of the MOSFIRE DRP.

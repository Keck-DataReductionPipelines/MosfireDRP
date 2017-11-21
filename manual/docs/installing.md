# Installation

## Requirements

The 2016 version of the pipeline no longer requires IRAF/PyRAF, so the installation should be simpler than previous versions.

The pipeline requires the following python modules:

* numpy
* astropy
* ccdproc
* scipy

## Installing python

### Using Anaconda Cloud and Conda Environments

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).

Now we will create a conda [environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) specifically for the MOSFIRE DRP.  Rather than specify everything on the command line, we will get the specification for the environment from the Anaconda Cloud service:

    conda env create joshwalawender/mospy_2016

Now we will invoke that environment:

    source activate mospy_2016

Now we will install the DRP itself.  From now on, if you want to run the DRP, first invoke the `mospy_2016` environment using `source activate mospy_2016`.


## Download and Install the DRP

Download the zip file of the released master branch from the [github page](https://github.com/Keck-DataReductionPipelines/MosfireDRP) for the DRP, or directly from [this link](https://github.com/Keck-DataReductionPipelines/MosfireDRP/archive/master.zip).

Move the zip file to a location on your computer where you want the source code to reside, then unzip the file:

    unzip MosfireDRP-master.zip

Change in to the resulting ```MosfireDRP-master/``` directory:

    cd MosfireDRP-master

Run the install program:

    python setup.py install

The executable `mospy` should now be in your path.  If you used the Anaconda based install, it will be in the Anaconda bin directory (e.g. `~/anaconda2/bin/mospy`).


## Alternate Methods of Installing Python

Note, these are no longer the recommended methods of installing the DRP as they do not gauranttee that the various package versions are compatible with the DRP.

### Using the Anaconda Distribution

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).  The pipeline currently only runs on python 2.7, so download and install that version, not the python 3.x version.

Once anaconda is installed, you can use the `conda` command line tool to get the other packages you will need.  Begin by updating conda itself by running:

    conda update conda

If you like, you can now update all packages which are out of date by running:

    conda update --all

Install [ccdproc](http://ccdproc.readthedocs.io/en/latest/index.html) (an "astropy affiliated" package) using the astropy channel:

    conda install -c astropy ccdproc

You should now have all the requirements to run the 2016 version of the MOSFIRE DRP.

Note: An Anaconda-based install has been tested on both Mac OS X 10.11.5 ("El Capitan") and on linux: Cent OS 7 (64-bit).

### Using Other Python Install Methods

The DRP support group recommends the anaconda python install and has tested the DRP using that installer, but if an appropriate version of python (e.g. python 2.7) is installed via some other package manager (e.g. apt-get, brew, yum, etc.), then you should be able to install the python package dependencies using either that package manager (if they are available via that package manager) or using `pip`.  For example:

    pip install numpy
    pip install astropy
    pip install ccdproc


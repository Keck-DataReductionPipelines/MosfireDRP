# Installation

## Requirements

The 2016 version of the pipeline no longer requires IRAF/PyRAF, so the installation should be simpler than previous versions.

The pipeline requires the following python modules:

* numpy
* astropy
* ccdproc
* scipy

## Installing Python

### Using Anaconda Cloud and Conda Environments

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).

Now we will create a conda [environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) specifically for the MOSFIRE DRP.  Rather than specify everything on the command line, we will get the specification for the environment from the Anaconda Cloud service.  There are two specifications, one for linux (tested on a CentOS 7 system) and one for macOS (tested on macOS 10.12.6).  Get the one appropriate for your system using one of the commands below:

    conda env create KeckObservatory/mospy_2016_linux

or

    conda env create KeckObservatory/mospy_2016_macOS

Now we will invoke that environment:

    source activate mospy_2016_linux

or

    source activate mospy_2016_macOS

Now we will install the DRP itself.  From now on, if you want to run the DRP, first invoke the appropriate environment using `source activate mospy_2016_linux` or `source activate mospy_2016_macOS`.


## Download and Install the DRP

Download the zip file of the released version from [GitHub](https://github.com/Keck-DataReductionPipelines/MosfireDRP/releases/download/Release2016A/MosfireDRP-2016release.zip).

Move the zip file to a location on your computer where you want the source code to reside, then unzip the file:

    unzip MosfireDRP-2016release.zip

Change in to the resulting ```MosfireDRP-2016release/``` directory:

    cd MosfireDRP-2016release

Run the install program:

    python setup.py install

The executable `mospy` should now be in your path.  If you used the Anaconda based install, it will be in the Anaconda bin directory (e.g. `~/anaconda2/bin/mospy`).


## Alternate Methods of Installing Python

Note, these are no longer the recommended methods of installing the DRP as they do not gauranttee that the various package versions are compatible with the DRP.

### Using the Anaconda Distribution

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).  The pipeline currently (2016 release) only runs on python 2.7, so download and install that version, not the python 3.x version.

To generate an environment similar to the one in the recommended anaconda cloud based install, you can use the following command:

```
conda create --no-default-packages -c astropy -n mospy_2016 python=2.7.13 astropy=2.0.2 ccdproc=1.3.0 ipython=5.4.1 numpy=1.11.2 scipy=0.18.1 PyQt=5.6.0
```

You should now have all the requirements to run the MOSFIRE DRP.  This should work on any anaconda install, even if the pre-packaged linux and macOS environments are incompatible with your machine.

### Using Other Python Install Methods

The DRP support group recommends the anaconda python install and has tested the DRP using that installer, but if an appropriate version of python (e.g. python 2.7) is installed via some other package manager (e.g. apt-get, brew, yum, etc.), then you should be able to install the python package dependencies using either that package manager (if they are available via that package manager) or using `pip`.  For example:

    pip install numpy
    pip install astropy
    pip install ccdproc


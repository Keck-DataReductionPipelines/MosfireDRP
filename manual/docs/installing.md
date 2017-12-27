# Installation

## Requirements

The pipeline requires the following python modules:

* numpy
* astropy
* ccdproc
* scipy

## Installing Python

### Using Anaconda Cloud and Conda Environments

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).

Now we will create a conda [environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) specifically for the MOSFIRE DRP.  Rather than specify everything on the command line, we will get the specification for the environment from the Anaconda Cloud service.  There are two specifications, one for linux (tested on a CentOS 7 system) and one for macOS (tested on macOS 10.12.6).  Get the one appropriate for your system using one of the commands below:

    conda env create KeckObservatory/mospy_2017_linux

or

    conda env create KeckObservatory/mospy_2017_macos

Now we will invoke that environment:

    source activate mospy_2017_linux

or

    source activate mospy_2017_macOS

Now we will install the DRP itself.  From now on, if you want to run the DRP, first invoke the appropriate environment using `source activate mospy_2017_linux` or `source activate mospy_2017_macos`.


## Download and Install the DRP

Download the zip file of the released version from [GitHub](https://github.com/Keck-DataReductionPipelines/MosfireDRP/releases/download/Release2017/MosfireDRP-2017release.zip).

Move the zip file to a location on your computer where you want the source code to reside, then unzip the file:

    unzip MosfireDRP-2017release.zip

Change in to the resulting ```MosfireDRP-2017release/``` directory:

    cd MosfireDRP-2017release

Run the install program:

    python setup.py install

The executable `mospy` should now be in your path.  If you used the Anaconda based install, it will be in the Anaconda bin directory (e.g. `~/anaconda/envs/mospy_2017_macos/bin/mospy`).


## Alternate Methods of Installing Python

Note, these are no longer the recommended methods of installing the DRP as they do not gauranttee that the various package versions are compatible with the DRP.

### Using the Anaconda Distribution

Install Anaconda as per the instructions on the [Anaconda web site](https://www.continuum.io/downloads).  The pipeline currently (2016 release) only runs on python 2.7, so download and install that version, not the python 3.x version.

To generate an environment similar to the one in the recommended anaconda cloud based install, you can use the following command:

```
conda create --no-default-packages -c astropy -n mospy_2017_macos python=3.6.3 astropy=2.0.3 ccdproc=1.3.0 ipython=6.2.1 numpy=1.13.3 scipy=1.0.0 PyQt=5.6.0
```

You should now have all the requirements to run the MOSFIRE DRP.  This should work on any anaconda install, even if the pre-packaged linux and macOS environments are incompatible with your machine.

### Using Other Python Install Methods

The DRP support group recommends the anaconda python install and has tested the DRP using that installer, but if an appropriate version of python is installed via some other package manager (e.g. apt-get, brew, yum, etc.), then you should be able to install the python package dependencies using either that package manager (if they are available via that package manager) or using `pip`.  For example:

    pip install numpy
    pip install astropy
    pip install ccdproc


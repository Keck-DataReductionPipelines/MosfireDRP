#!/usr/stsci/pyssg/Python-2.7/bin/python
from distutils.core import setup

setup(name='MOSFIRE',
    version='1.0',
    description='Python Distribution Utilities',
    author='Nick Konidaris',
    author_email='npk@astro.caltech.edu',
    url='http://mosfire.googlecode.com',
    py_modules=['Background', 'Bspline', 'Combine', 'CSU', 'Detector',
    'Filters', 'Fit', 'Flats', 'IO', 'Longslit', 'Options', 'Rectify', 'Util',
    'Wavelength','MosfireDrpLog'])



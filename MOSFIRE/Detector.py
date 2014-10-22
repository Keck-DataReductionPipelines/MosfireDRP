'''
MOSFIRE Detector Code

File contains values associated with the Hawaii-2RG detector.

pixelsize of 18 micron is provided by Teledyne.

'''

import numpy as np

mm = 1.0

pixelsize = 0.018 * mm
npix = (2048, 2048)

gain = 2.15 # From MOSFIRE pre ship review pg. 125
RN = 21.0



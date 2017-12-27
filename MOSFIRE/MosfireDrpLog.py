import logging, sys

LOG_FILENAME = "mosfire_DRP.log"


formatter = logging.Formatter('%(asctime)s - %(module)12s.%(funcName)20s - %(levelname)s: %(message)s')
# set up logging to STDOUT for all levels DEBUG and higher
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)

# set up logging to a file for all levels DEBUG and higher
fh = logging.FileHandler(LOG_FILENAME)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)

# create Logger object
mylogger = logging.getLogger('MyLogger')
mylogger.setLevel(logging.DEBUG)
mylogger.addHandler(sh)    # enabled: stdout
mylogger.addHandler(fh)    # enabled: file


# create shortcut functions
debug = mylogger.debug
info = mylogger.info
warning = mylogger.warning
error = mylogger.error
critical = mylogger.critical


#logger.info ("log started")

# Add log entries for versions of numpy, matplotlib, astropy, ccdproc
info(sys.version)
info('python version = {}.{}.{}'.format(sys.version_info.major,
                                        sys.version_info.minor,
                                        sys.version_info.micro))
import numpy as np
info('numpy version = {}'.format(np.__version__))
import matplotlib
info('matplotlib version = {}'.format(matplotlib.__version__))
import astropy
info('astropy version = {}'.format(astropy.__version__))
import ccdproc
info('ccdproc version = {}'.format(ccdproc.__version__))

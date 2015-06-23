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

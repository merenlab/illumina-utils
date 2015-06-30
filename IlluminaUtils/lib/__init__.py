import os
import pkg_resources

try:
    __version__ = pkg_resources.require("illumina-utils")[0].version
except:
    # maybe it is not installed but being run from the codebase dir?
    try:
        __version__ = open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/../../VERSION').read().strip()
    except:
        __version__ = 'unknown'

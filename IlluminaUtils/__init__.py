import os
import sys
import pkg_resources

# Make sure the Python environment hasn't changed since the installation (happens more often than you'd think
# on systems working with multiple Python installations that are managed through modules):
try:
    if sys.version_info < (2, 7, 5):
        v =  '.'.join([str(x) for x in sys.version_info[0:3]])
        sys.stderr.write("Your active Python version is '%s'. Anything less than '2.7.5' will not do it for illumina-utils :/\n" % v)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(illumina-utils failed to learn about your Python version, but it will pretend as if nothing happened)\n\n")


def set_version():
    try:
        __version__ = pkg_resources.require("illumina-utils")[0].version
    except:
        # maybe it is not installed but being run from the codebase dir?
        try:
            __version__ = open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/../VERSION').read().strip()
        except:
            __version__ = 'unknown'

    return __version__


__version__ = set_version()

def print_version():
    print "Illumina-utils v%s" % __version__


if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()

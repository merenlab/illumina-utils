import sys

illumina_utils_version = "2.11"

# Make sure the Python environment hasn't changed since the installation (happens more often than you'd think
# on systems working with multiple Python installations that are managed through modules):
try:
    if sys.version_info.major != 3:
        sys.stderr.write("Your active Python major version ('%d') is not compatible with what illumina-utils expects :/ We recently switched to Python 3.\n" % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(illumina-utils failed to learn about your Python version, but it will pretend as if nothing happened)\n\n")

__version__ = illumina_utils_version

def print_version():
    print("v%s" % illumina_utils_version)


if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()

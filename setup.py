import os
import uuid
import glob
from setuptools import setup, find_packages

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

init_py_path = os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/IlluminaUtils/__init__.py'
version_string = [l.strip() for l in open(init_py_path).readlines() if l.strip().startswith('illumina_utils_version')][0]
illumina_utils_version = version_string.split('=')[1].strip().strip("'").strip('"')

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

install_reqs = parse_requirements('requirements.txt', session=uuid.uuid1())
try:
    reqs = [str(ir.requirement) for ir in install_reqs]
except AttributeError:
    reqs = [str(ir.req) for ir in install_reqs]

setup(
    name = "illumina-utils",
    version = illumina_utils_version,
    description = "A library and collection of scripts to work with Illumina paired-end data (for CASAVA 1.8+).",
    author = u"Illumina Utils Authors",
    author_email = "a.murat.eren@gmail.com",
    license = "GPLv3+",
    url = "https://github.com/meren/illumina-utils",
    packages = find_packages(),

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
    ],

    scripts = [script for script in glob.glob('scripts/*') if not script.endswith('-OBSOLETE')],

    include_package_data = True,
    package_data={'': ['examples']},

    install_requires=reqs,
)

import os
import uuid
import glob
from setuptools import setup, find_packages

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

install_reqs = parse_requirements('requirements.txt', session=uuid.uuid1())
reqs = [str(ir.req) for ir in install_reqs]

setup(
    name = "illumina-utils",
    version = open('VERSION').read().strip(),
    description = "A library and collection of scripts to work with Illumina paired-end data (for CASAVA 1.8+).",
    author = u"A. Murat Eren",
    author_email = "meren@mbl.edu",
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

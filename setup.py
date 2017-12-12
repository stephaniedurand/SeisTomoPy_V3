"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from distutils.ccompiler import CCompiler
from distutils.errors import DistutilsExecError, CompileError
from distutils.unixccompiler import UnixCCompiler
from setuptools.extension import Extension
import urllib
import zipfile

import inspect
import os
from subprocess import Popen, PIPE
import sys


here = path.abspath(path.dirname(__file__))

DIR_OR = os.path.dirname(os.path.abspath(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup_config = dict(
    name='SeisTomoPy',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.3.0',

    description='Python tools for using global tomographic models',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/stephaniedurand/SeisTomoPy',

    # Author details
    author='S. Durand, R. Abreu & C. Thomas',
    author_email='durand@uni-muenster.de',

    # Choose your license
    license='GPLv2',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        #'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Operating System :: Unix",
        "Operating System :: MacOS",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    packages=find_packages(),
    package_data={
        DIR_OR + "/SeisTomoPy":
        [os.path.join("gui","qt_window.ui")]},
            # ["gui/qt_window.ui"]},
    # packages = ['SeisTomoPy'],
    # entry_points = {
    #     'console_scripts': ['seistomopy = seistomopy.__main__:main']
        # },

    # What does your project relate to?
    keywords='global tomography',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    #packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy >= 1.11.0','obspy >= 1.0.1','scipy >= 0.17.0','pyqtgraph >= 0.9.10','basemap >= 1.0.7','matplotlib >= 2.0.2'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    #extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)

if __name__ == "__main__":
    setup(**setup_config)

    DIR_OR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(DIR_OR + '/SeisTomoPy/fortran_files/src')
    os.system('make clean')
    os.system('make all')

    os.chdir('../')
    urllib.urlretrieve("http://earth.uni-muenster.de/~durand/models.zip","models.zip")
    zip_ref = zipfile.ZipFile('models.zip', 'r')
    zip_ref.extractall()
    zip_ref.close()
    os.remove('models.zip')
    os.chdir(DIR_OR)
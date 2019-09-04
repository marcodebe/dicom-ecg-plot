#!/usr/bin/env python

# Use setuptools if we can
try:
    from setuptools.core import setup
except ImportError:
    from distutils.core import setup
from ecg import __version__

setup(name='dicomecg_convert',
      version=__version__,
      description='Dicom ECG Conversion',
      long_description=open('README.md').read(),
      author='Marco De Benedetto',
      author_email='debe@galliera.it',
      url='https://github.com/marcodebe/dicomecg_convert',
      download_url='https://github.com/marcodebe/dicomecg_convert/downloads',
      license='MIT',
      packages=['ecg', ],
      )

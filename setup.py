#!/usr/bin/env python
import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

with open("VERSION", "r") as fh:
    version = fh.read().strip('\n')


setuptools.setup(
        name='dicom-ecg-plot',
        version=version,
        description='Plot Dicom ECG Waveforms',
        long_description=long_description,
        long_description_content_type="text/markdown",
        author='Marco De Benedetto',
        author_email='debe@galliera.it',
        url='https://github.com/marcodebe/dicomecg_convert',
        packages=setuptools.find_packages(),
        scripts=['dicom-ecg-plot'],
        install_requires=[
            'pydicom>=1.0.1',
            'numpy',
            'matplotlib',
            'scipy',
            'docopt',
            'requests',
            ],
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'Intended Audience :: Healthcare Industry',
            ],
        )

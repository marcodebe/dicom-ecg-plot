[logo]: https://raw.github.com/marcodebe/dicomecg_convert/master/images/logo.png
![ECG Dicom Convert][logo]

Dicom ECG Conversion
====================
Convert a Dicom ECG (waveform) file to PDF

Install
-------
The python library dependencies are:
* dicom
* numpy
* matplotlib
* scipy

You can install the corresponding packages from your distribution or in a virtualenv.

### Without virtualenv
```bash
apt-get install python-matplotlib python-dicom python-scipy python-numpy
git clone git@github.com:marcodebe/dicomecg_convert.git
```

### Inside a virtualenv

Installing the dependencies inside the virtualenv could be long and not smooth.
I had to install system libraries and a fortran compiler.

```bash
git clone git@github.com:marcodebe/dicomecg_convert.git
virtualenv dicomecg_convert
. dicomecg_convert/bin/activate
sudo apt-get install libblas-dev
sudo apt-get install liblapack-dev 
sudo apt-get install gfortran
pip install pydicom
pip install numpy
pip install matplotlib
pip install cython
pip install git+http://github.com/scipy/scipy/
```

Usage
-----
```bash
python ecgconvert.py sample_files/anonymous_ecg.dcm
```

The output is in out.pdf

The sample file is a 12-lead ECG produced by Mortara equipment.

The signals are filtered using a bandpass (0.05-40 Hz) butterworth filter of order 1.

Work in progress, we need:
 * save in more formats (it's about to change the extension of output file name)
 * print textual info (patient name, etc.) in the header
 * different layouts
 * exception handling
 * ...

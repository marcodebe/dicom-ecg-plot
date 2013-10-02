[logo]: https://raw.github.com/marcodebe/dicomecg_convert/master/images/logo.png
![ECG Dicom Convert][logo]

Dicom ECG Conversion
====================
Convert Dicom ECG (waveform) to PDF, PNG, etc.

Usage
-----
```bash
python convert.py <inputfile> [--layout=LAYOUT] [--output=FILE|--format=FMT]
python convert.py <stu> <ser> <obj> [--layout=LAYOUT] [--output=FILE|--format=FMT]
python convert.py --help
```

Input can be a file or the triplet studyUID, seriesUID, objectUID. In the latter
case dicom file is downloaded via [WADO](http://medical.nema.org/Dicom/2011/11_18pu.pdf).

The ouput format is deduced from the extension of the filename.

If output file is not given ```--format``` must be defined.

Supported output formats are: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.

The signals are filtered using a bandpass (0.05-40 Hz) butterworth filter of order 1.

```LAYOUT``` can be one of: 3X4\_1 (that is 3 rows for 4 colums plus 1 row), 3X4, 12X1 (default: 3X4_1)

New layouts can be defined adding the corresponding matrix in LAYOUT dictionary in ```config.py```.

Work in progress, we need:
 * print textual info (patient and wave info)
 * different layouts
 * exception handling
 * ...

The sample file is a 12-lead ECG anonymized  dicom file produced
by Mortara equipment and so are all the ECG files I have to test the program.

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
sudo apt-get install python-matplotlib python-dicom python-scipy python-numpy
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

Structure of DICOM Waveform
---------------------------
```
       Series
          |
 contains | (1,n)
          |
          |
          v
       Waveform         * Waveform Attributes
          |                 Time of Acquisition
 contains | (1,n)           Acquisition Context
          |                 Annotation
          |
          v
    Multiplex Group     * Multiplex Group Attributes
          |                 Number of Channels
 contains | (1,n)           Sampling Frequenciy
          |                 Timing
          |
          v
       Channel          * Channel Definition Attributes
          |                 Channel Source
          |                   Metric
          |                   Anatomic Location(s)
          |                   Function
 contains | (1,n)             Tecnique
          |                 Channel Sensitivity
          |                 Baseline
          |                 Skew
          |                 Filter Characteristics
          |
          v
        Sample

```

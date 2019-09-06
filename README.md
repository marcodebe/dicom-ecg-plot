[logo]: https://raw.github.com/marcodebe/dicomecg_convert/master/images/logo.png
![ECG Dicom Convert][logo]

# Dicom ECG plot
A python tool to plot Dicom ECG.

The DICOM file can also be specified as `studyUID seriesUID objectUID` and 
retrieved from your WADO server.

Github repository: [here](https://github.com/marcodebe/dicomecg_convert)

**THE PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT WITHOUT ANY WARRANTY OF ANY KIND.**


## Online demo
**[demo site](https://ecg.galliera.it)** 
You can convert your own DICOM files or use preloaded sample files from different modality models.

## Install
```bash
python3 -m venv ecg
. ecg/bin/activate
pip install dicom-ecg-plot
```

## Usage of `dicom-ecg-plot` tool
```bash
dicom-ecg-plot <inputfile> [--layout=LAYOUT] [--output=FILE|--format=FMT] --minor-grid
dicom-ecg-plot <stu> <ser> <obj> [--layout=LAYOUT] [--output=FILE|--format=FMT] --minor-grid
dicom-ecg-plot --help
```
Examples:
```bash
dicom-ecg-plot anonymous_ecg.dcm -o anonymous_ecg.pdf
dicom-ecg-plot anonymous_ecg.dcm --layout 6x2 --output anonymous_ecg.png
dicom-ecg-plot anonymous_ecg.dcm --format svg > anonymous_ecg.svg
```

The input can be a (dicom ecg) file or the triplet `studyUID, seriesUID,
objectUID`. In the latter case dicom file is downloaded via
[WADO](http://medical.nema.org/Dicom/2011/11_18pu.pdf).

If `--output` is given the ouput format is deduced from the extension of the `FILE`.
If the output file is not given `--format` must be defined.
Supported output formats are: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.

By default the 5mm grid is drawn, `--minor-grid` add the minor grid (1mm).

The signals are filtered using a lowpass (40 Hz)
[butterworth filter](http://en.wikipedia.org/wiki/Butterworth_filter) 
of order 2.

`LAYOUT` can be one of: 3x4\_1 (that is 3 rows for 4 columns plus 1 row), 3x4, 6x2, 12x1 (default: 3x4_1).
New layouts can be defined adding the corresponding matrix in LAYOUT dictionary in `config.py`.



## References
 * http://medical.nema.org/Dicom/supps/sup30_lb.pdf
 * http://dicomlookup.com/html/03_03PU.html#LinkTarget_229354
 * http://libir.tmu.edu.tw/bitstream/987654321/21661/1/B09.pdf
 * [Mortara ECG Conformance Statement](http://www.mortara.com/fileadmin/user_upload/global/Products/Healthcare/DICOM/ELI%20Electrocardiographs%20DICOM%20Conformance%20Statement.pdf)

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
dicom-ecg-plot <inputfile> [--layout=LAYOUT] [--paper=PAPER] [--output=FILE|--format=FMT] [--minor-grid] [--interpretation]
dicom-ecg-plot <stu> <ser> <obj> [--layout=LAYOUT] [--paper=PAPER] [--output=FILE|--format=FMT] [--minor-grid] [--interpretation]
dicom-ecg-plot --help
```

Options:
- `--layout=LAYOUT` — lead arrangement (default: `3x4_1`)
- `--paper=PAPER` — paper format: `a4` or `letter` (default: `a4`)
- `--output=FILE` — output file; format deduced from extension
- `--format=FMT` — explicit output format (used when no output file is given)
- `--minor-grid` — draw 1mm minor grid in addition to the default 5mm grid
- `--interpretation` — show automated ECG interpretation text (if present in DICOM)

Examples:
```bash
dicom-ecg-plot anonymous_ecg.dcm -o anonymous_ecg.pdf
dicom-ecg-plot anonymous_ecg.dcm --layout 6x2 --output anonymous_ecg.png
dicom-ecg-plot anonymous_ecg.dcm --format svg > anonymous_ecg.svg
dicom-ecg-plot anonymous_ecg.dcm --paper letter --minor-grid -o anonymous_ecg.pdf
dicom-ecg-plot anonymous_ecg.dcm --interpretation -o anonymous_ecg.pdf

# WADO retrieval
dicom-ecg-plot <studyUID> <seriesUID> <objectUID> --format pdf > anonymous_ecg.pdf
```

The input can be a DICOM ECG file or the triplet `studyUID seriesUID objectUID`.
In the latter case the DICOM file is downloaded via WADO from the server configured
in `ecgconfig.py`.

If `--output` is given the output format is deduced from the extension of `FILE`.
If no output file is given, `--format` must be specified and the result is written to stdout.
Supported output formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.

By default the 5mm grid is drawn; `--minor-grid` adds the 1mm minor grid.

The signals are filtered using a lowpass (40 Hz)
[Butterworth filter](https://en.wikipedia.org/wiki/Butterworth_filter) of order 2.

`LAYOUT` can be one of: `3x4_1` (3 rows × 4 columns + 1 rhythm strip), `3x4`, `6x2`, `12x1` (default: `3x4_1`).
Custom layouts can be defined in `ecgconfig.py` by adding entries to the `LAYOUT` dictionary.



## References
 * http://medical.nema.org/Dicom/supps/sup30_lb.pdf
 * http://dicomlookup.com/html/03_03PU.html#LinkTarget_229354
 * http://libir.tmu.edu.tw/bitstream/987654321/21661/1/B09.pdf
 * [Mortara ECG Conformance Statement](http://www.mortara.com/fileadmin/user_upload/global/Products/Healthcare/DICOM/ELI%20Electrocardiographs%20DICOM%20Conformance%20Statement.pdf)

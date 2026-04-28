[logo]: https://raw.github.com/marcodebe/dicomecg_convert/master/images/logo.png
![ECG Dicom Convert][logo]

# dicom-ecg-plot

A Python tool to render 12-lead ECG waveforms from DICOM files.

The input can be a local `.dcm` file or a WADO triplet (`studyUID seriesUID objectUID`);
in the latter case the DICOM object is downloaded from the server configured in `ecgconfig.py`.

**THE PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT WITHOUT ANY WARRANTY OF ANY KIND.**

---

## About

ECG devices and PACS systems store electrocardiograms in the DICOM Waveform format (Supplement 30),
but most viewers either require proprietary software or render the tracing as a static image embedded
in the DICOM file itself. `dicom-ecg-plot` reads the raw waveform samples directly, applies standard
signal processing, and produces a publication-quality plot you can embed in reports, share as a PDF,
or pipe into automated workflows.

It is designed for use in clinical informatics and research environments where ECGs need to be
extracted from a PACS/RIS and converted to an open, portable format without depending on the
originating device's software.

## Features

- **Reads raw DICOM waveform data** — decodes waveform samples, lead labels, and metadata directly
  from DICOM tags, independent of the device manufacturer
- **Multiple lead layouts** — `3x4_1` (standard 12-lead + rhythm strip), `3x4`, `6x2`, `12x1`;
  custom layouts configurable in `ecgconfig.py`
- **WADO retrieval** — fetches DICOM objects from a PACS via WADO-URI using a
  `studyUID / seriesUID / objectUID` triplet
- **Signal filtering** — 2nd-order 40 Hz lowpass Butterworth filter to remove high-frequency noise
- **Standard ECG grid** — 5 mm major grid with optional 1 mm minor grid
- **Automated interpretation** — optionally renders the machine-generated interpretation text stored
  in the DICOM file
- **Flexible output** — writes PDF, PNG, SVG, TIFF and other Matplotlib-supported formats; streams
  to stdout for pipeline use
- **Configurable branding** — institution name can be read from the DICOM tag or overridden in config

---

## Online demo

**[ecg.galliera.it](https://ecg.galliera.it)** — upload your own DICOM files or explore preloaded samples from different ECG device models.

---

## Install

```bash
python3 -m venv ecg
. ecg/bin/activate
pip install dicom-ecg-plot
```

---

## Usage

```
dicom-ecg-plot <inputfile> [options]
dicom-ecg-plot <studyUID> <seriesUID> <objectUID> [options]
dicom-ecg-plot --help
```

### Options

| Option | Description | Default |
|---|---|---|
| `--layout=LAYOUT` | Lead arrangement | `3x4_1` |
| `--paper=PAPER` | Paper size: `a4` or `letter` | `a4` |
| `--output=FILE` | Output file; format inferred from extension | — |
| `--format=FMT` | Explicit output format (used when no output file is given) | — |
| `--minor-grid` | Draw 1 mm minor grid in addition to the default 5 mm grid | off |
| `--interpretation` | Include automated ECG interpretation text, if present in DICOM | off |

### Output formats

`eps`, `jpeg`, `jpg`, `pdf`, `pgf`, `png`, `ps`, `raw`, `rgba`, `svg`, `svgz`, `tif`, `tiff`

When `--output` is given the format is inferred from the file extension.  
When no output file is given, `--format` must be specified and the result is written to stdout.

### Layouts

| Layout | Description |
|---|---|
| `3x4_1` | 3 rows × 4 columns + 1 rhythm strip (default) |
| `3x4` | 3 rows × 4 columns |
| `6x2` | 6 rows × 2 columns |
| `12x1` | 12 rows × 1 column |

Custom layouts can be added to `ecgconfig.py` by extending the `LAYOUT` dictionary.

---

## Python API

The `ECG` class can be used directly in Python scripts. It supports the context manager protocol for automatic resource cleanup.

```python
from ecg import ECG

# Convert a local DICOM file to PNG
with ECG('anonymous_ecg.dcm') as ecg:
    ecg.draw('3x4_1')
    ecg.save('anonymous_ecg.png')
```

```python
# Render to bytes (e.g. for a web response or in-memory processing)
with ECG('anonymous_ecg.dcm') as ecg:
    ecg.draw('3x4_1')
    png_bytes = ecg.save(outformat='png')
```

```python
# Custom layout, letter paper, minor grid
with ECG('anonymous_ecg.dcm', paper='letter') as ecg:
    ecg.draw('6x2', mm_mv=10, minor_axis=True)
    ecg.save('anonymous_ecg.pdf')
```

```python
# WADO retrieval from a PACS
with ECG({'stu': '<studyUID>', 'ser': '<seriesUID>', 'obj': '<objectUID>'}) as ecg:
    ecg.draw('3x4_1', interpretation=True)
    ecg.save('anonymous_ecg.pdf')
```

`ECG.draw()` parameters:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `layoutid` | `str` | — | Layout key: `3x4_1`, `3x4`, `6x2`, `12x1` |
| `mm_mv` | `float` | `10.0` | Amplitude scale in mm/mV |
| `minor_axis` | `bool` | `False` | Draw 1 mm minor grid |
| `interpretation` | `bool` | `False` | Show automated interpretation text |

---

## Examples

```bash
# Local file → PDF
dicom-ecg-plot anonymous_ecg.dcm --output anonymous_ecg.pdf

# Custom layout → PNG
dicom-ecg-plot anonymous_ecg.dcm --layout 6x2 --output anonymous_ecg.png

# Pipe to stdout
dicom-ecg-plot anonymous_ecg.dcm --format svg > anonymous_ecg.svg

# Letter paper with minor grid
dicom-ecg-plot anonymous_ecg.dcm --paper letter --minor-grid --output anonymous_ecg.pdf

# Include automated interpretation
dicom-ecg-plot anonymous_ecg.dcm --interpretation --output anonymous_ecg.pdf

# WADO retrieval
dicom-ecg-plot <studyUID> <seriesUID> <objectUID> --format pdf > anonymous_ecg.pdf
```

---

## Configuration

Copy or edit `ecgconfig.py` to customise the tool:

```python
# WADO server URL
WADOSERVER = "http://your-wado-server/"

# Institution name displayed on the plot (None = read from DICOM tag InstitutionName)
INSTITUTION = None

# Custom layouts: list of rows, each row is a list of lead indices (0–11)
LAYOUT = { ... }
```

---

## Signal processing

Signals are filtered with a 2nd-order 40 Hz lowpass
[Butterworth filter](https://en.wikipedia.org/wiki/Butterworth_filter).
The grid follows standard ECG conventions: 5 mm major grid, optional 1 mm minor grid.

---

## References

- [DICOM Supplement 30 — Waveforms](http://medical.nema.org/Dicom/supps/sup30_lb.pdf)
- [DICOM lookup — Waveform modules](http://dicomlookup.com/html/03_03PU.html#LinkTarget_229354)
- [ECG signal encoding in DICOM](http://libir.tmu.edu.tw/bitstream/987654321/21661/1/B09.pdf)
- [Mortara ECG DICOM Conformance Statement](http://www.mortara.com/fileadmin/user_upload/global/Products/Healthcare/DICOM/ELI%20Electrocardiographs%20DICOM%20Conformance%20Statement.pdf)

---

## License

MIT — see [LICENSE](LICENSE).

## Authors

Marco De Benedetto, Simone Ferretti, Francesco Formisano

dcm_ecg2pdf
===========

Convert a Dicom ECG (waveform) file to PDF

Usage:
```bash
python ecgconvert.py sample_files/anonymous_ecg.dcm
```

The output is in out.pdf

The sample file is a 12-lead ECG produced by Mortara equipement.

The signals are filtered using a bandpass (0.05-40 Hz) butterworth filter of order 1.

Work in progress, we need:
 * print textual info (patient name, etc.) in the header
 * different layouts
 * exception handling
 * ...

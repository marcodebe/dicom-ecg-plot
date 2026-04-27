# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

---

## [1.4.1] — 2026-04-28

### Fixed
- `processor.py`: `ChannelSensitivityCorrectionFactor` tag now accessed with a safe fallback (`1.0`) instead of a bare attribute access that would crash on non-standard DICOM files
- `processor.py`: `ChannelSensitivityUnitsSequence` now handled with a fallback to `'uV'` when the sequence is absent
- `plotter.py`: `ChannelSourceSequence` now accessed safely; missing sequence falls back to a generic label (`Ch{n}`)
- `plotter.py`: `PatientName` tag now read via `dicom.get()` to avoid `AttributeError` on files where the tag is absent
- `plotter.py`: invalid layout ID now raises a clear `ValueError` instead of a cryptic `KeyError`
- `reader.py`: WADO requests now include a 30 s timeout and call `raise_for_status()`, with a readable error message on failure
- `ecg.py`: added public `close()` method and context manager support (`with ECG(...) as ecg:`)
- `dicom-ecg-plot`: removed Python 2 stdout fallback (project is Python 3 only)
- `i18n.py`: locale path is now absolute and relative to the module file, independent of the working directory

### Documentation
- Rewrote README in English: added About section, Features list, options table, layout table, configuration section, and authors

---

## [1.4.0] — 2026-03-14

Major refactor: the monolithic `ECG` class has been split into focused modules and the `scipy` dependency has been removed.

### Added
- `ecg/reader.py` — DICOM reading and metadata extraction
- `ecg/processor.py` — waveform signal extraction and scaling
- `ecg/plotter.py` — Matplotlib rendering logic
- `ecg/dsp.py` — pure-Python DSP (Butterworth lowpass filter), replacing `scipy`
- `--paper letter` option in the CLI
- Graceful error handling in the CLI (invalid files, missing tags)

### Changed
- Lead labels now resolved via DICOM `CodeValue` lookup instead of text manipulation
- DICOM date and patient age parsing made more robust (handles missing or malformed tags)
- All internal comments and docstrings translated from Italian to English

### Removed
- `scipy` dependency

---

## [1.3.5] — 2025-07-25

### Fixed
- Plot area height calculation was incorrect in some layouts
- 6×2 layout was rendering lead V1 twice instead of the correct lead pair

---

## [1.3.3] — 2022-01-28

### Fixed
- Unknown or missing `VRate` tag no longer raises an unhandled exception

---

## [1.3.2] — 2021-04-13

### Fixed
- Acquisition date/time missing or malformed in DICOM no longer causes a crash

---

## [1.3.1] — 2020-11-11

### Fixed
- `papertype` argument removed from `savefig` call (deprecated in Matplotlib)
- `MaxNLocator` / `AutoMinorLocator` `numticks` argument coerced to `int`
- WADO response decoded as binary string correctly
- Broken pip package from previous release

---

## [1.3.0] — 2020-05-11

### Added
- Ventricular rate read from the `VRate` DICOM tag

---

## [1.2.0] — 2019-09-06

### Changed
- Project renamed from `dicomecg_convert` to `dicom-ecg-plot` and published on PyPI
- Migrated to `pydicom` (dropped the legacy `dicom` import)
- Binary output to stdout handled correctly on Python 3

---

## [1.1.0] — 2018-02-26

### Changed
- Full Python 2 → Python 3 migration

---

## [1.0.0] — 2016-07-18

### Fixed
- WADO `contentType` set to `application/dicom`

---

## Earlier history (2013–2015)

- **2015** — Strip microseconds from acquisition date; handle invalid DICOM gracefully
- **2014** — Optional automated interpretation text (`--interpretation`); i18n/l10n support
- **2013** — Initial release: DICOM waveform decoding, 3×4+1 / 3×4 / 6×2 layouts, WADO retrieval, Butterworth lowpass filter, configurable institution name, multiple output formats

[Unreleased]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/compare/v1.4.1...HEAD
[1.4.1]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/compare/v1.4.0...v1.4.1
[1.4.0]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/compare/v1.3.5...v1.4.0
[1.3.5]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/compare/v1.3.3...v1.3.5
[1.3.3]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/compare/v1.3.0...v1.3.3
[1.3.2]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.3.2
[1.3.1]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.3.1
[1.3.0]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.3.0
[1.2.0]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.2.0
[1.1.0]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.1.0
[1.0.0]: https://gitlab.galliera.it/debe/dicom-ecg-plot/-/tags/v1.0.0

# -*- coding: utf-8 -*-
"""DICOM ECG reader: loads and validates a DICOM ECG waveform file."""

from __future__ import annotations

import io
from typing import BinaryIO

import requests
import pydicom

try:
    from ecgconfig import WADOSERVER
except ImportError:
    from .constants import DEFAULT_WADOSERVER as WADOSERVER


class ECGReadFileError(pydicom.filereader.InvalidDicomError):
    pass


class DicomECGReader:
    """Loads a DICOM ECG file and exposes its waveform data.

    Accepted source types:
    - str or file-like: path or buffer passed directly to pydicom.dcmread()
    - dict with keys 'stu', 'ser', 'obj': fetched from a WADO server
    """

    def __init__(self, source: str | BinaryIO | dict[str, str]) -> None:
        self.dicom = self._load(source)
        self._validate()
        self._extract_waveform()

    # ------------------------------------------------------------------
    # Loading

    def _load(self, source):
        if isinstance(source, dict):
            inputdata = self._fetch_wado(source)
        elif isinstance(source, str) or hasattr(source, 'getvalue'):
            inputdata = source
        else:
            raise ECGReadFileError(
                "'source' must be a path/to/file string, a buffer, "
                "or a dict with keys 'stu', 'ser', 'obj'."
            )
        try:
            return pydicom.dcmread(inputdata)
        except pydicom.filereader.InvalidDicomError:
            source_name = source if isinstance(source, str) else "input"
            raise ECGReadFileError(
                f"'{source_name}' is not a valid DICOM file."
            )

    def _fetch_wado(self, source):
        """Retrieve a DICOM object from a WADO server."""
        if set(source.keys()) != {'stu', 'ser', 'obj'}:
            raise ECGReadFileError(
                "WADO source dict must have exactly keys 'stu', 'ser', 'obj'."
            )
        payload = {
            'requestType': 'WADO',
            'contentType': 'application/dicom',
            'studyUID': source['stu'],
            'seriesUID': source['ser'],
            'objectUID': source['obj'],
        }
        headers = {'content-type': 'application/json'}
        resp = requests.get(WADOSERVER, params=payload, headers=headers)
        return io.BytesIO(resp.content)

    # ------------------------------------------------------------------
    # Validation

    def _validate(self):
        """Validate that the loaded DICOM contains a valid ECG waveform."""
        if not hasattr(self.dicom, 'WaveformSequence') or not self.dicom.WaveformSequence:
            raise ECGReadFileError("DICOM file has no WaveformSequence.")

        seq = self.dicom.WaveformSequence[0]

        if seq.WaveformSampleInterpretation != 'SS':
            raise ECGReadFileError(
                f"Unsupported WaveformSampleInterpretation: "
                f"'{seq.WaveformSampleInterpretation}'. Expected 'SS'."
            )
        if seq.WaveformBitsAllocated != 16:
            raise ECGReadFileError(
                f"Unsupported WaveformBitsAllocated: "
                f"{seq.WaveformBitsAllocated}. Expected 16."
            )

        for idx, definition in enumerate(seq.ChannelDefinitionSequence):
            if definition.WaveformBitsStored != 16:
                raise ECGReadFileError(
                    f"Channel {idx}: unsupported WaveformBitsStored "
                    f"{definition.WaveformBitsStored}. Expected 16."
                )

    # ------------------------------------------------------------------
    # Waveform data extraction

    def _extract_waveform(self):
        """Populate waveform attributes from the first WaveformSequence item."""
        seq = self.dicom.WaveformSequence[0]
        self.channel_definitions = seq.ChannelDefinitionSequence
        self.waveform_data = seq.WaveformData
        self.channels_no = seq.NumberOfWaveformChannels
        self.samples = seq.NumberOfWaveformSamples
        self.sampling_frequency = float(seq.SamplingFrequency)
        self.duration = self.samples / self.sampling_frequency

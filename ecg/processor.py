# -*- coding: utf-8 -*-
"""ECG signal processing: Butterworth filter and waveform extraction."""

from __future__ import annotations

import struct
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .reader import DicomECGReader
from .dsp import butter, lfilter
from .constants import FILTER_HIGHCUT, FILTER_ORDER, MILLIVOLTS


def butter_lowpass(highcut: float, sampfreq: float, order: int) -> tuple[list[float], list[float]]:
    """Compute Butterworth lowpass filter coefficients."""
    nyquist_freq = 0.5 * sampfreq
    high = highcut / nyquist_freq
    num, denom = butter(order, high, btype='lowpass')
    return num, denom


def butter_lowpass_filter(data: np.ndarray, highcut: float, sampfreq: float, order: int) -> np.ndarray:
    """Apply a Butterworth lowpass filter to waveform data."""
    num, denom = butter_lowpass(highcut, sampfreq, order=order)
    return lfilter(num, denom, data)


def extract_signals(reader: DicomECGReader) -> np.ndarray:
    """Extract, scale and filter waveform signals from a DicomECGReader.

    Returns a 2-D numpy array of shape (channels_no, samples) in millivolts.
    """
    factor = np.zeros(reader.channels_no) + 1
    baseln = np.zeros(reader.channels_no)
    units = []

    for idx in range(reader.channels_no):
        definition = reader.channel_definitions[idx]

        if definition.get('ChannelSensitivity'):
            factor[idx] = (
                float(definition.ChannelSensitivity) *
                float(definition.ChannelSensitivityCorrectionFactor)
            )

        if definition.get('ChannelBaseline'):
            baseln[idx] = float(definition.ChannelBaseline)

        units.append(
            definition.ChannelSensitivityUnitsSequence[0].CodeValue
        )

    unpack_fmt = f'<{int(len(reader.waveform_data) / 2)}h'
    unpacked = struct.unpack(unpack_fmt, reader.waveform_data)
    signals = np.asarray(unpacked, dtype=np.float32).reshape(
        reader.samples, reader.channels_no
    ).transpose()

    for channel in range(reader.channels_no):
        signals[channel] = (signals[channel] + baseln[channel]) * factor[channel]

    for i, signal in enumerate(signals):
        signals[i] = butter_lowpass_filter(
            np.asarray(signal),
            FILTER_HIGHCUT,
            reader.sampling_frequency,
            order=FILTER_ORDER,
        ) / MILLIVOLTS[units[i]]

    return signals

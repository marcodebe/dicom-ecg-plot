# -*- coding: utf-8 -*-
"""ECG (waveform) Dicom module

Read and plot images from DICOM ECG waveforms.
"""
from __future__ import annotations

from typing import BinaryIO


__author__ = "Marco De Benedetto and Simone Ferretti"
__license__ = "MIT"
__credits__ = ["Marco De Benedetto", "Simone Ferretti", "Francesco Formisano"]
__email__ = "debe@galliera.it"

try:
    from ecgconfig import LAYOUT
except ImportError:
    from .constants import DEFAULT_LAYOUT as LAYOUT

from .constants import DEFAULT_PAPER
from .reader import DicomECGReader, ECGReadFileError
from .processor import extract_signals
from .plotter import ECGPlotter


class ECG:
    """Orchestrates DICOM ECG loading, processing and rendering.

    Public API:
        ecg = ECG(source)
        ecg.draw(layoutid, mm_mv, minor_axis, interpretation)
        ecg.save(outputfile, outformat)
    """

    def __init__(self, source: str | BinaryIO | dict[str, str], paper: str = DEFAULT_PAPER) -> None:
        """Load and prepare an ECG from source.

        @param source: path string, file-like buffer, or dict with
                       keys 'stu', 'ser', 'obj' (WADO query).
        @param paper: paper format ('a4' or 'letter'), default 'a4'.
        """
        reader = DicomECGReader(source)
        signals = extract_signals(reader)
        self._plotter = ECGPlotter(reader, signals, LAYOUT, paper=paper)
        self.dicom = reader.dicom

    def draw(self, layoutid: str, mm_mv: float = 10.0, minor_axis: bool = False, interpretation: bool = False) -> None:
        """Draw grid, signals and patient info onto the figure.

        @param layoutid: layout key (e.g. '3x4_1', '3x4', '6x2', '12x1')
        @param mm_mv: amplitude scale in mm/mV
        @param minor_axis: draw 1-mm minor grid
        @param interpretation: show automated ECG interpretation text
        """
        self._plotter.draw(layoutid, mm_mv, minor_axis, interpretation)

    def save(self, outputfile: str | None = None, outformat: str | None = None) -> bytes | None:
        """Save the rendered figure.

        @param outputfile: output file path; if None, returns raw bytes.
        @param outformat: image format (pdf, png, svg, …).
        """
        return self._plotter.save(outputfile, outformat)

    def __del__(self):
        if hasattr(self, '_plotter'):
            self._plotter.close()

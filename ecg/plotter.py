# -*- coding: utf-8 -*-
"""ECG plotter: matplotlib-based rendering of ECG waveforms."""

from __future__ import annotations

import io
import re
import os
from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .reader import DicomECGReader

from matplotlib import use
if os.environ.get('DISPLAY', '') == '':
    use('Agg')

from matplotlib import pylab as plt

from . import i18n
from .constants import (
    PAPER_SIZES, DEFAULT_PAPER,
    PLOT_WIDTH, PLOT_HEIGHT, MARGIN_BOTTOM,
    GRID_COLOR, GRID_LINEWIDTH,
    SIGNAL_LINEWIDTH,
    TEXT_FONTSIZE,
    INFO_TOP_X, INFO_TOP_Y, LEGEND_TOP_X, INTERPRETATION_TOP_X,
    DURATION_BOTTOM_X, INSTITUTION_BOTTOM_X, FILTER_BOTTOM_X, INFO_BOTTOM_Y,
    LEAD_LABEL_H_OFFSET, LEAD_LABEL_V_RATIO, SEPARATOR_V_RATIO, SEPARATOR_HEIGHT,
    LEAD_LABELS,
)

try:
    from ecgconfig import INSTITUTION
except ImportError:
    from .constants import DEFAULT_INSTITUTION as INSTITUTION


def _parse_dicom_date(value: str) -> str:
    """Parse a DICOM DA value to a human-readable string.

    Handles full (YYYYMMDD), partial (YYYYMM, YYYY), and empty values.
    """
    value = (value or '').strip()
    formats = [('%Y%m%d', None), ('%Y%m', '%b %Y'), ('%Y', '%Y')]
    for fmt, out in formats:
        try:
            dt = datetime.strptime(value, fmt)
            if fmt == '%Y%m%d':
                return f"{dt.day} {dt.strftime('%b %Y')}"
            return dt.strftime(out)
        except ValueError:
            continue
    return value


def _parse_dicom_datetime(value: str) -> str:
    """Parse a DICOM DT value (YYYYMMDDHHMMSS.FFFFFF&ZZXX) to a display string.

    Strips optional fractional seconds and UTC offset before parsing.
    """
    value = (value or '').strip()
    value = re.sub(r'\.\d+', '', value)       # remove fractional seconds
    value = re.sub(r'[+-]\d{4}$', '', value)  # remove UTC offset
    formats = [
        ('%Y%m%d%H%M%S', '%d %b %Y %H:%M'),
        ('%Y%m%d%H%M',   '%d %b %Y %H:%M'),
        ('%Y%m%d',       '%d %b %Y'),
    ]
    for fmt, out in formats:
        try:
            return datetime.strptime(value, fmt).strftime(out)
        except ValueError:
            continue
    return value


def _parse_patient_age(value: str) -> str:
    """Return a numeric age string from a DICOM AS value (e.g. '045Y' → '45')."""
    value = (value or '').strip()
    if value and value[-1] in ('Y', 'M', 'W', 'D'):
        return str(int(value[:-1]))
    return value


def _lead_label(channel_source) -> str:
    """Return a short lead label from a ChannelSourceSequence item.

    Prefers a stable CodeValue lookup (MDC / SCPECG schemes) over text
    manipulation, which is vendor-dependent.  Falls back to stripping the
    common 'Lead' prefix and '(Einthoven)' suffix from CodeMeaning so that
    non-standard files still render something sensible.
    """
    label = LEAD_LABELS.get(getattr(channel_source, 'CodeValue', None))
    if label:
        return label
    meaning = getattr(channel_source, 'CodeMeaning', '') or ''
    return re.sub(r'\(Einthoven\)', '', meaning, flags=re.IGNORECASE) \
              .replace('Lead', '').replace('lead', '').strip()


class ECGPlotter:
    """Renders ECG signals onto a matplotlib figure."""

    def __init__(self, reader: DicomECGReader, signals: np.ndarray, layout_config: dict[str, list[list[int]]], paper: str = DEFAULT_PAPER) -> None:
        if paper not in PAPER_SIZES:
            raise ValueError(f"Unknown paper format '{paper}'. Valid options: {list(PAPER_SIZES)}")
        self.reader = reader
        self.signals = signals
        self.layout_config = layout_config
        self.mm_mv = None
        self._init_paper(paper)
        self.fig, self.axis = self._create_figure()

    def _init_paper(self, paper: str) -> None:
        """Compute paper-dependent layout coordinates from the chosen format."""
        spec = PAPER_SIZES[paper]
        paper_w, paper_h = spec['w_mm'], spec['h_mm']
        margin_left = 0.5 * (paper_w - PLOT_WIDTH)
        self._plot_left = margin_left / paper_w
        self._plot_right = self._plot_left + PLOT_WIDTH / paper_w
        self._plot_bottom = MARGIN_BOTTOM / paper_h
        self._plot_top = self._plot_bottom + PLOT_HEIGHT / paper_h
        self._fig_size = (spec['w_in'], spec['h_in'])

    # ------------------------------------------------------------------
    # Figure setup

    def _create_figure(self):
        fig, axes = plt.subplots()
        fig.subplots_adjust(
            left=self._plot_left, right=self._plot_right,
            top=self._plot_top, bottom=self._plot_bottom,
        )
        axes.set_ylim([0, PLOT_HEIGHT])
        axes.set_xlim([0, self.reader.samples - 1])
        return fig, axes

    # ------------------------------------------------------------------
    # Grid

    def _draw_grid(self, minor_axis):
        """Draw the ECG background grid."""
        if minor_axis:
            self.axis.xaxis.set_minor_locator(
                plt.LinearLocator(int(PLOT_WIDTH + 1))
            )
            self.axis.yaxis.set_minor_locator(
                plt.LinearLocator(int(PLOT_HEIGHT + 1))
            )

        self.axis.xaxis.set_major_locator(
            plt.LinearLocator(int(PLOT_WIDTH / 5 + 1))
        )
        self.axis.yaxis.set_major_locator(
            plt.LinearLocator(int(PLOT_HEIGHT / 5 + 1))
        )

        for axe in 'x', 'y':
            for which in 'major', 'minor':
                self.axis.grid(
                    which=which, axis=axe,
                    linestyle='-',
                    linewidth=GRID_LINEWIDTH[which],
                    color=GRID_COLOR[which],
                )
                self.axis.tick_params(
                    which=which, axis=axe,
                    color=GRID_COLOR[which],
                    bottom=False, top=False, left=False, right=False,
                )

        self.axis.set_xticklabels([])
        self.axis.set_yticklabels([])

    # ------------------------------------------------------------------
    # Signal plotting

    def _plot_signals(self, layoutid, mm_mv):
        """Plot all signal chunks according to the requested layout."""
        self.mm_mv = mm_mv
        if layoutid not in self.layout_config:
            raise ValueError(
                f"Unknown layout '{layoutid}'. Valid options: {list(self.layout_config)}"
            )
        layout = self.layout_config[layoutid]
        rows = len(layout)

        for numrow, row in enumerate(layout):
            columns = len(row)
            row_height = PLOT_HEIGHT / rows
            h_delta = self.reader.samples / columns

            v_delta = round(
                PLOT_HEIGHT * (1.0 - 1.0 / (rows * 2)) -
                numrow * (PLOT_HEIGHT / rows)
            )
            # Align origin to a multiple of 5 mm
            v_delta = (v_delta + 2.5) - (v_delta + 2.5) % 5

            chunk_size = int(self.reader.samples / len(row))
            for numcol, signum in enumerate(row):
                self._plot_signal_chunk(
                    numrow, numcol, signum, mm_mv,
                    chunk_size, h_delta, v_delta, row_height,
                )

        self.fig.set_size_inches(*self._fig_size)

    def _plot_signal_chunk(self, numrow, numcol, signum, mm_mv,
                           chunk_size, h_delta, v_delta, row_height):
        """Plot a single signal chunk and its lead label."""
        left = numcol * chunk_size
        right = (1 + numcol) * chunk_size

        signal = v_delta + mm_mv * self.signals[signum][left:right]
        self.axis.plot(
            list(range(left, right)),
            signal,
            clip_on=False,
            linewidth=SIGNAL_LINEWIDTH,
            color='black',
            zorder=10,
        )

        cseq = getattr(self.reader.channel_definitions[signum], 'ChannelSourceSequence', None)
        meaning = _lead_label(cseq[0]) if cseq else f"Ch{signum + 1}"

        h = h_delta * numcol
        v = v_delta + row_height * SEPARATOR_V_RATIO
        self._draw_lead_separator(h, v)

        self.axis.text(
            h + LEAD_LABEL_H_OFFSET,
            v_delta + row_height * LEAD_LABEL_V_RATIO,
            meaning,
            zorder=50,
            fontsize=TEXT_FONTSIZE,
        )

    def _draw_lead_separator(self, h, v):
        """Draw a short vertical line at the start of a lead chunk."""
        plt.plot(
            [h, h],
            [v - SEPARATOR_HEIGHT, v],
            lw=SIGNAL_LINEWIDTH,
            color='black',
            zorder=50,
        )

    # ------------------------------------------------------------------
    # Patient / acquisition info

    def _print_info(self, interpretation):
        """Render patient info, ECG metrics and signal parameters onto the figure."""
        plt.figtext(INFO_TOP_X, INFO_TOP_Y, self._format_patient_info(), fontsize=TEXT_FONTSIZE)
        plt.figtext(LEGEND_TOP_X, INFO_TOP_Y, self._legend(), fontsize=TEXT_FONTSIZE)

        if interpretation:
            plt.figtext(INTERPRETATION_TOP_X, INFO_TOP_Y, self._interpretation(), fontsize=TEXT_FONTSIZE)

        plt.figtext(DURATION_BOTTOM_X, INFO_BOTTOM_Y, self._format_signal_info(), fontsize=TEXT_FONTSIZE)
        plt.figtext(INSTITUTION_BOTTOM_X, INFO_BOTTOM_Y, self._format_institution_info(), fontsize=TEXT_FONTSIZE)
        plt.figtext(FILTER_BOTTOM_X, INFO_BOTTOM_Y, self._format_filter_info(), fontsize=TEXT_FONTSIZE)

    def _format_patient_info(self):
        dicom = self.reader.dicom

        try:
            pat_surname, pat_firstname = str(dicom.get('PatientName', '')).split('^')
        except ValueError:
            pat_surname = str(dicom.get('PatientName', ''))
            pat_firstname = ''

        pat_name = ' '.join((pat_surname, pat_firstname.title())).strip()
        pat_age = _parse_patient_age(dicom.get('PatientAge', ''))
        pat_id = dicom.get('PatientID', '')
        pat_sex = dicom.get('PatientSex', '')
        pat_bdate = _parse_dicom_date(dicom.get('PatientBirthDate', ''))
        acquisition_date = _parse_dicom_datetime(dicom.get('AcquisitionDateTime', ''))

        return (
            f"{pat_name}\n"
            f"{i18n.pat_id}: {pat_id}\n"
            f"{i18n.pat_sex}: {pat_sex}\n"
            f"{i18n.pat_bdate}: {pat_bdate} ({pat_age} {i18n.pat_age})\n"
            f"{i18n.acquisition_date}: {acquisition_date}"
        )

    def _format_signal_info(self):
        return f"{i18n.duration}: {self.reader.duration} s {i18n.sampling_frequency}: {self.reader.sampling_frequency} Hz"

    def _format_institution_info(self):
        institution = INSTITUTION
        if not institution:
            institution = self.reader.dicom.get('InstitutionName', '')
        return institution

    def _format_filter_info(self):
        mm_s = PLOT_WIDTH / self.reader.duration if self.reader.duration else 0
        mm_mv = self.mm_mv if self.mm_mv is not None else 0
        return f"{mm_s:.4g} mm/s {mm_mv:.4g} mm/mV 0.05-40 Hz"

    # ------------------------------------------------------------------
    # ECG annotation data

    def _legend(self):
        """Return a formatted string of ECG interval/axis measurements."""
        dicom = self.reader.dicom
        if not hasattr(dicom, 'WaveformAnnotationSequence'):
            return ''

        ecgdata = {}
        for was in dicom.WaveformAnnotationSequence:
            if was.get('ConceptNameCodeSequence'):
                cncs = was.ConceptNameCodeSequence[0]
                if cncs.CodeMeaning in (
                    'QT Interval', 'QTc Interval', 'RR Interval', 'VRate',
                    'QRS Duration', 'QRS Axis', 'T Axis', 'P Axis', 'PR Interval',
                ):
                    ecgdata[cncs.CodeMeaning] = str(was.NumericValue)

        try:
            vrate = float(ecgdata.get('VRate'))
        except (TypeError, ValueError):
            try:
                vrate = f"{60.0 / self.reader.duration * self.reader.samples / float(ecgdata.get('RR Interval')):.1f}"
            except (TypeError, ValueError, ZeroDivisionError):
                vrate = "(unknown)"

        ret_str = (
            f"{i18n.ventr_freq}: {vrate} BPM\n"
            f"{i18n.pr_interval}: {ecgdata.get('PR Interval', '')} ms\n"
            f"{i18n.qrs_duration}: {ecgdata.get('QRS Duration', '')} ms\n"
            f"{i18n.qt_qtc}: {ecgdata.get('QT Interval', '')}/{ecgdata.get('QTc Interval', '')} ms\n"
            f"{i18n.prt_axis}: {ecgdata.get('P Axis', '')} {ecgdata.get('QRS Axis', '')} {ecgdata.get('T Axis', '')}"
        )
        return ret_str

    def _interpretation(self):
        """Return the automatic ECG interpretation text."""
        dicom = self.reader.dicom
        if not hasattr(dicom, 'WaveformAnnotationSequence'):
            return ''

        lines = []
        for note in dicom.WaveformAnnotationSequence:
            if hasattr(note, 'UnformattedTextValue') and note.UnformattedTextValue:
                lines.append(note.UnformattedTextValue)
        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Public interface

    def draw(self, layoutid: str, mm_mv: float, minor_axis: bool, interpretation: bool) -> None:
        """Draw grid, signals and info onto the figure."""
        self._draw_grid(minor_axis)
        self._plot_signals(layoutid, mm_mv)
        self._print_info(interpretation)

    def save(self, outputfile: str | None = None, outformat: str | None = None) -> bytes | None:
        """Save figure to a file or return as bytes."""
        def _save(output):
            plt.savefig(output, dpi=300, format=outformat, orientation='landscape')

        if outputfile:
            _save(outputfile)
        else:
            output = io.BytesIO()
            _save(output)
            return output.getvalue()

    def close(self) -> None:
        """Release matplotlib resources."""
        plt.cla()
        plt.clf()
        plt.close()

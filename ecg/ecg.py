# -*- coding: utf-8 -*-
"""ECG (waveform) Dicom module

Read and plot images from DICOM ECG waveforms.
"""

"""
The MIT License (MIT)

Copyright (c) 2013 Marco De Benedetto <debe@galliera.it>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
import numpy as np
import pydicom as dicom
import struct
import io
import os
import requests
from . import i18n
import re
from datetime import datetime
from matplotlib import use
from scipy.signal import butter, lfilter

# python2 fails if DISPLAY is not defined with:
# _tkinter.TclError: no display name and no $DISPLAY environment variable
if os.environ.get('DISPLAY', '') == '':
    use('Agg')

from matplotlib import pylab as plt

try:
    from ecgconfig import WADOSERVER, LAYOUT, INSTITUTION
except ImportError:
    WADOSERVER = "http://example.com"
    LAYOUT = {'3x4_1': [[0, 3, 6, 9],
                        [1, 4, 7, 10],
                        [2, 5, 8, 11],
                        [1]],
              '3x4': [[0, 3, 6, 9],
                      [1, 4, 7, 10],
                      [2, 5, 8, 11]],
              '6x2': [[0, 6],
                      [1, 7],
                      [3, 8],
                      [4, 9],
                      [5, 10],
                      [6, 11]],
              '12x1': [[0],
                       [1],
                       [2],
                       [3],
                       [4],
                       [5],
                       [6],
                       [7],
                       [8],
                       [9],
                       [10],
                       [11]]}

    # If INSTITUTION is set to None the value of the tag InstitutionName is
    # used

    INSTITUTION = None

__author__ = "Marco De Benedetto and Simone Ferretti"
__license__ = "MIT"
__credits__ = ["Marco De Benedetto", "Simone Ferretti", "Francesco Formisano"]
__email__ = "debe@galliera.it"


def butter_lowpass(highcut, sampfreq, order):
    """Supporting function.

    Prepare some data and actually call the scipy butter function.
    """

    nyquist_freq = .5 * sampfreq
    high = highcut / nyquist_freq
    num, denom = butter(order, high, btype='lowpass')
    return num, denom


def butter_lowpass_filter(data, highcut, sampfreq, order):
    """Apply the Butterworth lowpass filter to the DICOM waveform.

    @param data: the waveform data.
    @param highcut: the frequencies from which apply the cut.
    @param sampfreq: the sampling frequency.
    @param order: the filter order.
    """

    num, denom = butter_lowpass(highcut, sampfreq, order=order)
    return lfilter(num, denom, data)


class ECG(object):
    """The class representing the ECG object
    """

    paper_w, paper_h = 297.0, 210.0

    # Dimensions in mm of plot area
    width = 250.0
    height = 170.0
    margin_left = margin_right = .5 * (paper_w - width)
    margin_bottom = 10.0

    # Normalized in [0, 1]
    left = margin_left / paper_w
    right = left + width / paper_w
    bottom = margin_bottom / paper_h
    top = bottom + height / paper_h

    def __init__(self, source):
        """The ECG class constructor.

        @param source: the ECG source, it could be a filename, a buffer
        or a dict of study, serie, object info (to query
        a WADO server).
        @type source: C{str} or C{dict}.
        """

        def err(msg):
            raise Exception

        def wadoget(stu, ser, obj):
            """Query the WADO server.

            @return: a buffer containing the DICOM object (the WADO response).
            @rtype: C{cStringIO.StringIO}.
            """
            payload = {
                'requestType': 'WADO',
                'contentType': 'application/dicom',
                'studyUID': stu,
                'seriesUID': ser,
                'objectUID': obj
            }
            headers = {'content-type': 'application/json'}

            resp = requests.get(WADOSERVER, params=payload, headers=headers)
            return io.BytesIO(resp.content)

        if isinstance(source, dict):
            # dictionary of stu, ser, obj
            if set(source.keys()) == set(('stu', 'ser', 'obj')):
                inputdata = wadoget(**source)
            else:
                err("source must be a dictionary of stu, ser and obj")
        elif isinstance(source, str) or getattr(source, 'getvalue'):
            # it is a filename or a (StringIO or cStringIO buffer)
            inputdata = source
        else:
            # What is it?
            err("'source' must be a path/to/file.ext string\n" +
                "or a dictionary of stu, ser and obj")

        try:
            self.dicom = dicom.read_file(inputdata)
            """@ivar: the dicom object."""
        except dicom.filereader.InvalidDicomError as err:
            raise ECGReadFileError(err)

        sequence_item = self.dicom.WaveformSequence[0]

        assert (sequence_item.WaveformSampleInterpretation == 'SS')
        assert (sequence_item.WaveformBitsAllocated == 16)

        self.channel_definitions = sequence_item.ChannelDefinitionSequence
        self.wavewform_data = sequence_item.WaveformData
        self.channels_no = sequence_item.NumberOfWaveformChannels
        self.samples = sequence_item.NumberOfWaveformSamples
        self.sampling_frequency = sequence_item.SamplingFrequency
        self.duration = self.samples / self.sampling_frequency
        self.mm_s = self.width / self.duration
        self.signals = self._signals()
        self.fig, self.axis = self.create_figure()

    def __del__(self):
        """
        Figures created through the pyplot interface
        (`matplotlib.pyplot.figure`) are retained until explicitly
        closed and may consume too much memory.
        """

        plt.cla()
        plt.clf()
        plt.close()

    def create_figure(self):
        """Prepare figure and axes"""

        # Init figure and axes
        fig = plt.figure(tight_layout=False)
        axes = fig.add_subplot(1, 1, 1)

        fig.subplots_adjust(left=self.left, right=self.right, top=self.top,
                            bottom=self.bottom)

        axes.set_ylim([0, self.height])

        # We want to plot N points, where N=number of samples
        axes.set_xlim([0, self.samples - 1])
        return fig, axes

    def _signals(self):
        """
        Retrieve the signals from the DICOM WaveformData object.

        sequence_item := dicom.dataset.FileDataset.WaveformData[n]

        @return: a list of signals.
        @rtype: C{list}
        """

        factor = np.zeros(self.channels_no) + 1
        baseln = np.zeros(self.channels_no)
        units = []
        for idx in range(self.channels_no):
            definition = self.channel_definitions[idx]

            assert (definition.WaveformBitsStored == 16)

            if definition.get('ChannelSensitivity'):
                factor[idx] = (
                    float(definition.ChannelSensitivity) *
                    float(definition.ChannelSensitivityCorrectionFactor)
                )

            if definition.get('ChannelBaseline'):
                baseln[idx] = float(definition.get('ChannelBaseline'))

            units.append(
                definition.ChannelSensitivityUnitsSequence[0].CodeValue
            )

        unpack_fmt = '<%dh' % (len(self.wavewform_data) / 2)
        unpacked_waveform_data = struct.unpack(unpack_fmt, self.wavewform_data)
        signals = np.asarray(
            unpacked_waveform_data,
            dtype=np.float32).reshape(
            self.samples,
            self.channels_no).transpose()

        for channel in range(self.channels_no):
            signals[channel] = (
                (signals[channel] + baseln[channel]) * factor[channel]
            )

        high = 40.0

        # conversion factor to obtain millivolts values
        millivolts = {'uV': 1000.0, 'mV': 1.0}

        for i, signal in enumerate(signals):
            signals[i] = butter_lowpass_filter(
                np.asarray(signal),
                high,
                self.sampling_frequency,
                order=2
            ) / millivolts[units[i]]

        return signals

    def draw_grid(self, minor_axis):
        """Draw the grid in the ecg plotting area."""

        if minor_axis:
            self.axis.xaxis.set_minor_locator(
                plt.LinearLocator(int(self.width + 1))
            )
            self.axis.yaxis.set_minor_locator(
                plt.LinearLocator(int(self.height + 1))
            )

        self.axis.xaxis.set_major_locator(
            plt.LinearLocator(int(self.width / 5 + 1))
        )
        self.axis.yaxis.set_major_locator(
            plt.LinearLocator(int(self.height / 5 + 1))
        )

        color = {'minor': '#ff5333', 'major': '#d43d1a'}
        linewidth = {'minor': .1, 'major': .2}

        for axe in 'x', 'y':
            for which in 'major', 'minor':
                self.axis.grid(
                    which=which,
                    axis=axe,
                    linestyle='-',
                    linewidth=linewidth[which],
                    color=color[which]
                )

                self.axis.tick_params(
                    which=which,
                    axis=axe,
                    color=color[which],
                    bottom=False,
                    top=False,
                    left=False,
                    right=False
                )

        self.axis.set_xticklabels([])
        self.axis.set_yticklabels([])

    def legend(self):
        """A string containing the legend.

        Auxiliary function for the print_info method.
        """

        if not hasattr(self.dicom, 'WaveformAnnotationSequence'):
            return ''

        ecgdata = {}
        for was in self.dicom.WaveformAnnotationSequence:
            if was.get('ConceptNameCodeSequence'):
                cncs = was.ConceptNameCodeSequence[0]
                if cncs.CodeMeaning in (
                        'QT Interval',
                        'QTc Interval',
                        'RR Interval',
                        'VRate',
                        'QRS Duration',
                        'QRS Axis',
                        'T Axis',
                        'P Axis',
                        'PR Interval'
                ):
                    ecgdata[cncs.CodeMeaning] = str(was.NumericValue)

        # If VRate is not defined we calculate ventricular rate from
        # RR interval
        try:
            vrate = float(ecgdata.get('VRate'))
        except (TypeError, ValueError):
            try:
                vrate = "%.1f" % (
                    60.0 / self.duration *
                    self.samples / float(ecgdata.get('RR Interval'))
                )
            except (TypeError, ValueError, ZeroDivisionError):
                vrate = "(unknown)"

        ret_str = "%s: %s BPM\n" % (i18n.ventr_freq, vrate)
        ret_str_tmpl = "%s: %s ms\n%s: %s ms\n%s: %s/%s ms\n%s: %s %s %s"

        ret_str += ret_str_tmpl % (
            i18n.pr_interval,
            ecgdata.get('PR Interval', ''),
            i18n.qrs_duration,
            ecgdata.get('QRS Duration', ''),
            i18n.qt_qtc,
            ecgdata.get('QT Interval', ''),
            ecgdata.get('QTc Interval', ''),
            i18n.prt_axis,
            ecgdata.get('P Axis', ''),
            ecgdata.get('QRS Axis', ''),
            ecgdata.get('T Axis', '')
        )

        return ret_str

    def interpretation(self):
        """Return the string representing the automatic interpretation
        of the study.
        """

        if not hasattr(self.dicom, 'WaveformAnnotationSequence'):
            return ''

        ret_str = ''
        for note in self.dicom.WaveformAnnotationSequence:
            if hasattr(note, 'UnformattedTextValue'):
                if note.UnformattedTextValue:
                    ret_str = "%s\n%s" % (
                        ret_str,
                        note.UnformattedTextValue
                    )

        return ret_str

    def print_info(self, interpretation):
        """Print info about the patient and about the ecg signals.
        """

        try:
            pat_surname, pat_firstname = str(self.dicom.PatientName).split('^')
        except ValueError:
            pat_surname = str(self.dicom.PatientName)
            pat_firstname = ''

        pat_name = ' '.join((pat_surname, pat_firstname.title()))
        pat_age = self.dicom.get('PatientAge', '').strip('Y')

        pat_id = self.dicom.PatientID
        pat_sex = self.dicom.PatientSex
        try:
            pat_bdate = datetime.strptime(
                self.dicom.PatientBirthDate, '%Y%m%d').strftime("%e %b %Y")
        except ValueError:
            pat_bdate = ""

        # Strip microseconds from acquisition date
        regexp = r"\.\d+$"
        acquisition_date_no_micro = re.sub(
            regexp, '', self.dicom.AcquisitionDateTime)

        try:
            acquisition_date = datetime.strftime(
                datetime.strptime(
                    acquisition_date_no_micro, '%Y%m%d%H%M%S'),
                '%d %b %Y %H:%M'
            )
        except ValueError:
            acquisition_date = ""

        info = "%s\n%s: %s\n%s: %s\n%s: %s (%s %s)\n%s: %s" % (
            pat_name,
            i18n.pat_id,
            pat_id,
            i18n.pat_sex,
            pat_sex,
            i18n.pat_bdate,
            pat_bdate,
            pat_age,
            i18n.pat_age,
            i18n.acquisition_date,
            acquisition_date
        )

        plt.figtext(0.08, 0.87, info, fontsize=8)

        plt.figtext(0.30, 0.87, self.legend(), fontsize=8)

        if interpretation:
            plt.figtext(0.45, 0.87, self.interpretation(), fontsize=8)

        info = "%s: %s s %s: %s Hz" % (
            i18n.duration, self.duration,
            i18n.sampling_frequency,
            self.sampling_frequency
        )

        plt.figtext(0.08, 0.025, info, fontsize=8)

        info = INSTITUTION
        if not info:
            info = self.dicom.get('InstitutionName', '')

        plt.figtext(0.38, 0.025, info, fontsize=8)

        # TODO: the lowpass filter 0.05-40 Hz will have to became a parameter
        info = "%s mm/s %s mm/mV 0.05-40 Hz" % (self.mm_s, self.mm_mv)
        plt.figtext(0.76, 0.025, info, fontsize=8)

    def save(self, outputfile=None, outformat=None):
        """Save the plot result either on a file or on a output buffer,
        depending on the input params.

        @param outputfile: the output filename.
        @param outformat: the ouput file format.
        """

        def _save(output):
            plt.savefig(
                output, dpi=300, format=outformat,
                orientation='landscape'
            )

        if outputfile:
            _save(outputfile)
        else:
            output = io.BytesIO()
            _save(output)
            return output.getvalue()

    def plot(self, layoutid, mm_mv):
        """Plot the ecg signals inside the plotting area.
        Possible layout choice are:
        * 12x1 (one signal per line)
        * 6x2 (6 rows 2 columns)
        * 3x4 (4 signal chunk per line)
        * 3x4_1 (4 signal chunk per line. on the last line
        is drawn a complete signal)
        * ... and much much more

        The general rule is that the layout list is formed
        by as much lists as the number of lines we want to draw into the
        plotting area, each one containing the number of the signal chunk
        we want to plot in that line.

        @param layoutid: the desired layout
        @type layoutid: C{list} of C{list}
        """

        self.mm_mv = mm_mv

        layout = LAYOUT[layoutid]
        rows = len(layout)

        for numrow, row in enumerate(layout):

            columns = len(row)
            row_height = self.height / rows

            # Horizontal shift for lead labels and separators
            h_delta = self.samples / columns

            # Vertical shift of the origin
            v_delta = round(
                self.height * (1.0 - 1.0 / (rows * 2)) -
                numrow * (self.height / rows)
            )

            # Let's shift the origin on a multiple of 5 mm
            v_delta = (v_delta + 2.5) - (v_delta + 2.5) % 5

            # Lenght of a signal chunk
            chunk_size = int(self.samples / len(row))
            for numcol, signum in enumerate(row):
                left = numcol * chunk_size
                right = (1 + numcol) * chunk_size

                # The signal chunk, vertical shifted and
                # scaled by mm/mV factor
                signal = v_delta + mm_mv * self.signals[signum][left:right]
                self.axis.plot(
                    list(range(left, right)),
                    signal,
                    clip_on=False,
                    linewidth=0.6,
                    color='black',
                    zorder=10)

                cseq = self.channel_definitions[signum].ChannelSourceSequence
                meaning = cseq[0].CodeMeaning.replace(
                    'Lead', '').replace('(Einthoven)', '')

                h = h_delta * numcol
                v = v_delta + row_height / 2.6
                plt.plot(
                    [h, h],
                    [v - 3, v],
                    lw=.6,
                    color='black',
                    zorder=50
                )

                self.axis.text(
                    h + 40,
                    v_delta + row_height / 3,
                    meaning,
                    zorder=50,
                    fontsize=8
                )

        # A4 size in inches
        self.fig.set_size_inches(11.69, 8.27)

    def draw(self, layoutid, mm_mv=10.0, minor_axis=False, interpretation=False):
        """Draw grid, info and signals"""

        self.draw_grid(minor_axis)
        self.plot(layoutid, mm_mv)
        self.print_info(interpretation)


class ECGReadFileError(dicom.filereader.InvalidDicomError):
    pass

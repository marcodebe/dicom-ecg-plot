#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ECG Conversion Tool

Usage:
    ecgconvert.py <inputfile> [--layout=LAYOUT] [--usetex] [--output=FILE|--format=FMT]
    ecgconvert.py <stu> <ser> <obj> [--layout=LAYOUT] [--usetex] [--output=FILE|--format=FMT]

Options:
    -h, --help                 This help.
    <inputfile>                Input dicom file.
    <stu> <ser> <obj>          studyUID seriesUID objectUID
                               UIDs for WADO download.
    -l LAYOUT --layout=LAYOUT  Layout [default: 3X4_1].
    -o FILE --output=FILE      Output file (format deduced by extension).
    -f FMT --format=FMT        Output format.
    --usetex                   Use TeX for better fonts (25% slower).

Valid layouts are: 3X4_1, 3X4, 12X1

The output format is deduced from the extension of the filename.

Supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba,
                   svg, svgz, tif, tiff.
"""

from docopt import docopt
from datetime import datetime
from matplotlib import pylab as plt
from scipy.signal import butter, lfilter
import numpy as np
import dicom
import struct
import cStringIO
import requests

from config import WADOSERVER, LAYOUT

__author__ = "Marco De Benedetto"
__email__ = "debe@galliera.it"


def butter_bandpass(lowcut, highcut, size, order=5):
    nyquist_freq = .5 * size
    low = lowcut / nyquist_freq
    high = highcut / nyquist_freq
    num, denom = butter(order, [low, high], btype='band')
    return num, denom


def butter_bandpass_filter(data, lowcut, highcut, size, order):
    num, denom = butter_bandpass(lowcut, highcut, size, order=order)
    return lfilter(num, denom, data)


class ECG(object):

    paper_w, paper_h = 297.0, 210.0

    # Dimensions in mm of plot area
    width = 250.0
    height = 170.0
    margin_left = margin_right = .5 * (paper_w - width)
    margin_bottom = 10

    # Normalized in [0, 1]
    left = margin_left/paper_w
    right = left+width/paper_w
    bottom = margin_bottom/paper_h
    top = bottom+height/paper_h

    def __init__(self, filename=None, stu=None, ser=None, obj=None):

        def wadoget(stu, ser, obj):

            payload = {'requestType': 'WADO', 'studyUID': stu,
                       'seriesUID': ser, 'objectUID': obj}
            headers = {'content-type': 'application/json'}

            resp = requests.get(WADOSERVER, params=payload, headers=headers)
            return cStringIO.StringIO(resp.content)

        if filename:
            inputdata = filename
        else:
            inputdata = wadoget(stu, ser, obj)

        self.dicom = dicom.read_file(inputdata)

        sequence_item = self.dicom.WaveformSequence[0]

        self.channel_definitions = sequence_item.ChannelDefinitionSequence
        self.wavewform_data = sequence_item.WaveformData
        self.channels_no = sequence_item.NumberOfWaveformChannels
        self.samples = sequence_item.NumberOfWaveformSamples

        self.signals = self._signals()
        self.fig, self.axis = self.create_figure()

    def create_figure(self):
        """
        Prepare figure and axes
        """

        # Init figure and axes
        fig = plt.figure(tight_layout=False)
        axes = fig.add_subplot(1, 1, 1)
        # axes.set_frame_on(False)
        fig.subplots_adjust(left=self.left, right=self.right,
                            top=self.top, bottom=self.bottom)
        axes.set_ylim([0, self.height])
        # we want to plot N points, where N=number of samples
        axes.set_xlim([0, self.samples-1])
        return fig, axes

    def _signals(self):
        """
        sequence_item := dicom.dataset.FileDataset.WaveformData[n]
        Return a list of signals.
        """

        factor = np.zeros(self.channels_no) + 1
        for idx in range(self.channels_no):
            definition = self.channel_definitions[idx]
            if definition.get('ChannelSensitivity'):
                factor[idx] = (
                    float(definition.ChannelSensitivity) *
                    float(definition.ChannelSensitivityCorrectionFactor))

        signals = np.asarray(
            struct.unpack('<' + str(len(self.wavewform_data)/2) + 'h',
                          self.wavewform_data), dtype=np.float32).reshape(
                              self.samples, self.channels_no).transpose()

        for channel in range(self.channels_no):
            signals[channel] *= factor[channel]

        low = .05
        high = 40.0
        for i, signal in enumerate(signals):
            # micro to milli volt
            signals[i] = butter_bandpass_filter(np.asarray(signal),
                                                low, high, 1000, order=1)/1000

        return signals

    def draw_grid(self):

        #self.axis.xaxis.set_minor_locator(plt.LinearLocator(self.width+1))
        self.axis.xaxis.set_major_locator(plt.LinearLocator(self.width/5+1))
        #self.axis.yaxis.set_minor_locator(plt.LinearLocator(self.height+1))
        self.axis.yaxis.set_major_locator(plt.LinearLocator(self.height/5+1))

        color = {'minor': '#ff5333', 'major': '#d43d1a'}
        linewidth = {'minor': .2, 'major': .6}
        alpha = 1

        for axe in 'x', 'y':
            for which in 'major', 'minor':
                self.axis.grid(which=which, axis=axe,
                               linewidth=linewidth[which],
                               linestyle='-', color=color[which], alpha=alpha)
                self.axis.tick_params(which=which, axis=axe,
                                      color=color[which],
                                      bottom='off', top='off',
                                      left='off', right='off')

        self.axis.set_xticklabels([])
        self.axis.set_yticklabels([])

    def legend(self):

        ecgdata = {}
        for was in self.dicom.WaveformAnnotationSequence:
            if was.get('ConceptNameCodeSequence'):
                cncs = was.ConceptNameCodeSequence[0]
                if cncs.CodeMeaning in ('QT Interval',
                                        'QTc Interval',
                                        'RR Interval',
                                        'QRS Duration',
                                        'QRS Axis',
                                        'T Axis',
                                        'P Axis',
                                        'PR Interval'):
                    ecgdata[cncs.CodeMeaning] = str(was.NumericValue)

        return 'Ventr. Freq.: ' + \
               str(60000 / int(ecgdata['RR Interval'])) + ' BPM\n' + \
               'PR Interval: ' + ecgdata['PR Interval'] + ' ms\n' + \
               'QRS Duration: ' + ecgdata['QRS Duration'] + ' ms\n' + \
               'QT/QTc: ' + \
               ecgdata['QT Interval'] + '/' + \
               ecgdata['QTc Interval'] + ' ms\n' + \
               'P-R-T Axis: ' + \
               ecgdata['P Axis'] + ' ' + \
               ecgdata['QRS Axis'] + ' ' + \
               ecgdata['T Axis'] + ' ms'

    def print_info(self):

        try:
            pat_surname, pat_firstname = self.dicom.PatientName.split('^')
        except ValueError:
            pat_surname = self.dicom.PatientName
            pat_firstname = ''

        pat_name = ' '.join((pat_surname, pat_firstname.title()))

        pat_id = self.dicom.PatientID
        pat_sex = self.dicom.PatientSex
        text_y = self.height+18

        ecg_date_str = (self.dicom.InstanceCreationDate +
                        self.dicom.InstanceCreationTime)
        ecg_date = datetime.strftime(datetime.strptime(ecg_date_str,
                                                       '%Y%m%d%H%M%S'),
                                     '%d %b %Y %H:%M')
        patient_str = "%s (%s) sex: %s" % (pat_name, pat_id, pat_sex)
        self.axis.text(0, text_y, patient_str, fontsize=12)
        self.axis.text(0, text_y-5, "ECG date: " + ecg_date, fontsize=10)

    def save(self, outputfile=None, outformat=None):

        def _save(output):
            plt.savefig(output, dpi=300, format=outformat,
                        papertype='a4', orientation='landscape')

        if outputfile:
            _save(outputfile)
        else:
            output = cStringIO.StringIO()
            _save(output)
            return output.getvalue()

    def plot(self, layout):

        rows = len(layout)

        for numrow, row in enumerate(layout):
            columns = len(row)
            h_delta = self.samples / columns
            signal = np.ndarray(0)
            row_height = self.height / rows
            v_delta = round(self.height * (1 - 1.0/(rows*2)) -
                            numrow*(self.height/rows))

            v_delta = (v_delta + 2.5) - (v_delta + 2.5) % 5
            chunk_size = int(self.samples/len(row))

            for numcol, signum in enumerate(row):
                left = numcol*chunk_size
                right = (1+numcol)*chunk_size
                h = h_delta * numcol
                plt.plot([h, h], [v_delta-row_height/2.6,
                                  v_delta+row_height/2.6],
                         'k-', lw=1, color='blue', zorder=50)
                signal = np.concatenate((
                    signal,
                    10.0*self.signals[signum][left:right])
                )
                cseq = self.channel_definitions[signum].ChannelSourceSequence
                meaning = cseq[0].CodeMeaning.replace(
                    'Lead', '').replace('(Einthoven)', '')
                self.axis.text(h+40, v_delta+row_height/3,
                               meaning, color='b', zorder=50)

            self.axis.text(4000, self.height+2, self.legend(),
                           fontsize=10, color='k', zorder=50)

            self.axis.plot(signal+v_delta, linewidth=.6, color='black',
                           antialiased=True, zorder=10)

        # A4 size in inches
        self.fig.set_size_inches(11.69, 8.27)


def convert(source, layout, outformat=None, outputfile=None):

    def err(msg):
        raise(Exception(msg))

    if isinstance(source, dict):
        if len(source) == 3:
            ecg = ECG(**source)
        else:
            err("source must be a dictionary of stu, ser and obj")
    elif isinstance(source, basestring):
        ecg = ECG(filename=source)
    else:
        err("`source´ must be a path/to/file.ext string\n" +
            "or a dictionary of stu, ser and obj")

    ecg.draw_grid()
    ecg.print_info()
    ecg.plot(layout)
    return ecg.save(outformat=outformat, outputfile=outputfile)

if __name__ == '__main__':

    arguments = docopt(__doc__, version='ECG Convert 0.1')
    inputfile = arguments['<inputfile>']
    stu = arguments['<stu>']
    ser = arguments['<ser>']
    obj = arguments['<obj>']
    outputfile = arguments['--output']
    usetex = arguments['--usetex']
    outformat = arguments['--format']
    layout = LAYOUT[arguments['--layout']]

    if usetex:
        plt.rc('text', usetex=True)

    if inputfile:
        source = inputfile
        output = convert(inputfile, layout, outformat, outputfile)
    else:
        source = {'stu': stu, 'ser': ser, 'obj': obj}

    output = convert(source, layout, outformat, outputfile)

    if output:
        print output

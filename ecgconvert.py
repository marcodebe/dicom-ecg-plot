#!/usr/bin/env python

"""ECG Conversion Tool

Usage:
    ecgconvert.py <inputfile> -o <outputfile>

Options:
    -h, --help      This help
    <inputfile>     Input dicom file
    -o <outpufile>  Output file

The output format is deduced from the extension of the filename.

Supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba,
                   svg, svgz, tif, tiff.
"""

from docopt import docopt
from matplotlib import pylab as plt
from scipy.signal import butter, lfilter
import numpy as np
import dicom

__author__ = "Marco De Benedetto"
__email__ = "debe@galliera.it"


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyquist_freq = 0.5 * fs
    low = lowcut / nyquist_freq
    high = highcut / nyquist_freq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def norm(x, max_value):
    sgn = lambda x: lambda y: x-y > 0 and 1 or -1
    norm = lambda y: lambda x: x % (sgn(y)(x) * y)
    return norm(max_value)(x)


class ECG(object):

    paper_w, paper_h = 297.0, 210.0

    # Dimensions in mm of plot area
    width = 250
    height = 170
    margin_left = margin_right = .5 * (paper_w - width)
    margin_bottom = 10

    # Normalized in [0, 1]
    left = margin_left/paper_w
    right = left+width/paper_w
    bottom = margin_bottom/paper_h
    top = bottom+height/paper_h

    def create_figure(self):
        """
        Prepare figure and axes
        """

        # Init figure and axes
        fig = plt.figure(tight_layout=False)
        ax = fig.add_subplot(1, 1, 1)
        # ax.set_frame_on(False)
        fig.subplots_adjust(left=self.left, right=self.right,
                            top=self.top, bottom=self.bottom)
        ax.set_ylim([0, self.height])
        # we want to plot N points, where N=number of samples
        ax.set_xlim([0, self.samples-1])
        return fig, ax

    def get_signals(self, sequence_item):
        """
        sequence_item := dicom.dataset.FileDataset.WaveformData[n]
        Return a list of signals.
        """

        channel_definitions = sequence_item.ChannelDefinitionSequence
        wavewform_data = sequence_item.WaveformData
        channels_no = sequence_item.NumberOfWaveformChannels

        self.samples = sequence_item.NumberOfWaveformSamples

        signals = []
        factor = np.zeros(channels_no) + 1
        for idx in range(channels_no):
            signals.append([])

            definition = channel_definitions[idx]
            if definition.get('ChannelSensitivity'):
                factor[idx] = (
                    float(definition.ChannelSensitivity) *
                    float(definition.ChannelSensitivityCorrectionFactor))

        for idx in xrange(0, len(wavewform_data), 2):
            channel = int((idx/2) % channels_no)
            signals[channel].append(factor[channel] *
                                    norm(ord(wavewform_data[idx]) +
                                         256 * ord(wavewform_data[idx+1]),
                                         2**15))

        low = .05
        high = 40.0
        for i, s in enumerate(signals):
            # micro to milli volt
            signals[i] = butter_bandpass_filter(np.asarray(s),
                                                low, high, 1000, order=1)/1000

        return signals

    def __init__(self, filename):
        self.filename = filename
        self.dicom = dicom.read_file(self.filename)
        self.signals = self.get_signals(self.dicom.WaveformSequence[0])
        self.fig, self.ax = self.create_figure()

    def draw_grid(self):

        self.ax.xaxis.set_minor_locator(plt.LinearLocator(self.width+1))
        self.ax.xaxis.set_major_locator(plt.LinearLocator(self.width/5+1))
        self.ax.yaxis.set_minor_locator(plt.LinearLocator(self.height+1))
        self.ax.yaxis.set_major_locator(plt.LinearLocator(self.height/5+1))

        color = {'minor': '#ff5333', 'major': '#d43d1a'}
        linewidth = {'minor': .2, 'major': .6}

        for axis in 'x', 'y':
            for which in 'major', 'minor':
                self.ax.grid(which=which, axis=axis,
                             linewidth=linewidth[which],
                             linestyle='-', color=color[which], alpha=.5)
                self.ax.tick_params(which=which, axis=axis, color=color[which],
                                    bottom='off', top='off',
                                    left='off', right='off')

        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])

    def print_info(self):

        pat_name = self.dicom.PatientName.replace('^', ' ')
        pat_id = self.dicom.PatientID
        pat_sex = self.dicom.PatientSex
        text_y = self.height+23

        from datetime import datetime
        ecg_date_str = (self.dicom.InstanceCreationDate +
                        self.dicom.InstanceCreationTime)
        ecg_date = datetime.strftime(datetime.strptime(ecg_date_str,
                                                       '%Y%m%d%H%M%S'),
                                     '%d %b %Y %H:%M')
        patient_str = "%s (%s) sex: %s" % (pat_name, pat_id, pat_sex)
        self.ax.text(0, text_y, patient_str, fontsize=12)
        self.ax.text(0, text_y-5, "ECG date: " + ecg_date, fontsize=10)

    def plot(self, outputfile):

        self.signals.reverse()

        premax = delta = 0
        delta = 3
        for num, signal in enumerate(self.signals):
            delta += 1 + 10 * (premax - signal.min())
            self.ax.plot(10.0*signal+delta, linewidth=0.6, color='black',
                         antialiased=True, zorder=10)
            premax = signal.max()

        # A4 size in inches
        self.fig.set_size_inches(11.69, 8.27)

        plt.savefig(outputfile, dpi=300)

if __name__ == '__main__':

    arguments = docopt(__doc__, version='ECG Convert 0.1')
    ecg = ECG(arguments['<inputfile>'])
    ecg.draw_grid()
    ecg.print_info()
    ecg.plot(arguments['-o'])

# -*- coding: utf-8 -*-
"""Constants for ECG plotting layout and rendering."""

# Supported paper formats (landscape orientation)
# w_mm/h_mm: dimensions in millimetres; w_in/h_in: dimensions in inches
PAPER_SIZES: dict = {
    'a4':     {'w_mm': 297.0,  'h_mm': 210.0,  'w_in': 11.69, 'h_in': 8.27},
    'letter': {'w_mm': 279.4,  'h_mm': 215.9,  'w_in': 11.0,  'h_in': 8.5},
}
DEFAULT_PAPER = 'a4'

# ECG plot area dimensions (mm) — fixed regardless of paper format
PLOT_WIDTH = 250.0
PLOT_HEIGHT = 170.0
MARGIN_BOTTOM = 10.0

# Butterworth lowpass filter
FILTER_HIGHCUT = 40.0
FILTER_ORDER = 2

# ECG grid colors and linewidths
GRID_COLOR = {'minor': '#ff5333', 'major': '#d43d1a'}
GRID_LINEWIDTH = {'minor': 0.1, 'major': 0.2}

# Signal rendering
SIGNAL_LINEWIDTH = 0.6

# Unit conversion to millivolts
MILLIVOLTS = {'uV': 1000.0, 'mV': 1.0}

# Text rendering
TEXT_FONTSIZE = 8
INFO_TOP_X = 0.08
INFO_TOP_Y = 0.87
LEGEND_TOP_X = 0.30
INTERPRETATION_TOP_X = 0.45
DURATION_BOTTOM_X = 0.08
INSTITUTION_BOTTOM_X = 0.38
FILTER_BOTTOM_X = 0.76
INFO_BOTTOM_Y = 0.025

# Lead label positioning relative to row height
LEAD_LABEL_H_OFFSET = 40
LEAD_LABEL_V_RATIO = 1 / 3
SEPARATOR_V_RATIO = 1 / 2.6
SEPARATOR_HEIGHT = 3

# Default ECG layouts: each layout is a list of rows,
# each row is a list of lead indices (0-based).
DEFAULT_LAYOUT = {
    '3x4_1': [[0, 3, 6, 9],
               [1, 4, 7, 10],
               [2, 5, 8, 11],
               [1]],
    '3x4':   [[0, 3, 6, 9],
               [1, 4, 7, 10],
               [2, 5, 8, 11]],
    '6x2':   [[0, 6],
               [1, 7],
               [2, 8],
               [3, 9],
               [4, 10],
               [5, 11]],
    '12x1':  [[0], [1], [2], [3], [4], [5],
               [6], [7], [8], [9], [10], [11]],
}

DEFAULT_WADOSERVER = "http://example.com"
DEFAULT_INSTITUTION = None

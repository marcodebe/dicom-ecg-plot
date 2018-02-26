# -*- coding: utf-8 -*-
"""ECG Conversion Tool

Usage:
    convert.py <inputfile> [--layout=LAYOUT] [--output=FILE|--format=FMT] [--minor-grid] [--interpretation]
    convert.py <stu> <ser> <obj> [--layout=LAYOUT] [--output=FILE|--format=FMT] [--minor-grid] [--interpretation]
    convert.py --help

Options:
    -h, --help                 This help.
    <inputfile>                Input dicom file.
    <stu> <ser> <obj>          studyUID seriesUID objectUID
                               UIDs for WADO download.
    -l LAYOUT --layout=LAYOUT  Layout [default: 3x4_1].
    -o FILE --output=FILE      Output file (format deduced by extension).
    -f FMT --format=FMT        Output format.
    --minor-grid               Draw minor axis grid (1mm).
    --interpretation           Show "Automated ECG interpretation"

Valid layouts are: 3x4_1, 3x4, 12x1

The output format is deduced from the extension of the filename, if present, or
from --format option when filename is not specified.

Valid formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz,
               tif, tiff.
"""
from __future__ import unicode_literals, print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import open
from ecg import ECG
from docopt import docopt
from io import BytesIO


def convert(source, layout, outformat, outputfile,
            minor_axis=False, interpretation=False):

    ecg = ECG(source)
    ecg.draw(layout, 10, minor_axis, interpretation=interpretation)
    return ecg.save(outformat=outformat, outputfile=outputfile)

if __name__ == '__main__':

    arguments = docopt(__doc__, version='ECG Convert 0.1')
    inputfile = arguments['<inputfile>']
    stu = arguments['<stu>']
    ser = arguments['<ser>']
    obj = arguments['<obj>']
    outputfile = arguments['--output']
    outformat = arguments['--format']
    layout = arguments['--layout']
    minor_axis = arguments['--minor-grid']
    interpretation = arguments['--interpretation']

    if inputfile:
        source = BytesIO(open(inputfile, mode='rb').read())
    else:
        source = {'stu': stu, 'ser': ser, 'obj': obj}

    output = convert(source, layout, outformat, outputfile, minor_axis, interpretation)

    if output:
        print(output)

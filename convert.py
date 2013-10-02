# -*- coding: utf-8 -*-
"""ECG Conversion Tool

Usage:
    convert.py <inputfile> [--layout=LAYOUT] [--output=FILE|--format=FMT]
    convert.py <stu> <ser> <obj> [--layout=LAYOUT] [--output=FILE|--format=FMT]
    convert.py --help

Options:
    -h, --help                 This help.
    <inputfile>                Input dicom file.
    <stu> <ser> <obj>          studyUID seriesUID objectUID
                               UIDs for WADO download.
    -l LAYOUT --layout=LAYOUT  Layout [default: 3X4_1].
    -o FILE --output=FILE      Output file (format deduced by extension).
    -f FMT --format=FMT        Output format.

Valid layouts are: 3X4_1, 3X4, 12X1

The output format is deduced from the extension of the filename, if present, or
from --format option when filename is not specified.

Valid formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz,
               tif, tiff.
"""

from ecg import ECG
from docopt import docopt


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
    outformat = arguments['--format']
    layout = arguments['--layout']

    if inputfile:
        source = inputfile
        output = convert(inputfile, layout, outformat, outputfile)
    else:
        source = {'stu': stu, 'ser': ser, 'obj': obj}

    output = convert(source, layout, outformat, outputfile)

    if output:
        print output

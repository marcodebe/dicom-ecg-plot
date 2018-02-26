from __future__ import unicode_literals, print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()
WADOSERVER = "http://example.com/"

LAYOUT = {'3x4_1': [[0, 3, 6, 9],
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
          '12x1':  [[0],
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

# If INSTITUTION is set to None the value of the tag InstitutionName is used
INSTITUTION = None

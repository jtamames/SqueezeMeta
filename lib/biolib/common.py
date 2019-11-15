###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import errno
import sys
import logging
import ntpath
import re
import gzip

from numpy import (abs as np_abs,
                   array as np_array)


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input

    Parameters
    ----------
    question : str
        Prompt presented to the user.
    default : str
        Presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    Returns
    -------
    boolean
        True for "yes", False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def is_float(s):
    """Check if a string can be converted to a float.

    Parameters
    ----------
    s : str
        String to evaluate.

    Returns
    -------
    boolean
        True if string can be converted, else False.
    """

    try:
        float(s)
    except ValueError:
        return False

    return True


def find_nearest(array, value):
    """Find nearest array element to a given value.

    Parameters
    ----------
    array : iterable
        List of values.
    value : float
        Target value

    Returns
    -------
    array element
        Closest element in 'array' to 'value'.

    """
    idx = (np_abs(np_array(array) - value)).argmin()
    return array[idx]


def concatenate_files(input_files, output_file, common_header=False):
    """Concatenate several files into a single file.

    Creates a compressed file if the extension of
    the output file ends with .gz.

    Parameters
    ----------
    input_files : iterable
        Files to concatenate.
    output_file : str
        Name of output file.
    common_header : boolean
        If True, only write first line of first file.
    """

    if output_file.endswith('.gz'):
        open_file = gzip.open
    else:
        open_file = open

    with open_file(output_file, "wb") as outfile:
        for file_no, f in enumerate(input_files):
            with open(f, "rb") as infile:
                for line_no, line in enumerate(infile):
                    if common_header and line_no == 0:
                        if file_no == 0:
                            outfile.write(line)
                        else:
                            pass
                    else:
                        outfile.write(line)
                        
                # force newline between files
                if line[-1] != '\n':
                    outfile.write('\n')

                        
def alphanumeric_sort(l):
    """Sorts the given iterable alphanumerically.

    http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python

    Parameters
    ----------
    l : iterable
        The iterable to be sorted.

    Returns
    -------
    iterable
        Iterable sorted alphanumerically.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def check_file_exists(input_file):
    """Check if file exists."""
    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exists: ' + input_file + '\n')
        sys.exit()


def check_dir_exists(input_dir):
    """Check if directory exists."""
    if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
        logger = logging.getLogger('timestamp')
        logger.error('Input directory does not exists: ' + input_dir + '\n')
        sys.exit()


def make_sure_path_exists(path):
    """Create directory if it does not exist."""

    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created: ' + path + '\n')
            sys.exit()


def remove_extension(filename, extension=None):
    """Remove extension from filename.

    A specific extension can be specified, otherwise
    the extension is taken as all characters after the
    last period.
    """
    f = ntpath.basename(filename)

    if extension and f.endswith(extension):
        f = f[0:f.rfind(extension)]
    else:
        f = os.path.splitext(f)[0]

    if f[-1] == '.':
        f = f[0:-1]

    return f

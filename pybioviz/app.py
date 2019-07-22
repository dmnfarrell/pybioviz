#!/usr/bin/env python
"""
    Implements sequence utilities for pybioviz
    Created June 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import os,sys,subprocess
import pandas as pd
from . import viewers, utils

def main():
    "Run the application"

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-b", "--bam-viewer", dest="bam",  action="store_true",
                        default=False, help="Run bam viewer")
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Do tests")
    opts, remainder = parser.parse_args()
    if opts.bam is True:
        cmd = 'panel serve'
    elif opts.test is True:
        cmd = 'panel serve notebooks/tests.ipynb'
        subprocess.check_output(cmd, shell=True)
if __name__ == '__main__':
    main()

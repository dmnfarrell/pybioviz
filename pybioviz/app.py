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
from . import dashboards, utils, __version__

def run_server(name, path='.', port=8000):
    
    if name == 'test':
        app = dashboards.test_app()
    elif name == 'align':
        app = dashboards.sequence_alignment_viewer(filename)
    elif name == 'bam':
        app = bam_view('wt_mbovis.bam', 'Mbovis-AF212297.fa', 'Mbovis_AF212297.gff', width=1000)    
    else:
        app = dashboards.view_features(path)    
    from bokeh.server.server import Server
    def modify_doc(doc):               
        return app.server_doc(doc=doc, title='pybioviz: %s' %name)
    
    print('Opening application on http://localhost:%s/' %port)
    server = Server({'/': modify_doc}, port=port)
    server.start()
    server.show('/')
    server.run_until_shutdown()
    return

def main():
    "Run the application"

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-a", "--align", dest="seq",  action="store_true",
                        default=False, help="Run seq alignment viewer")
    parser.add_option("-b", "--bam-viewer", dest="bam",  action="store_true",
                        default=False, help="Run bam viewer")    
    parser.add_option("-p", "--port", dest="port", default=8000,
                        help="Port for web app, default 8000")    
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="Show test plots")
    opts, remainder = parser.parse_args()
    if opts.bam is True:
        cmd = 'panel serve'
             
    elif opts.test is True:
        run_server('test')
        
if __name__ == '__main__':
    main()

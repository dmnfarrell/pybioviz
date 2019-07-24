#!/usr/bin/env python
"""
    Implements sequence viewers for pybioviz
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

import os,sys,random
import numpy as np
import pandas as pd
from . import utils
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn

def test_plot(rows=20, cols=100):
    """Bokeh random colors plot"""

    x = np.arange(1, cols)
    y = np.arange(0, rows, 1)    
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    
    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(len(gx))]
    source = ColumnDataSource(dict(x=gx, y=gy, color=colors))    
    plot_height = 50
    x_range = Range1d(0,10, bounds='auto')

    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width=800, plot_height=200,
                    tools="pan")
    rects = Rect(x="x", y="y",  width=1, height=1, fill_color="color", line_width=0)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  
    #p.title('test bokeh plot')
    p = gridplot([[p]], toolbar_location='below')
    return p
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

import os,sys
import numpy as np
import pandas as pd
from . import utils, plotters
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn
import panel.widgets as pnw


def view_features(features=None):
    """gene feature viewer app"""
    
    #features = utils.gff_to_features('Mbovis_AF212297.gff')
    gff_input = pnw.TextInput(name='gff file',value='')
    loc_input = pnw.TextInput(name='location',value='',width=200)
    gene_input = pnw.TextInput(name='find_gene',value='',width=200)
    zoomout_btn = pnw.Button(name='-',width=40, button_type='primary')
    #load_btn = pn.widgets.FileInput()
    slider = pnw.IntRangeSlider(start=0,end=1000000,step=10,value=(1,20000),width=900)
    feature_pane = pn.pane.Bokeh(height=100,margin=10)
    found = None
    if features is None:
        features = utils.gff_to_features(gff_input.value)
    
    def load_file(event):
        nonlocal features
        features = utils.gff_to_features(gff_input.value)
        #feature_pane.object = plot_features(features,preview=False,x_range=xrange, plot_width=900)
        update(event)
        
    def find_gene(event):
        gene = gene_input.value         
        df = utils.features_to_dataframe(features).fillna('-')
        found = df[df.gene.str.contains(gene)].iloc[0]        
        loc = (found.start-200,found.end+200)
        slider.value = loc
        #feature_pane.object = view_features(features,preview=False,x_range=loc, plot_width=900)
        return
    
    def zoomout(event):
        
        return
    
    def update(event):    
        xrange = slider.value
        loc_input.value = str(xrange[0])+':'+str(xrange[1])
        #p1 = annot_pane.object = preview(features)
        feature_pane.object = plotters.plot_features(features,preview=False,x_range=xrange, plot_width=900)
        return

    slider.param.watch(update,'value')
    slider.param.trigger('value')
    gene_input.param.watch(find_gene,'value')
    gff_input.param.watch(load_file,'value')
    zoomout_btn.param.watch(zoomout,'clicks')
    test_pane = pn.pane.Str(50,width=180,style={'margin': '4pt'})
    
    top=pn.Row(gff_input,loc_input,gene_input,zoomout_btn,test_pane)
    main = pn.Column(feature_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main)
    return app



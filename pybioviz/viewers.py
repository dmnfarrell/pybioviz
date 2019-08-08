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

def view_sequence_alignment(aln, fontsize="8pt", plot_width=800):
    """Bokeh sequence alignment viewer.
    Args:
        aln: biopython Multiple Sequence Alignment        
    """
    
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    #ids=range(len(seqs))
    text = [i for s in list(seqs) for i in s]
    colors = utils.get_sequence_colors(seqs)    
    cons = utils.get_cons(aln)
    N = len(seqs[0])
    S = len(seqs)    
    width=.4
    
    x = np.arange(1, N+1)
    y = np.arange(0,S,1)
    print (y[:20])
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy+.5
    h= 1/S
    print (N,S)
    #print (text)
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"
    
    #box_select = BoxSelectTool(callback=callback_select)

    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width=plot_width, plot_height=50, x_range=x_range, y_range=(0,S), tools=tools, 
                    min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    #p.xaxis.visible = False
    p.yaxis.visible = False
    p.grid.visible = False  
    
    #sequence text view, zoom fixed
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height, x_range=view_range, y_range=ids, tools="xpan,reset", 
                    min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black", text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)
            
    p1.grid.visible = False
    p.toolbar.logo = None
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    
    source2 = ColumnDataSource(dict(x=x, cons=cons))
    
    p3 = figure(title=None, plot_width=plot_width, plot_height=30, x_range=p1.x_range, y_range=(Range1d(min(cons),.5)), tools="xpan")
    rects2 = Rect(x="x", y=0,  width=1, height="cons", fill_color="gray", line_color=None, fill_alpha=0.7)
    p3.add_glyph(source2, rects2)
    
    p3.xaxis.visible = False
    p3.yaxis.visible = False
    p3.grid.visible = False    
    p3.background_fill_color = "beige"

    #p.js_on_change('selected', callback_select)
    jscode="""    
    var start = cb_obj.value;    
    x_range.setv({"start": start, "end": start+l})   
    """
    callback = CustomJS(
        args=dict(x_range=p1.x_range,l=viewlen), code=jscode)
    slider = Slider (start=0, end=N, value=1, step=10)
    slider.js_on_change('value', callback)
    
    p = gridplot([[p],[slider],[p3],[p1]], toolbar_location='below')
    return p

def view_features():
    """Genome features viewer app"""
        
    gff_input = pnw.TextInput(name='gff file',value='Mbovis_AF212297.gff')
    loc_input = pnw.TextInput(name='location',value='',width=200)
    gene_input = pnw.TextInput(name='find_gene',value='',width=200)
    next_btn = pnw.Button(name='\u25b6',width=40, button_type='primary')
    #load_btn = pn.widgets.FileInput()
    slider = pnw.IntRangeSlider(start=0,end=1000000,step=10,value=(1,20000),width=900)
    feature_pane = pn.pane.Bokeh(height=100,margin=10)
    found = None
    features = utils.gff_to_features(gff_input.value)
    
    def load_file(event):
        features = utils.gff_to_features(gff_input.value)
        update(event)
        
    def find_gene(event):
        gene = gene_input.value         
        df = utils.features_to_dataframe(features).fillna('-')
        found = df[df.gene.str.contains(gene)].iloc[0]        
        loc = (found.start-200,found.end+200)
        slider.value = loc
        #feature_pane.object = view_features(features,preview=False,x_range=loc, plot_width=900)
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
    
    top = pn.Row(gff_input,loc_input,gene_input,next_btn)
    main = pn.Column(feature_pane, sizing_mode='stretch_width')
    app = pn.Column(top,slider,main)
    return app


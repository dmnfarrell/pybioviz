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
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool, NumeralTickFormatter, Arrow, NormalHead
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn
from . import utils

def test_plot(rows=20, cols=100):
    """Bokeh random colors plot"""

    x = np.arange(1, cols)
    y = np.arange(0, rows, 1)    
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    
    colors = utils.random_colors(len(gx))
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

def plot_coverage(df, plot_width=800):
    """Plot a bam coverage dataframe returned from get_coverage"""
    
    if df is None:
        return figure(plot_width=plot_width,plot_height=60,tools="")
    source = ColumnDataSource(df)
    x_range=(df.pos.min(),df.pos.max())
    top = df.coverage.max()
    p = figure(title=None, plot_width=plot_width, plot_height=60,
               x_range=x_range, y_range=(0,top), tools="xwheel_zoom",
               min_border=0, toolbar_location='right')
    rects = Rect(x="pos", y=0, width=1, height="coverage", fill_color="gray", fill_alpha=0.3)
    p.grid.visible = False
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.xaxis.visible = False
    return p

def plot_sequence_alignment(aln, fontsize="8pt", plot_width=800):
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

def plot_features(features, preview=True, x_range=None, fontsize="8pt", plot_width=800, plot_height=150):
    """Bokeh sequence alignment view"""
    
    df = utils.features_to_dataframe(features)#, cds=True)
    df = df[df.type!='region']    
    df['length'] = df.end-df.start
    df['level'] = 1
    df['color'] = utils.random_colors(len(df)) #'green'
    df['x'] = df.start+df.length/2
    df['end_x'] = df.start+50
    df['y'] = [random.choice(range(1,6)) for i in range(len(df))]
    
    #print (df[:3])
    text = df.gene
    S = df.start.min()
    N = df.end.max()+10        
    x = list(df.start+df.length/2)
    h = 20

    source = ColumnDataSource(df)    
        
    viewlen=5000
    if x_range == None:
        x_range = (0,viewlen)
    else:
        #viewlen = int(x_range[1])-int(x_range[0])
        print (x_range,viewlen)    

    hover = HoverTool(
        tooltips=[            
            ("gene", "@gene"),     
            ("locus_tag", "@locus_tag"),
            ("protein_id", "@protein_id"), 
            ("length", "@length"),             
        ],
        #names=['rects']
    )  
    tools=[hover,"xpan, xwheel_zoom, save"]
    
    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height, x_range=x_range,
                y_range=(0,10), tools=tools, min_border=0, toolbar_location='right')#, lod_factor=1)
    if viewlen<30000:
        tags = Text(x="x", y="y", y_offset=-10, text="gene", text_align='center',text_color="black", 
                     text_font="monospace",text_font_size=fontsize, name="genetext")
    rects = Rect(x="x", y="y", width="length", height=.4, fill_color="color", fill_alpha=0.4, name='rects')
    arr = Arrow(source=source, x_start="start", x_end="end_x", y_start="y", y_end="y", 
                line_color="black", name='arrows', end=NormalHead(size=10))    
    p1.add_glyph(source, rects)
    p1.add_glyph(source, tags)
    p1.add_layout(arr)
    
    p1.grid.visible = False
    p1.yaxis.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p1.toolbar.logo = None
    p1.xaxis.formatter = NumeralTickFormatter(format="(0,0)")
    
    if preview == True:
        #entire sequence view (no text, with zoom)
        p = figure(title=None, plot_width=plot_width, plot_height=100, x_range=x_range, y_range=(-2,2), tools=tools, 
                        min_border=0, toolbar_location='below')
        rects = Rect(x="x", y="strand", width="length", height=.4, fill_color="colors", line_color='black', fill_alpha=0.6)
        p.add_glyph(source, rects)
        p.yaxis.visible = False
        p.grid.visible = False

        jscode="""    
        var start = cb_obj.value;    
        x_range.setv({"start": start, "end": start+l})   
        """
        callback = CustomJS(
            args=dict(x_range=p1.x_range,l=viewlen), code=jscode)
        slider = Slider (start=1, end=N, value=1, step=100)
        slider.js_on_change('value', callback)
        p = gridplot([[p],[slider],[p1]], toolbar_location='below')
    else:
        p = p1
    return p
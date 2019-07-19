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
import pandas as pd
from . import utils
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn

def view_sequence_alignment(aln, fontsize="8pt"):
    """Bokeh sequence alignment view"""
    
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    #ids=range(len(seqs))
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)    
    cons = get_cons(aln)
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
    p = figure(title=None, plot_width=800, plot_height=50, x_range=x_range, y_range=(0,S), tools=tools, 
                    min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    #p.xaxis.visible = False
    p.yaxis.visible = False
    p.grid.visible = False  
    
    #sequence text view, zoom fixed
    p1 = figure(title=None, plot_width=800, plot_height=plot_height, x_range=view_range, y_range=ids, tools="xpan,reset", 
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
    
    p3 = figure(title=None, plot_width=800, plot_height=30, x_range=p1.x_range, y_range=(Range1d(min(cons),.5)), tools="xpan")
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

def view_features(features, preview=True, fontsize="8pt", plot_width=800, plot_height=100):
    """Bokeh sequence alignment view"""
    
    df = utils.features_to_dataframe(features)#, cds=True)
    df = df[df.type!='region']
    #df['gene'] = df.gene.apply(lambda x: x.locus_tag)
    df['length'] = df.end-df.start
    df['level'] = 1
    df['color'] = 'green'
    df['x'] = df.start+df.length/2
    
    #print (df[:3])
    text = df.gene
    S = df.start.min()
    N = df.end.max()+10
        
    x = list(df.start+df.length/2)
    widths = df.end-df.start
    #print (x,widths)
    h = 20

    source = ColumnDataSource(df)
    x_range = Range1d(S,N, bounds='auto')
        
    viewlen=3000
    view_range = (0,viewlen)
    #tools="xpan, xwheel_zoom, reset, save"

    hover = HoverTool(
        tooltips=[            
            ("gene", "@gene"),     
            ("locus_tag", "@locus_tag"),
            ("protein_id", "@protein_id"), 
            ("length", "@length"),             
        ]
    )  
    tools=[hover,"xpan, xwheel_zoom, save"]
    
    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height, x_range=view_range,
                y_range=(-2,2), tools=tools, min_border=0, toolbar_location='right')#, lod_factor=1)          
    glyph = Text(x="x", y="strand", y_offset=5, text="gene", text_align='center',text_color="black", 
                 text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="strand", width="length", height=.5, fill_color="color", fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)
  
    p1.grid.visible = False
    p1.yaxis.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p1.toolbar.logo = None
    
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

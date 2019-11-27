#!/usr/bin/env python
"""
    Implements viewers/panel apps for pybioviz
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

import os,sys,io
import numpy as np
import pandas as pd
from . import utils, plotters
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, Plot, LinearAxis, Grid, Range1d,CustomJS, Slider, HoverTool)
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column
import panel as pn
import panel.widgets as pnw

def test_app():
    """Test dashboard"""
    
    def refresh(event):
        plot1.object = plotters.test_plot(cols=col_sl.value,rows=row_sl.value, plot_width=600)
        plot2.object = plotters.dummy_plot(rows=row_sl.value, plot_width=600)
        return
    from . import __version__
    title = pn.pane.Markdown('# pybioviz v%s test plots' %__version__)
    plot1 = pn.pane.Bokeh()   
    plot2 = pn.pane.Bokeh()
    col_sl = pnw.IntSlider(name='cols',value=30,start=5,end=200,step=1)
    col_sl.param.watch(refresh, 'value')
    row_sl = pnw.IntSlider(name='rows',value=10,start=5,end=100,step=1)
    row_sl.param.watch(refresh, 'value')
    col_sl.param.trigger('value')    
    app = pn.Column(title,col_sl,row_sl,plot1,plot2)
    return app

def sequence_alignment_viewer(filename=None):
    """Sequence alignment viewer"""
    
    title = pn.pane.Markdown('### Sequence aligner: %s' %filename)
    aln_btn = pnw.Button(name='align',width=100,button_type='primary')
    file_input = pnw.FileInput(name='load file',width=100,accept='.fa,.fasta,.faa')
    aligner_sel = pnw.Select(name='aligner',value='muscle',options=['muscle','clustal','maaft'],width=100)
    highlight_sel = pnw.Select(name='highlight mode',value='default',options=['default',''],width=100)
    seq_pane = pn.pane.HTML(name='sequences',height=200,css_classes=['scrollingArea'])
    bokeh_pane = pn.pane.Bokeh()

    def update_file(event):
        nonlocal seqtext
        seqtext = file_input.value.decode('utf-8')
        title.object = file_input.filename
        #print(file_input.filename)
        #sequences = SeqIO.parse(filename,format='fasta')
        #s = '<p>'.join([rec.format("fasta") for rec in sequences])
        #seq_pane.object = '<div class=monospace>'+s+'</div>'
        return

    def align(event):
        #this function does the alignment
        nonlocal seqtext        
        if seqtext is not None:
            sequences = SeqIO.parse(io.StringIO(seqtext),format='fasta')
        elif filename is not None:    
            sequences = SeqIO.parse(filename, format='fasta')
        else:      
            return
        sequences = list(sequences)
        #print (sequences)
        aligner = aligner_sel.value
        if aligner == 'muscle':
            aln = utils.muscle_alignment(sequences)    
        elif aligner == 'clustal':
            aln = utils.clustal_alignment(sequences)
        elif aligner == 'mafft':
            aln = utils.mafft_alignment(sequences)
        if aln is None:
            bokeh_pane.object = plotters.plot_empty('aligner not found',900)
        else:
            #the bokeh pane is then updated with the new figure
            bokeh_pane.object = plotters.plot_sequence_alignment(aln)  
        return 
   
    seqtext = None
    file_input.param.watch(update_file,'value')
    aln_btn.param.watch(align, 'clicks')
    aln_btn.param.trigger('clicks')
    side = pn.Column(aln_btn,file_input,aligner_sel,highlight_sel,seq_pane,css_classes=['form'],width=200,margin=20)
    app = pn.Column(title,pn.Row(side, bokeh_pane), sizing_mode='stretch_width',width_policy='max',margin=20)
    return app

def view_features(features=None):
    """Genome feature viewer app"""
    
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



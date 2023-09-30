[![PyPI version shields.io](https://img.shields.io/pypi/v/pybioviz.svg)](https://pypi.python.org/pypi/pybioviz/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/pybioviz/badge/?version=latest)](https://pybioviz.readthedocs.io/en/latest/?badge=latest)

# pybioviz

<img align="right" src=img/logo.svg width=150px>

Bioinformatics visualization tools with PyViz Panel and Bokeh. This is a demonstration of bioinformatic dashbaords with panel and bokeh. These could be re-used inside Jupyter notebooks as part of bioinformatic workflows or deployed as local web apps. Not under active development so mainly proof of concept. Use code as needed.

Tools implemented:

* sequence alignment viewer
* genome feature viewer
* bam alignment viewer

<img src=https://github.com/dmnfarrell/pybioviz/raw/master/doc/source/sequence_align_plot.png width=600px>

## Installation

```
pip install -e git+https://github.com/dmnfarrell/pybioviz.git#egg=pybioviz
```

If using JupyterLab you may need to run:
```
jupyter labextension install @bokeh/jupyter_bokeh
jupyter labextension install @pyviz/jupyterlab_pyviz
```

## Notebook

Try the basics.ipynb notebook in the notebook folder to see how it works.

## Usage 

```python
from bokeh.io import show, output_notebook
output_notebook()
import panel as pn
pn.extension()
from pybioviz import dashboards, plotters, utils
```

### Sequence alignment

```python
from Bio import AlignIO, SeqIO

aln_file = os.path.join(utils.datadir, 'test.aln')
aln = AlignIO.read(aln_file,'clustal')
seqview = plotters.plot_sequence_alignment(aln)
main = pn.Column(m,seqview)
main
```

### Plot bam coverage

```python
bam_file='wt_mbovis.bam'
chr='NC_002945.4'
start=1000
end=3000
df = utils.get_coverage(bam_file, chr, start, end)
p = plotters.plot_coverage(df)
```

## Dashboards 

Dashboards are small interctive apps that are created by using the plotting tools detailed above, combined with widgets that allow user interactivity. They can be created inside a notebook or launched as standalone apps in a web browser using a command in the terminal.

### Sequence alignment dashboard

This dashboard is for used for viewing and manipulating multiple sequence alignments. You can load files from local directories and re-align, zoom in and out of the sequence. Other features to add and remove sequences have yet to be added.

To use this inside a notebook:

```python
app = dashboards.sequence_alignment_viewer('file.fa')
```

<img src=https://github.com/dmnfarrell/pybioviz/raw/master/doc/source/sequence_aligner_dashboard.gif width=600px>


### Genomic feature viewer

```python
app = dashboards.genome_features_viewer('Mbovis_AF212297.gff','Mbovis-AF212297.fa')
```

<img src=https://github.com/dmnfarrell/pybioviz/raw/master/doc/source/genome_features_dashboard.gif  width=600px>

### Links

* https://panel.pyviz.org/
* https://docs.bokeh.org/en/latest/

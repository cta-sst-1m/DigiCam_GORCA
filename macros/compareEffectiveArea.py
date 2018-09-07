#!/usr/bin/env python

import sys, os, glob
from plotAreaEfficiencies import *

files = os.path.normpath(os.path.expandvars(os.path.expanduser(sys.argv[1])))
outname = sys.argv[2] if len(sys.argv)>2 else None
if files[-4:]!=".npz" and files.find('*') < 0:
    files = os.path.join(files,'*')
if files.find('*') >= 0:
    files = sorted(glob.glob(files))
if not isinstance(files, list):
    files = [files]

colors = ['b','r','k','g','m','c','y']
markers = ['o','s','h','H','*']

fig_and_plot = None
for i,fname in enumerate(files):
    with np.load(fname) as f:
        E, Aeff, AeffError = f['E'], f['Aeff'], f['AeffError']
    label = os.path.basename(fname)
    label = " ".join(label.split("_")[:-1])
    color, marker = colors[i%len(colors)], markers[i/len(colors)]
    latest = fname==files[-1]
    fig_and_plot = EffectiveArea(E, Aeff, AeffError, outname=outname, fig_and_plot=fig_and_plot, latest=latest, color=color, marker=marker, label=label)
raw_input(" => Press enter to close... <=")


    

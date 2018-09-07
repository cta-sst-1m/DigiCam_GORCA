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

for i,fname in enumerate(files):
    with np.load(fname) as f:
        Rs, trigg, triggErrors, Rerrors, N, (thisLge,nextLge) = f['Rs'], f['trigg'], f['triggErrors'], f['Rerrors'], f['N'], f['Lge']
    Refficiency(Rs, trigg, triggErrors, Rerrors, N, thisLge, nextLge, outname)
raw_input(" => Press enter to close... <=")


    

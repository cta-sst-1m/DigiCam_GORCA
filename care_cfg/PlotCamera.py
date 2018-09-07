#!/usr/bin/env python

import sys, os, docopt, numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import camera_descriptor

version = "Version: 0.1 (I. Al Samarai, imen.alsamarai@gmail.com)"
doc = """

Tool to plot the camera pixels and compare configuration files on demand. Eg. python PlotCamera.py 'camera_config1.dat camera_config2.dat camera_config....dat'

Usage:
  %s (-h | <input>) [--identif ]  

Argument:
  <input>                 File containing pixel coordinates 
Options:
  -h --help               Print this help
  --identif               Plot pixel ID
""" % os.path.basename(sys.argv[0])
args = docopt.docopt(doc)


# Absolute path
path = "./"
files = args['<input>'].split()

fig = plt.figure()


def PlotPix(filename, name = None ):

    print filename
    data_pix = numpy.genfromtxt(filename, usecols=(4, 5, 6), skip_header=4, max_rows=1296, unpack=True)
    #data_patch = numpy.genfromtxt(filename, usecols=(4, 5, 6, 7), skip_header=3907, max_rows=433, unpack=True)


    ax = fig.add_subplot(111, aspect='equal', title = '')
    ax.set_xlim([min(data_pix[1, :]) - 25, max(data_pix[1, :]) + 25])
    ax.set_ylim([min(data_pix[2, :]) - 25, max(data_pix[2, :]) + 25])
    for p in [patches.RegularPolygon((data_pix[1, i], data_pix[2, i]), 6, 13.4, fill=False, edgecolor="green")
        for i in range(len(data_pix[0, :]))
    ]:
        ax.add_patch(p)


    if args["--identif"]  :
        for i in range(len(data_pix[0, :])):
            plt.text(data_pix[1, i], data_pix[2, i], '%d' % data_pix[0, i],
                                                    horizontalalignment='center', verticalalignment='center',
                                                    fontsize=5)


for num, file in enumerate(files):
    PlotPix(path+file, name = file)

plt.show()

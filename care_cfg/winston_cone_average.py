#!/usr/bin/env python

import os, sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as pyp
import docopt

doc="""
This program plot the Winstone Cone efficiency as a function of the incident angle
and calculate its integral average (with and without fresnel losses)

Usage:
  %s (-h | <input>) [<savename> -f FRESNEL_LOSSES -l LIGHTGUIDE] [-s...]

Argument:
  <input>                      File name of the lookup table with the Winsont Cone
                               efficiency as a function of the angle
  <savename>                   Save name of the plot. Is not extension is given, '.eps' is used

Options:
  -f --fresnel FRESNEL_LOSSES  Fraction of the fresnel losses. It must be a number [default: 0.]
  -h --help                    Print this help
  -l --lightguide LIGHTGUIDE   Set the lightguide weights file
  -s --show                    Do plots. If <savename> is set, they will be saved (with .eps extension).
                               If called more than once the weighted plot will be shown, renormalized to the maximum
""" % os.path.basename(sys.argv[0])

args = docopt.docopt(doc)

rDish = 200. #cm; diameter of the outside circle is 4m => radius 200cm
hCamera = 44. #cm; camera height = 90cm, but active is 88cm => h of the camera hexagon = 44cm
dDishCamera = 560. #cm; it is the focal lenght!
try:
    fresnel_losses = float(args["--fresnel"])
except ValueError:
    print "\n*** --fresnel must be a (float) number! ***"
    print doc
    raise SystemExit
fresnel_accept = 1-fresnel_losses

maxAngle = round(np.arctan((rDish+hCamera*2/np.sqrt(3))/dDishCamera)*180/np.pi,2) # 22.3
print "The maximum accepted angle is %s\n" % maxAngle

ang,eff = np.loadtxt(args["<input>"],unpack=True)
w_ang,weights = np.loadtxt(args["--lightguide"],unpack=True) if args["--lightguide"] is not None else ([], [])

# x values
x = np.append(ang[ang<maxAngle],maxAngle)
# Interpolation for the y value corresponding to x=maxAngle
dim = len(x)
frac = (x[-1]-x[-2])/(ang[dim-1]-x[-2])
minEff = frac*(eff[dim-1]-eff[dim-2])+eff[dim-2]
#y values
y = np.append(eff[ang<maxAngle],minEff)
# dx
xdiffs = np.diff(x)
# central bin values and dx
xcentral = (x[:dim-1]+x[1:dim])/2.
ycentral = (y[:dim-1]+y[1:dim])/2.
# ...weighting
if len(w_ang)>0:
    w_ang = (w_ang[:-1]+w_ang[1:])/2.
    weights = (weights[:-1]+weights[1:])/2.
    w_ang, weights = w_ang[w_ang>0.], weights[w_ang>0.]
    #select = w_ang>0.
    #w_ang = w_ang[select]
    #weights = weights[select]
    w_x = np.append(w_ang[w_ang<maxAngle],maxAngle)
    w_dim = len(w_x)
    w_frac = (w_x[-1]-w_x[-2])/(w_ang[w_dim-1]-w_x[-2])
    endw = w_frac*(weights[w_dim-1]-weights[w_dim-2])+weights[w_dim-2]
    w_y = np.append(weights[w_ang<maxAngle],endw)
    # Find knots, B-spline coefficients, and spline degree and set the interpolated trasmittance
    tck = interpolate.splrep(w_x,w_y)
    ww = interpolate.splev(xcentral,tck)
    xdiffs *= ww # 
    ycentral *= ww # weighting the y values
    if args['--show']:
        fw,pw = pyp.subplots()
        fw.suptitle(r"angle probaility (sum = 1)", weight='bold', size='xx-large')
        pw.plot(xcentral, ww/np.sum(ww), 'r-')
        pw.set_xlabel(r"Angle [$\mathbf{^o\!}$]", weight='bold')
        pw.set_ylabel(r"a.u.", weight='bold')
        fw.show()
    
# winstone cone (weighted) average
effAverage = np.sum(ycentral*np.diff(x))/np.sum(xdiffs)
print "The integral average of the efficiency is %1.2f" % effAverage
print "  ...with a fresnel losses of %s%%, the Winston Cone average efficiency is %1.2f\n" % (100*fresnel_losses,effAverage*fresnel_accept)

# Plot efficiency
if args['--show']:
    f,p = pyp.subplots()
    title = r"Funnel efficiency"# (average = %1.2f)" % effAverage
    f.suptitle(title, weight='bold', size='xx-large')
    p.plot((ang[:-1]+ang[1:])/2.,(eff[:-1]+eff[1:])/2., 'r-', mec='b', label=r"not weighted")
    #p.plot(ang,eff,"ob",mec='b',label="original values")
    #p.plot(x,y,'or',mec='r',ms=3,label="cut values")
    p.set_xlabel(r"Angle [$\mathbf{^o\!}$]", weight='bold')
    p.set_ylabel(r"Efficiency", weight='bold')
    p.vlines(maxAngle,0,1,colors='g',linestyle='-.')
    p.annotate("max angle:\n%s"%maxAngle,xy=(maxAngle,0.7),xytext=(maxAngle+5,0.9),size='x-large',ha="center",
               arrowprops=dict(arrowstyle="fancy",fc="g", ec="none",connectionstyle="angle3,angleA=90,angleB=0"))
    if args['--show']>1 and len(w_ang)>0:
        p.plot(xcentral,ycentral*np.max(y)/np.max(ycentral), 'b--', mec='r', label="weighted")
        p.legend(loc="right",bbox_to_anchor=(1.03,.5), handletextpad=0.)#, title="normalized\nto the same\nmaximum", fontsize='x-large')
    #p.set_xlim(0.,np.max(ang))
    #p.set_ylim(.0,1.)
    f.show()
# save
if args["<savename>"] is not None:
    name = args["<savename>"]
    base,ext = os.path.splitext(name)
    name = base+".eps"
    if args['--show']:
        f.savefig(name)
        if len(w_ang)>0:
            wname = base+'_weights.eps'
            fw.savefig(wname)
raw_input("\n\t==> Press any key to close <==\n")







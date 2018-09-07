#!/usr/bin/env python

import os, sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as pyp
from matplotlib import rc, rcParams
import docopt

doc="""
Set the lookup table of the PDE for CARE configs. If the fresnel trasmittance of the
window is given, this will be used to correct the Q.E. of the PDE

Usage:
  %s (-h | <qe>) [<fresnel> -c CNAME -o OUTNAME -t TELTYPE -vw] [-s...]

Argument:
  <qe>                  File name of the lookup table of the QE
  <fresnel>             File name of the lookup table of the fresnel losses (as a function
                        of wavelengths)

Options:
  -h --help             Print this help
  -c --cherenkov CNAME  File with the cherenkov spectrum data
  -o --outname OUTNAME  Save name of the lookup table. Is not extension is given, '.dat' is used.
                        If not set, the lookup table will printed in the stdout
  -s --show             Do plots. If -o is set, they will be saved (with .eps extension).
                        If called twice (level 2) will show the <fresnel> effect, if defined
                        Set a level 3 (called three or more times) to show the weighted plot,
                        renormalized to the maximum
  -t --teltype TELTYPE  Define the telescope type ID (see CARE). Must be an int! [default: 0]
  -v --verbose          Print the lookup table in stdout
""" % os.path.basename(sys.argv[0])

args = docopt.docopt(doc)
print

rcParams['xtick.labelsize'] = 'xx-large'
rcParams['ytick.labelsize'] = 'xx-large'

# Check on inputs
try:
    teltype = float(args["--teltype"])
except ValueError:
    print "*** --teltype must be an int! ***"
    print doc
    raise SystemExit
name = args["--outname"]

# Load the qe file
wave,qe = np.loadtxt(args["<qe>"],unpack=True)
weights, extra_weight = np.ones_like(wave), 1.
if args["--show"]:
    f1,p1 = pyp.subplots()
    #title = "PDE"
    #f1.suptitle("PDE", size="xx-large",weight='bold')
    p1.set_xlabel(r"$\lambda$ [nm]",fontsize='xx-large')
    p1.xaxis.set_label_coords(.5, -0.06)
    p1.set_ylabel(r"Efficiency",fontsize='xx-large')
    p1.plot(wave,qe, 'b-', label="PDE")
    p1.axvline(300.,color='g',ls='-.')
    p1.axvline(550.,color='g',ls='-.')
# Eventual frensel file with qe adjustment with the window fresnel losses 
if args["<fresnel>"] is not None:
    wave_fresnel,transmittance = np.loadtxt(args["<fresnel>"],unpack=True)
    # Find knots, B-spline coefficients, and spline degree and set the interpolated trasmittance
    tck = interpolate.splrep(wave_fresnel,transmittance)
    limits = wave<=wave_fresnel[-1]
    new_trasm = interpolate.splev(wave[limits],tck)
    if wave[-1] > wave_fresnel[-1]:
        new_trasm = np.append(new_trasm,[0.]*(len(wave)-len(wave[limits])))
    qe *= new_trasm
    if args["--show"]:
        f2,p2 = pyp.subplots()
        f2.suptitle(r"Transmittance (funnels $\times$ window)", size="xx-large",weight='bold')
        p2.set_xlabel(r"$\lambda$ [nm]",fontsize='xx-large')
        p2.xaxis.set_label_coords(.5, -0.06)
        p2.set_ylabel(r"Transmittance",fontsize='xx-large')
        p2.plot(wave_fresnel, transmittance, 'b-')#, wave, new_trasm, 'ro', mec='r')
        #p2.legend(["Fresnel transmittance",r"points for $\lambda$ in QE"])#,
        #          #bbox_to_anchor=(0.5,0.2),loc='lower center')
        p2.axvline(300.,color='g',ls='-.')
        p2.axvline(550.,color='g',ls='-.')
        if args["--show"]>1:
            p1.plot(wave,qe, 'r--', label="with funnel &\nwindow trans.")
        f2.show()

if args["--cherenkov"] is not None:
    wave_cherenkov, flux_cherenkov = np.loadtxt(args["--cherenkov"],unpack=True,usecols=(0,1))
    # Find knots, B-spline coefficients, and spline degree and set the interpolated trasmittance
    tck_ch = interpolate.splrep(wave_cherenkov,flux_cherenkov)
    limits = wave<=wave_cherenkov[-1]
    weights = interpolate.splev(wave[limits],tck_ch)
    if wave[-1] > wave_cherenkov[-1]:
        weights = np.append(weights,[0.]*(len(wave)-len(wave[limits])))
    if np.min(wave)>300.:
        extra_weight = interpolate.splev((np.min(wave)+300.)/2.,tck_ch)
    #weights = new_flux*np.max(qe)/np.max(new_flux)
    if args["--show"]:
        f3,p3 = pyp.subplots()
        f3.suptitle("Cherenkov fluence", size="xx-large",weight='bold')
        p3.set_xlabel(r"$\lambda$ [nm]",fontsize='xx-large')
        p3.xaxis.set_label_coords(.5, -0.06)
        p3.set_ylabel(r"fluence [$\mathbf{ph\cdot\,}$nm$\mathbf{^{-1}\cdot\,}$m$\mathbf{^{-2}}$]",fontsize='xx-large')
        p3.plot(wave_cherenkov, flux_cherenkov, 'b-')#, wave, weights, 'ro', mec='r')
        #p3.legend(["Cherenkov flux",r"points for $\lambda$ in QE"])#,
        #          #bbox_to_anchor=(0.5,0.2),loc='lower center')
        p3.axvline(300.,color='g',ls='-.')
        p3.axvline(550.,color='g',ls='-.')
        if args["--show"]>2:
            p1.plot(wave,qe*weights*np.max(qe)/np.max(qe*weights), 'k:', label="weighted with\nCherenkov flux")
        #title += " $\otimes$ funnel & window transm."
        f3.show()

# Average
select = (wave>=300.)&(wave<=550.)
dx = np.diff(wave[select])
sub_qe = qe[select]
weights = weights[select]
central_w = (weights[:-1]+weights[1:])/2.
average = np.sum((sub_qe[:-1]+sub_qe[1:])*central_w*dx/2.)/(np.sum(dx*central_w)+((np.min(wave)-300.)*extra_weight if np.min(wave)>300. else 0.))
print " --> Q.E. total average = %s\n"%average
if args["--show"]:
    title = "PDE"+((" (average from 300 to 550 nm = %1.2f %%)" % (average*100)) if args["--show"]>2 else "")
    f1.suptitle("%s"%title, size="xx-large",weight='bold')

# Build the text of the qe efficiency
qe_txt = "* QUEFF %d %d\n" % (teltype,len(wave))
for i,w in enumerate(wave):
    qe_txt += "{:>3.2f} {:<1.6f}\n".format(w,qe[i])
if name is None or args['--verbose']:
    print "--> These are the Q.E.:\n"
    print qe_txt
if name is not None:
    name,ext = os.path.splitext(args["--outname"])
    if ext == "":
        ext = ".dat"
    with file(name+ext,'w') as f:
        f.write(qe_txt)
    print "...saved the Q.E. text as %s\n" % (name+ext)
if args["--show"]:
    if args["--show"]>1:
        p1.legend()
    if name is not None:
        f1.savefig(name+".eps")
        print "...saved the Q.E. plot as %s.eps\n" % name
        if args["<fresnel>"] is not None:
            f2.savefig(name+"_trans.eps")
            print "...saved the fresnel trasmittance plot as %s_trans.eps\n" % name
        if args["--cherenkov"] is not None:
            f3.savefig(name+"_ChFlux.eps")
            print "...saved the fresnel trasmittance plot as %s_ChFlux.eps\n" % name
    f1.show()
    raw_input("\n\t==> Press any key to close <==\n")

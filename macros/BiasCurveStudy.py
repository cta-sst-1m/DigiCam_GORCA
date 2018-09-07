#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys

    version = "Version: 0.1 (A. Porcelli, ale_led@yahoo.it)"
    version = "Version: 0.2 (I. Al Samarai)"
    doc = """
Plot NSB trigger rate from the log files in a folder

Usage:
  %(prog)s <input>...  [-g GAIN -o OUTPUT -y YLIM [-y YLIM] ] [-l LABEL...] [-V...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <input>                    Input folder containing the log file (gunzipped!) of the care BiasCurve runs or directly the log file.
                             If multiple <input>, each of them is considered a different rate to be plotted. Max supported inputs:
                             14 (are you sure you want more than 14 inputs? The plot it will too crowded and hard to be read...)

General options:
  -o --output OUTPUT         If defined, this name will be used to save the plot(s)
  -V --verbose               Level of verbose output

Setting options:
  -g --gain GAIN             Set the gain used to convert pe in ADC counts [default: 5.6]
  -l --label LABEL           String used as label. It is 1 to 1 with the <input> list. If not given 'input <N>' will be used, whene <N>
                             is the number of the positional <input> arguments to which the legend is referring (starting from 1).
                             Note: if only an <input> is given, no legend is plotted. 
  -y --ylim YLIM             Minimum and maximum value(s) of the y-axis [default: 100]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from commons import *

if __name__=='__main__':

    gain = float(args['--gain'])
    ylim = [float(yl) for yl in args['--ylim']]
    
    line  = lambda i: ['-','--','-.',':'][i%4]
    mark = lambda i: ['o','s','^','v','<','>','x'][i%7]
    color = lambda i: ['k','b','r','g','y','c','m'][i%7]
    fill = lambda i: color(i) if i//7==0 else 'w'

    hplot, wplot = 6-1.2,8*(1-0.125-0.1)
    dleft, dright, dtop, dbottom = 0.125*8., 8.*0.1, .75, 0.6
    figsize = (dleft+wplot+dright, dbottom+hplot+dtop)

    f, p = pyp.subplots(figsize=figsize)
    f.subplots_adjust(left=dleft/figsize[0],right=(dleft+wplot)/figsize[0],bottom=dbottom/figsize[1],top=(dbottom+hplot)/figsize[1])
    p.set_xlabel("threshold [ADC]")
    p.set_ylabel("Rate [Hz]")
    p.set_yscale('log')
    xmin, xmax = np.inf, -np.inf
    for i,path in enumerate(args['<input>']):
        if i>=14:
            print("Too many input! Skipped above the 14th...")
            break
        files = glob.glob(aspath(*((path,"*.log") if os.path.isdir(path) else (path,))))
        th, rt, er = getNSBRatesFromLogs(*files)
        sel = (rt!=0)&(er!=0)
        th = th[sel]
        if len(th)==0:
            print(" %s: all data are 0"%path)
            continue
        xmin, xmax = min(xmin, th[0]), max(xmax, th[-1])
        p.errorbar(th, rt[sel], er[sel], fmt=line(i), mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                   label=r"%s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
    p.set_xlim(xmin, xmax)
    p.set_ylim(*ylim)
    ppe = p.twiny()
    ppe.set_xlabel("p.e. in the cluster")
    ppe.set_xlim(xmin/gain, xmax/gain)
    if len(args['<input>'])>1:
        p.legend()
    if args['--output'] is not None:
        f.savefig(args['--output'])
    f.show()
    raw_input("\n == Press enter to close ==\n")

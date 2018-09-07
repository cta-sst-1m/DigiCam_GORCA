#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys

    version = "Version: 0.1 (A. Porcelli, ale_led@yahoo.it)"
    version = "Version: 0.2 (I. Al Samarai)"
    doc = """
    
Plot NSB trigger rate from the log files in a folder

Usage:
  %(prog)s <bcinput>...  [-g GAIN -o OUTPUT -u UNIT -X XMAX -x XMIN -y YLIM [-y YLIM] ] [-l LABEL...] [-p PINPUT...] [-V...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <bcinput>                  Input folder containing the log file (gunzipped!) of the care BiasCurve runs or directly the log file.
                             If multiple <bcinput>, each of them is considered a different rate to be plotted. Max supported inputs:
                             7 (are you sure you want more than 7 inputs? The plot it will too crowded and hard to be read...), beyond
                             that, color and markers will repeat and the result might be confusing...

General options:
  -o --output OUTPUT         If defined, this name will be used to save the plot(s)
  -V --verbose               Level of verbose output

Setting options:
  -g --gain GAIN             Set the gain used to convert pe in ADC counts [default: 5.6]
  -l --label LABEL           String used as label. It is 1 to 1 with the <input> list. If not given 'input <N>' will be used, whene <N>
                             is the number of the positional <input> arguments to which the legend is referring (starting from 1).
                             Note: if only an <input> is given, no legend is plotted.
  -p --proton PINPUT         Input of the proton rate input. It is 1 to 1 with the <bcinput> list. If not given nothing will be taken into account.
                             Note: if only an <input> is given, no legend is plotted.
  -u --unit UNIT             Unit of the x-axis [default: ADC]
  -x --xmin XMIN             Minimum value of the x-axis
  -X --xmax XMAX             Maximum value of the x-axis
  -y --ylim YLIM             Minimum and maximum value(s) of the y-axis [default: 100]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from commons import *

def findThreshold(bth, brt, pth, prt, unit='ADC'):
    if len(bth)<len(pth):
        tck = interpolate.splrep(bth, prt)
        brt = interpolate.splev(pth,tck)
        bth = pth
    elif len(bth)>len(pth) or not all(pth==bth):
        tck = interpolate.splrep(pth, prt)
        prt = interpolate.splev(bth,tck)
        pth = bth
    i0, i1 = np.where(brt>=prt)[0][-1], np.where(brt<prt)[0][0]
    pth0, pth1, bth0, bth1 = pth[i0], pth[i1], bth[i0], bth[i1]
    plgr0, plgr1, blgr0, blgr1 = np.log10(prt[i0]), np.log10(prt[i1]), np.log10(brt[i0]), np.log10(brt[i1])
    p_m, b_m = (plgr1-plgr0)/(pth1-pth0), (blgr1-blgr0)/(bth1-bth0)
    safe_th = (pth0*p_m - bth0*b_m + blgr0 - plgr0) / (p_m - b_m)
    x = np.linspace(pth0, pth1, 1000)
    y = plgr0 + (x-pth0)*p_m
    print ">>> SAFE THRESHOLD: %2.1f %s <<<" % (safe_th,unit)
    return safe_th, 10**y[x>=safe_th][0]


if __name__=='__main__':

    gain = float(args['--gain'])
    ylim = [float(yl) for yl in args['--ylim']]
    xmin = np.inf if args['--xmin'] is None else float(args['--xmin'])
    xmax = -np.inf if args['--xmax'] is None else float(args['--xmax'])
    unit = args['--unit']
    
    line  = lambda i: ['-','--','-.',':'][i]#[i%4]
    mark = lambda i: ['o','s','^','v','<','>','x'][i%7]
    color = lambda i: ['b','r','k','g','y','c','m'][i%7]
    fill = lambda i: color(i) if i//7==0 else 'w'

    hplot, wplot = 6-1.2,8*(1-0.125-0.1)
    dleft, dright, dtop, dbottom = 0.125*8., 8.*0.1, .75, 0.6
    figsize = (dleft+wplot+dright, dbottom+hplot+dtop)

    f, p = pyp.subplots(figsize=figsize)
    f.subplots_adjust(left=dleft/figsize[0],right=(dleft+wplot)/figsize[0],bottom=dbottom/figsize[1],top=(dbottom+hplot)/figsize[1])
    p.set_xlabel("threshold [%s]"%unit)
    p.xaxis.set_label_coords(.5, -0.06)
    p.set_ylabel("Rate [Hz]")
    p.set_yscale('log')
    SafeTh, YTH = [], []
    for i,path in enumerate(args['<bcinput>']):
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
        rtsel, ersel = rt[sel], er[sel]
        p.errorbar(th, rtsel, ersel, fmt=line(0), mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                   label=r"NSB: %s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
        if len(args['--proton'])>i:
            pth, prt, per = np.loadtxt(args['--proton'][i], unpack=True, skiprows=1)
            nans = (np.isnan(prt))&(np.isnan(per))
            pth, prt, per = pth[~nans], prt[~nans], per[~nans]
            p.errorbar(pth, prt, per, fmt=line(1), mec=color(i), mfc='w', marker=mark(i), color=color(i), lw=1.,
                       label=r"proton: %s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
            xmin, xmax = min(xmin, pth[0]), max(xmax, pth[-1])
            safe_th, y_th = findThreshold(th, rtsel, pth, prt, unit=unit)
            SafeTh += [safe_th]
            YTH += [y_th]
            p.axvline(safe_th, color=color(i), ls='-.', lw=1.)
    p.set_xlim(xmin, xmax)
    p.set_ylim(*ylim)
    ymin, ymax = p.get_ylim()
    ppe = p.twiny()
    ppe.set_xlabel("p.e. in the cluster")
    ppe.set_xlim(xmin/gain, xmax/gain)
    for i in range(len(args['--proton'])):
        p.annotate("safe th\n%2.1f %s"%(SafeTh[i],unit), xy=(SafeTh[i], YTH[i]), xytext=(SafeTh[i]*0.99,10**(np.log10(ymax)*0.99-len(args['--proton'])+i*1.5)),
                   size='x-large', ha="right", va='top', xycoords='data', color=color(i))
    if len(args['<bcinput>'])>1 or len(args['--proton'])>0:
        p.legend()
    if args['--output'] is not None:
        f.savefig(args['--output'])
    f.show()
    raw_input("\n == Press enter to close ==\n")

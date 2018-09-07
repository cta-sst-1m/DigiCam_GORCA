#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys

    version = "Version: 0.3 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Make the proton rate

Usage:
  %(prog)s <input>... [-b BIAS -c CONEVIEW -e EWIDTH -f -g GAIN -l LABEL -M MAXTH -m MINTH -n -o OUTPUT -r RWIDTH -s -t FILETYPE -u UNIT] [-V...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <input>                    Input directory name(s). The root files must be located in "<input>/Root/<file>.root".
                             All <file>s will be processed. The <input> must contain the pattern '_th<N>' where
                             <N> is a number. This number will be considered the threshold level to associate the
                             trigger rate calculated from all the <file>s in the folder. In <input> might be a '*'
                             to get multiple directories (i.e. TChain compatible). For example (typical use) the
                             "/my/output/foloder/MySims_th*" or /my/output/foloder/MySims_th* input will get all
                             the simulated trigger thresholds

General options:
  -s --savenpz               Save npz file of the variables/fit parameters, IF <output> is defined
  -V --verbose               Level of verbose output

Plot options:
  -e --energy-width EWIDTH   Width of the log10(energy) bins [default: 0.25] 
  -f --fit                   Use fitting functions
  -g --gain GAIN             <UNIT>/p.e. gain (to traslate x-axis in equivalent p.e.) [default: 5.6]
  -l --label LABEL           Label for the Effective Area plot
  -M --maxth MAXTH           Maximum threshold considered [default: np.inf]
  -m --minth MINTH           Minimum threshold considered [default: 0.]
  -n --nevents               Number of events are shown in the R-efficiency plots
  -o --output OUTPUT         If defined, this name will be used to save the plot(s) and, eventually npz file
  -r --radius-width RWIDTH   Width of the radius bins in meters [default: 100]
  -t --filetype FILETYPE     Type of file (npz and root are for now available; it is the file extension!)
                             [default: npz]
  -u --unit UNIT             x-axis unit [default: ADC]

Trigger rate:
  -b --biascurve BIAS        Txt (log) file or a folder with .log files where to find the biascurve data point.
                             If these data are given, the trigger threshold (in ADC counts) will be printed.
                             Enquoting BIAS, it is possible to use wildcards: multiple files will be merged,
                             but they must have the same threshold binning...
  -c --coneview CONEVIEW     Calculate the proton trigger rate with the as efficiency folded p-flux times
                             2pi*[1-cos(CONEVIEW/2.)]. CONEVIEW is the angulare aperture in deg [default: 11]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from TriggerInfo import *
import re

def plotTriggerRates(p_th, p_rate, p_rate_err, savename=None, show=True, fig_and_plot=None, biascurve=None, savenpz=False, unit='ADC', gain=5.6, **errorbar_kwd):
    f, p = pyp.subplots() if fig_and_plot is None else fig_and_plot
    marker, color, labels = errorbar_kwd.pop('marker','o'), errorbar_kwd.pop('color','b'), [errorbar_kwd.pop('label',None)]
    pplot = p.errorbar(p_th, p_rate, p_rate_err, marker=marker,fmt='b-',mec=color,mfc=color, **errorbar_kwd)
    handles = [pplot]
    p.set_yscale('log')
    p.set_xlabel(r"trigger threshold [%s]"%unit,fontsize='xx-large')
    p.set_ylabel(r"trigger rate [Hz]",fontsize='xx-large')
    if biascurve is not None:
        biascurve = glob.glob(aspath(*((biascurve,"*.log") if os.path.isdir(biascurve) else (biascurve,))))
        b_th, b_rate, b_rate_err = getNSBRatesFromLogs(*biascurve)#np.loadtxt(biascurve, unpack=True)
        bplot = p.errorbar(b_th, b_rate, b_rate_err, marker=marker,fmt='b--',mec=color,mfc='w', **errorbar_kwd)
        handles += [bplot]
        if labels[0] is None:
            labels[0] = 'protons'
            labels += ['NSB']
        else:
            labels += [labels[0]+' (NSB)']
            labels[0] += ' (protons)'
        if len(b_th)<len(p_th):
            tck = interpolate.splrep(b_th, p_rate)
            b_rate = interpolate.splev(p_th,tck)
            b_th = p_th
        elif len(b_th)>len(p_th) or not all(p_th==b_th):
            tck = interpolate.splrep(p_th, p_rate)
            p_rate = interpolate.splev(b_th,tck)
            p_th = b_th
        i0, i1 = np.where(b_rate>=p_rate)[0][-1], np.where(b_rate<p_rate)[0][0]
        p_th0, p_th1, b_th0, b_th1 = p_th[i0], p_th[i1], b_th[i0], b_th[i1]
        p_lgr0, p_lgr1, b_lgr0, b_lgr1 = np.log10(p_rate[i0]), np.log10(p_rate[i1]), np.log10(b_rate[i0]), np.log10(b_rate[i1])
        p_m, b_m = (p_lgr1-p_lgr0)/(p_th1-p_th0), (b_lgr1-b_lgr0)/(b_th1-b_th0)
        safe_th = (p_th0*p_m - b_th0*b_m + b_lgr0 - p_lgr0) / (p_m - b_m)
        print ">>> SAFE THRESHOLD: %2.1f %s <<<" % (safe_th,unit)
        p.axvline(safe_th,color='k',ls='-.', lw=1.)
        brate_sel = b_rate[b_rate>0.]
        y_pos_note = np.sqrt(brate_sel[0]*brate_sel[-1])#10**((np.log10(b_rate[0])+np.log10(b_rate[-1]))/2.)
        p.annotate("safe th\n%2.1f %s"%(safe_th,unit), xy=(safe_th,y_pos_note), xytext=(safe_th,y_pos_note), size='x-large', ha="left", xycoords='data')
    ppe = p.twiny()
    ppe.set_xlabel("p.e. in the cluster")
    ppe.set_xlim(xmin/gain, xmax/gain)
    if labels[0] is not None and show:
        p.legend(handles, labels, bbox_to_anchor=(.5,1.), loc="upper left")
    if savename is not None:
        name = basename(savename)
        savename += '_protonTriggerRate'
        f.suptitle("%s"%(name.replace("_"," ")), size="xx-large", weight='bold')
        if savenpz:
            np.savez("%s.npz"%savename, th=th, rate=rate, rate_err=rate_err)
    if show:
        if savename is not None:
            f.savefig(savename+'.eps')
        f.show()
    return f, p

if __name__=='__main__':
    folders = sorted([ifile for path in args['<input>'] for ifile in glob.glob(aspath(path))])
    savename = args['--output']
    if savename is not None:
        savename = aspath(savename)
        savedir = dirname(savename)
        if savedir != '':
            mkdirs(savedir)
    verbose = args['--verbose']
    lge_width = float(args['--energy-width'])
    r_width = float(args['--radius-width']) # m
    solid_angle = pi2*(1-np.cos(deg2rad*float(args['--coneview'])/2.))
    fittings = args['--fit']
    savenpz = args['--savenpz']
    gain = float(args['--gain'])
    unit = args['--unit']
    minth, maxth = eval(args['--minth']), eval(args['--maxth'])
    # calculate
    th, rate, rate_err = [], [], []
    for folder in folders:
        triggTh = re.findall('(?<=_th)\d+',folder)[-1]
        if len(triggTh)==0:
            continue
        triggTh = float(triggTh)
        if triggTh<minth or triggTh>maxth:
            continue
        th += [triggTh]
        print(" Threshold: %s"%th[-1])
        if savename is None:
            sname = None
        else:
            outfolder = aspath((savename+'_th%s'%triggTh))
            mkdirs(outfolder)
            sname = os.path.join(outfolder, 'crStudy')
        # Get data
        lge, r, r_triggered, _ = build_effArea_variables(os.path.join(folder,'*.%s'%args['--filetype']), sname if savenpz else None, verbose+1)
        infs = np.isinf(lge)
        E, Aeff, AeffError = calc_effectiveArea(lge[~infs], r[~infs], r_triggered[~infs], lge_width, r_width, None, args['--nevents'], sname, fittings, verbose, savenpz)
        del lge, r, r_triggered
        if verbose:
            print "...plot the effective area"
        prate, prate_err = plot_effArea(E, Aeff, AeffError, proton_flux, solid_angle, sname, factor=1.5, integral=True, label=args['--label'], savenpz=savenpz)
        del E, Aeff, AeffError
        if rate is not None:
            rate += [prate]
            rate_err += [prate_err]
        if verbose>1:
            raw_input(" press enter to continue...")
        pyp.close('all')
    th, rate, rate_err = np.array(th), np.array(rate), np.array(rate_err)
    thsel = (th>=minth)&(th<=maxth)
    #th, rate, rate_err = th[thsel], rate[thsel], rate_err[thsel]
    lines = "{:^9s} {:^12s} {:^12s}\n".format('Threshold', 'Rate', 'Rate Uncert.')
    for i,r in enumerate(rate):
        if not thsel[i]:
            continue
        lines += "{:^ 9.2f} {:^ 12.3f} {:^ 12.3f}\n".format(th[i], r, rate_err[i])#'nan' if rate_err[i] is None else rate_err[i])
    with file(savename+'_data.txt','w') as F:
        F.write(lines)
    if verbose:
        print lines
    if has_mpl:
        fig, plot = plotTriggerRates(th, rate, rate_err, biascurve=args['--biascurve'], savename=savename, show=True, label=args['--label'],
                                     gain=gain, unit=unit, savenpz=savenpz)
    if verbose:
        raw_input(" => Press enter to close... <=")
    else:
        time.sleep(2)
    pyp.close('all')

    
    

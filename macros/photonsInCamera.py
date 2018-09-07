#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys

    version = "Version: 0.2 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Make the proton rate

Usage:
  %(prog)s <input>...  [-c CUTPE -e ESLOPE -f FLUX -l LABEL -m MINPE -o OUTPUT -s -t FILETYPE -w WIDTH] [-V...]

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
  -o --output OUTPUT         If defined, this name will be used to save the plot(s) and, eventually npz file
  -V --verbose               Level of verbose output

Plot options:
  -c --cutpe CUTPE           Set a minimum number of pe to cut the data for calculate the average/std.dev. [default: 5]   
  -e --eslope ESLOPE         Slope used to simulate the energy. If not set, the counts of the events of a given energy
                             are used to balance the weights. This option is ignored if FLUX is 'none', 'None', or ''
  -f --flux FLUX             What flux weigths the events. Set 'none', 'None', or '' to not weight [default: crab]
  -l --label LABEL           Label for the Effective Area plot
  -m --minpe MINPE           Minumum number of pe to prduce all the distributions [default: 0]
  -s --savenpz               Save npz file of the variables/fit parameters, IF <output> is defined
  -t --filetype FILETYPE     Type of file (npz and root are for now available; it is the file extension!)
                             [default: npz]
  -w --width-bins WIDTH      Width size of the binning (in log scale) [default: 0.1]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from root2npz import *

if has_mpl:
    pyp.ion()

def get_binning(data, wbin=0.1, minpe=0.):
    pemaxXwbin = int(np.ceil(np.log10(data.max())/wbin))
    bins, peInPixels = np.linspace(1., pemaxXwbin*wbin, pemaxXwbin+1-int(np.ceil(1/wbin))), []
    bins = np.append([0.,np.log10(5)],bins)
    return bins[(bins>=np.log(minpe)) if minpe>0. else (~np.isinf(bins))]

def plotPEInPixel(pe, bins_kwd, wflux=None, savename=None, fig_and_plot=None, savenpz=False, extrainfo=None, verbose=False, **hist_kwd):
    weights = None if wflux is None else np.array(wflux).repeat(pe.shape[1])
    Bins = 10**get_binning(pe, **bins_kwd)#wbin=float(args['--width-bins']), minpe=minpe
    f, p = pyp.subplots() if fig_and_plot is None or len(fig_and_plot[:2])<2 else fig_and_plot
    histtype, color, ls, lw, fc = hist_kwd.pop('histtype','step'), hist_kwd.pop('color','b'), hist_kwd.pop('ls','solid'), hist_kwd.pop('lw',2.), hist_kwd.pop('fc',None)
    hist, bins, _ = p.hist(pe.flatten(), bins=Bins, weights=weights, normed=True, histtype=histtype,color=color, fc=fc, ls=ls, lw=lw, **hist_kwd)
    xlims = [np.floor(Bins[0]), np.ceil(Bins[-1])]
    ylims = [np.min(hist[hist>0.])/2., 1.]
    #if fig_and_plot is not None and len(fig_and_plot[:2])==2:
    #    ylims[0] = min(ylims[0], p.get_ylim()[0])
    f.suptitle("PE in pixel%s"%('' if extrainfo is None else (' (%s)'%extrainfo)), size="xx-large", weight='bold')
    if fig_and_plot is None or len(fig_and_plot[:2])<2:
        set_axis(p, 'x', r"p.e.", lims=xlims, scale='log', scalar_values=True)
        set_axis(p, 'y', r"fraction of p.e.s per pixel", lims=ylims, scale='log')
    if hist_kwd.get('label',None) is not None:
        p.legend(bbox_to_anchor=(0,0),loc='lower left')
    if savename is not None:
        saveplot(f, savename, extraname='_peInPixel', verbose='peInPixel' if verbose else False, **(dict(bins=bins, hist=hist) if savenpz else {}))
    pyp.pause(1e-20)
    return f, p

def plotPEInCamera(pe, bins_kwd, wflux=None, savename=None, fig_and_plot=None, savenpz=False, extrainfo=None, verbose=False, **hist_kwd):
    pe =  np.sum(pe, axis=1)
    Bins = 10**get_binning(pe, **bins_kwd)#wbin=float(args['--width-bins']), minpe=minpe
    f, p = pyp.subplots() if fig_and_plot is None or len(fig_and_plot[:2])<2 else fig_and_plot
    histtype, color, ls, lw, fc = hist_kwd.pop('histtype','step'), hist_kwd.pop('color','b'), hist_kwd.pop('ls','solid'), hist_kwd.pop('lw',2.), hist_kwd.pop('fc',None)
    hist, bins, _ = p.hist(pe, bins=Bins, weights=wflux, normed=True, histtype=histtype,color=color, fc=fc, ls=ls, lw=lw, **hist_kwd)
    xlims = [np.floor(Bins[0]), np.ceil(Bins[-1])]
    ylims = [np.min(hist[hist>0.])/2., 1.]
    #if fig_and_plot is not None and len(fig_and_plot[:2])==2:
    #    ylims[0] = min(ylims[0], p.get_ylim()[0])
    f.suptitle("PE in camera%s"%('' if extrainfo is None else (' (%s)'%extrainfo)), size="xx-large", weight='bold')
    if fig_and_plot is None or len(fig_and_plot[:2])<2:
        set_axis(p, 'x', r"p.e.", lims=xlims, scale='log', scalar_values=True)
        set_axis(p, 'y', r"fraction of p.e.s in camera", lims=ylims, scale='log')
    if hist_kwd.get('label',None) is not None:
        p.legend(bbox_to_anchor=(0,0),loc='lower left')
    if savename is not None:
        saveplot(f, savename, extraname='_peInCamera', verbose='peInCamera' if verbose else False, **(dict(bins=bins, hist=hist) if savenpz else {}))
    pyp.pause(1e-20)
    return f, p

def plotEDistribution(energies, bins, weights=None, title=None, savename=None, compare=None, **hist_kwd):
    f, p = pyp.subplots()
    histtype, color, ls, lw, fc = hist_kwd.pop('histtype','step'), hist_kwd.pop('color','b'), hist_kwd.pop('ls','solid'), hist_kwd.pop('lw',2.), hist_kwd.pop('fc',None)
    hist, bins, _ = p.hist(energies, bins=bins, weights=weights, histtype=histtype, color=color, fc=fc, ls=ls, lw=lw, **hist_kwd)
    if compare is not None:
        x = np.linspace(bins[0],bins[-1],1000)
        p.plot(x, compare(x), 'r--', label='reference')
        #p.hist(compare((bins[1:]-bins[:-1])/2.), bins=bins, label='reference', histtype=histtype, color='r', fc=fc, ls='dashed', lw=lw, **hist_kwd)
        p.legend()
    f.suptitle(r'%s'%title, weight='bold', size='xx-large')
    p.set_xlabel('E [TeV]', fontsize='xx-large')
    p.xaxis.set_label_coords(.5, -0.06)
    p.set_xscale('log')
    p.set_ylabel('a.u.', fontsize='xx-large')
    p.set_yscale('log')
    if savename is not None:
        f.savefig(savename+'.eps')
        if verbose:
            print " E distribution '%s' saved as %s"%(title, savename+'.eps')
    pyp.pause(1e-20)

def plotPEstats(energies, pe, stat, minpe=1., cutpe=5., extra=None, summed=False, savename=None, savenpz=False):
    stat = stat.lower()
    if stat not in ['mean','std']:
        print("In plotPEstats, stats is %s, which is not available..."%stat)
        return
    elif stat=='mean':
        fstat = '%s()'%stat
        uncertainties = lambda n, data2, toerrs: np.sqrt(toerrs/n) # toerrs = var (unbiased!)
    else:
        fstat = '%s(ddof=1)'%stat
        uncertainties = lambda n, data2, toerrs: np.sqrt(((toerrs - data2**2*(n-3)*(n-1)/n**2)/n) / (4*data2)) # toerrs=4-th central moment; data2 = data**2 (unbiased!)
    indeces = np.argsort(energies)
    energies = energies[indeces]
    pe = pe[indeces]
    if summed:
        pe = pe.sum(axis=1)
    e, ie = np.unique(energies, return_index=True)
    emin, ecut, data, dcut, to_error, to_ercut, Nmin, Ncut, n = [], [], [], [], [], [], [], [], 1
    if stat[:3]=='mean':
        uncertainty
    for i in ie:
        pe_selected = pe[i:ie[n] if len(ie)>n else len(energies)]
        pe_datamin, pe_datacut = pe_selected[pe_selected>=minpe], pe_selected[pe_selected>=cutpe]
        n+=1
        Nmin += [len(pe_datamin)]
        Ncut += [len(pe_datacut)]
        data += [eval('pe_datamin.%s'%fstat)]
        dcut += [eval('pe_datacut.%s'%fstat)]
        to_error += [pe_datamin.var(ddof=1) if stat=='mean' else np.sum((pe_datamin-pe_datamin.mean())**4)/Nmin[-1]]
        to_ercut += [pe_datacut.var(ddof=1) if stat=='mean' else np.sum((pe_datacut-pe_datacut.mean())**4)/Ncut[-1]]
    bad_min = (np.isnan(data))|(np.isnan(to_error))|(Nmin<2)
    bad_cut = (np.isnan(dcut))|(np.isnan(to_ercut))|(Ncut<2)
    emin, Nmin, data, to_error = e[~bad_min], np.array(Nmin)[~bad_min], np.array(data)[~bad_min], np.array(to_error)[~bad_min]
    ecut, Ncut, dcut, to_ercut = e[~bad_cut], np.array(Ncut)[~bad_cut], np.array(dcut)[~bad_cut], np.array(to_ercut)[~bad_cut]
    uncertmin, uncertcut = uncertainties(Nmin, data**2, to_error), uncertainties(Ncut, dcut**2, to_ercut)
    kwd_min = dict(fmt='-', marker='o', color='b', mfc='b', mec='b', label='all data' if minpe==0 else ('pe$\geq$%s pe'%int(minpe)))
    kwd_cut = dict(fmt='--', marker='s', color='r', mfc='r', mec='r', label=r'pe$\geq$%d pe'%cutpe)
    extra = '' if extra is None else ' (%s)'%extra
    ylabel = r'\langle pe\rangle' if stat=='mean' else r'\sqrt{Var(pe)}'
    values_xformat = lambda x, pos: (r"$%s$" % ((r"10^{%i}" % np.log10(x)) if x<=1e-3 or x>1000 else int(x) if x>=1 else float(x))).replace('-',u'\u2011')
    values_yformat = lambda x, pos: (r"$%s$" % ((r"10^{%i}" % np.log10(x)) if x<1e-3 or x>=1000 else int(x) if x>=1 else float(x))).replace('-',u'\u2011')
    f, p = pyp.subplots()
    ax_min = p.errorbar(emin, data, uncertmin, **kwd_min)
    ax_cut = p.errorbar(ecut, dcut, uncertcut, **kwd_cut)
    p.set_xscale('log')
    p.set_yscale('log')
    p.set_xlabel('energy [TeV]')
    p.set_ylabel(r'$\mathbf{%s}$%s'%(ylabel,extra))
    p.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(values_xformat))
    p.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(values_yformat))
    p.xaxis.set_label_coords(.5, -0.06)
    p.yaxis.set_label_coords(-0.09,.5)
    p.legend(bbox_to_anchor=(1.,0.), loc='lower right')
    if savename is not None:
        npz_kwds = dict(energy=e, data=data, data_cut=dcut, uncertaintis=uncertmin, uncertaintis_cut=uncertcut) if savenpz else {}
        saveplot(f, savename, extraname='_%s'%stat, ext='pdf', verbose='pe %s%s'%(stat,extra) if verbose else False, **npz_kwds)
    pyp.pause(1e-20)

if __name__=='__main__':
    files = sorted([ifile for path in args['<input>'] for ifile in glob.glob(aspath(path))])
    savename = args['--output']
    if savename is not None:
        savename = aspath(savename)
        savedir = dirname(savename)
        if savedir != '':
            mkdirs(savedir)
    verbose = args['--verbose']
    if args['--eslope'] is None:
        def wflatten(e):
            _,bounds = np.unique(e,True)
            bounds = np.append(bounds,[len(e)])
            sizes = np.array(np.diff(bounds),dtype=np.float_)
            weight = []
            for s in sizes:
                weight += [1./s]*int(s)
            #weight = np.array(weight)/float(len(sizes))#*float(len(e))
            #print weight.sum()
            return np.array(weight)#/float(len(sizes))#*float(len(e))
    else:
        wflatten = lambda e: e**(-float(args['--eslope']))
    savenpz = args['--savenpz']
    extralabel = '' if args['--label'] is None else args['--label']+': '
    minpe = float(args['--minpe'])
    cutpe = float(args['--cutpe'])
    # calculate
    lge, tiggeredTel, iPhotons, pixelPE = build_photonInfo(files, None, verbose)
    indeces = np.argsort(lge)
    lge, tiggeredTel, iPhotons, pixelPE = lge[indeces], tiggeredTel[indeces], iPhotons[indeces], pixelPE[indeces]
    selection = np.any(pixelPE>=0., axis=2)
    pe =  pixelPE[selection]
    #bins = get_binning(pe, wbin=float(args['--width-bins']), minpe=minpe)
    bins_kwd = dict(wbin=float(args['--width-bins']), minpe=minpe)
    if all(lge==lge[0]):
        multi_energy, wflux = False, None
        label = '%s$E = %1.2f$ TeV'%(extralabel, 10**lge[0])
        print ">>> Average energy = %1.2f TeV <<<" % (10**lge[0])
    else:
        multi_energy, e = True, np.round(10**(lge[selection.flatten()]),2)
        wflux = None if args['--flux'].lower() in ['none',''] else (wflatten(e)*eval(args['--flux']+'_flux')(e))
        label = '%s$%1.2f \leq E \leq %1.2f$ TeV%s'%(extralabel, e.min(), e.max(), ('' if wflux is None else ('\n(distributed as %s flux)'%args['--flux'])))
        print ">>> Average energy = %1.2f TeV" % (np.sum(e*wflux)/np.sum(wflux))
        try:
            energies = np.array(sorted([float(basename(f).split('_')[1][:-3]) for f in files]))
        except:
            steps = 10**np.linspace(-1,2.6, 37)
        else:
            steps = (energies[:-1] + energies[1:])/2.
            steps = np.append(steps, 2*energies[-1]-steps[-1])
            steps = np.append(2*energies[0]-steps[0], steps)
        #plotEDistribution(e, steps, weights=wflatten(e), title='Normalized to a flatten energy distribution...', savename=savename+'_Eflatten')
        #plotEDistribution(e, steps, weights=wflux, title='...and weighted to the crab', savename=savename+'_Ecrab',
        #                  compare=None if args['--flux'].lower() in ['none',''] else eval(args['--flux']+'_flux'))
        plotPEstats(e, pe, 'mean', minpe=minpe, cutpe=cutpe, savename=savename+'pixels', extra='per pixel', savenpz=savenpz)
        plotPEstats(e, pe, 'std', minpe=minpe, cutpe=cutpe, savename=savename+'pixels', extra='per pixel', savenpz=savenpz)
        plotPEstats(e, pe, 'mean', minpe=minpe, cutpe=cutpe, savename=savename+'camera', extra='on camera', summed=True, savenpz=savenpz)
        plotPEstats(e, pe, 'std', minpe=minpe, cutpe=cutpe, savename=savename+'camera', extra='on camera', summed=True, savenpz=savenpz)
        raw_input('...Press any key to continue...')
        pyp.close('all')
    extrainfo = None if wflux is None else ("folded with %s flux"%args['--flux'])
    fPix, pPix = plotPEInPixel(pe, bins_kwd, wflux=wflux, savename=savename, savenpz=savenpz, extrainfo=extrainfo, verbose=verbose, label=label)
    fCam, pCam = plotPEInCamera(pe, bins_kwd, wflux=wflux, savename=savename, savenpz=savenpz, extrainfo=extrainfo, verbose=verbose, label=label)
    if multi_energy:
        i=0
        c, ls = ['r','g','m'], ['--',':','-.']
        for ENE in energies:
            ENE = float(ENE)
            if 10**i<=ENE and 10**(i+1)>ENE:
                lbl = '%s$E = %1.2f$ TeV'%(extralabel, ENE)
                pe_data = pe[e==ENE]
                plotPEInPixel(pe_data, bins_kwd, fig_and_plot=(fPix,pPix), savename=savename, extrainfo=extrainfo, verbose=verbose, label=lbl, color=c[i], ls=ls[i])
                plotPEInCamera(pe_data, bins_kwd, fig_and_plot=(fCam,pCam), savename=savename, extrainfo=extrainfo, verbose=verbose, label=lbl, color=c[i], ls=ls[i])
                i+=1
    if has_mpl:
        raw_input(" => Press enter to close... <=")
        pyp.close('all')
        pyp.ioff()


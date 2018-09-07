#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys
    #import sys, os, glob
    #from root_numpy import tree2rec

    version = "Version: 0.3 (A. Porcelli, ale_led@yahoo.it)"
    version = "Version: 0.4 (I. Al Samarai)"
    doc = """
Make effective area

Usage:
  %(prog)s <input> [<output> -c CONEVIEW -d DIST... -e EWIDTH -f -l LABEL -n -r RWIDTH -isu] [-C | -G | -P] [-V...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <input>                    Input file(s). Format compatible to TChain for multiple parsing.
  <output>                   If defined, this name will be used to save plots and, eventually npz file

General options:
  -s --savenpz               Save npz file of the variables/fit parameters, IF <output> is defined
  -V --verbose               Level of verbose output

Plot options:
  -d --distance-cut DIST     Max distance in meters to select the the data for the trigger efficiency.
                             If not set, no selection will occur. Multiple settings allows a comparison
                             in the same plot (only the first 7 will be considered) and the value will
                             be used for the label   
  -e --energy-width EWIDTH   Width of the log10(energy) bins [default: 0.25] 
  -f --fit                   Use fitting functions
  -i --integral              Calculate the integral value
  -l --label LABEL           Label for the Effective Area plot
  -n --nevents               Number of events are shown in the R-efficiency plots
  -r --radius-width RWIDTH   Width of the radius bins in meters [default: 100]
  -u --fit-uncertainties     Plot the propagated uncertainty of the fits

Trigger rate:
  -c --coneview CONEVIEW     Calculate the proton trigger rate with the as efficiency folded p-flux times
                             2pi*[1-cos(CONEVIEW/2.)]. CONEVIEW is the angulare aperture in deg [default: 9]
  -C --cosmic-rays           To estimate approximately the CR differential trigger rate (use the 1.5 * proton flux)
  -G --gamma                 To estimate the gamma differential trigger rate (use the crab flux)
  -P --proton                To estimate the proton differential trigger rate (use the proton flux)

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from root2npz import *

def calc_effectiveArea(lge, r, r_triggered, lge_width=0.25, r_width=100, r_selections=None, xlabel='energy [TeV]', only_triggerEff=False,
                       nevent=False, savename=None, fittings=False, verbose=False, savenpz=False, **errorbar_kwd):
    drpi2 = pi2*r_width
    E, Etriggs, Ecounts, Aeff, AeffError, thisLge = [], [], [], [], [], np.floor(lge.min()/lge_width)*lge_width
    nextLge = thisLge
    while len(lge[lge>=nextLge])>0:
        nextLge += lge_width
        e_select = (lge>=thisLge)&(lge<nextLge)
        ene_bin = lge[e_select]
        if len(ene_bin)>1:
            rtriggered_ebin = r_triggered[e_select]
            rtrig_ebin = rtriggered_ebin[rtriggered_ebin>-1]
            if len(rtrig_ebin)>0:
                r_ebin = r[e_select]
                E += [ene_bin.mean()]
                if r_selections is None:
                    etrigg = [len(rtrig_ebin)]
                    ecount = [len(ene_bin)]
                else:
                    etrigg, ecount = [], []
                    for rsel in r_selections:
                        rselection = r_ebin<=float(rsel)
                        etrigg += [len(rtriggered_ebin[(rselection)&(rtriggered_ebin>-1)])]
                        ecount += [len(r_ebin[rselection])]
                Etriggs += [etrigg]
                Ecounts += [ecount]                    
                area, areaError, thisR, nextR, n_mult, iR, lenR = 0., 0., 0., 0., 1, 0, len(r_ebin)
                Rs, Rerrors, triggs, triggErrors, multi, N = [],[],[],[],[],[]
                while iR < lenR:
                    nextR += r_width
                    r_rbin = r_ebin[(r_ebin>=thisR)&(r_ebin<nextR)]
                    if len(r_rbin)>1:
                        Rs += [r_rbin.mean()]
                        Rerrors += [r_rbin.std(ddof=1)]
                        triggProb, triggError = clopper_pearson_binomial(len(rtrig_ebin[(rtrig_ebin>=thisR)&(rtrig_ebin<nextR)]), len(r_rbin))
                        triggs += [triggProb]
                        triggErrors += [triggError[0]]
                        multi += [n_mult]
                        N += [len(r_rbin)]
                        n_mult = 1
                        thisR = nextR
                    else:
                        n_mult += 1
                    iR += len(r_rbin)
                Rs, Rerrors, triggs, triggErrors, multi = np.array(Rs), np.array(Rerrors), np.array(triggs), np.array(triggErrors), np.array(multi)
                triggsErrs = np.swapaxes(triggErrors,0,1)
                if has_mpl and not only_triggerEff:
                    if verbose>1:
                        print "   ...plot the R-efficiency (%s<=log10(E)<%s)" % (thisLge, nextLge)
                    f,p = Refficiency(Rs, triggs, triggsErrs, Rerrors, thisLge, nextLge, N if nevent else None, verbose=verbose>1,
                                      savename=None if fittings else savename, show=verbose and not fittings, savenpz=savenpz)
                    if fittings:
                        # *--- TO BE FIXED
                        chi2norm, rtrigg_f = np.inf, None#rtrigg_distribution[0]
                        for i,func in enumerate(rtrigg_distribution):
                            result, thiscov, thischi2, thisndof = ChiSquareFit(eval(func), Rs, triggs, triggsErrs, init_rtrigg[i], zeroes=1e-3, bounds=bounds_rtrigg[i])
                            if thischi2 is None:
                                if func==rtrigg_distribution[-1] and rtrigg_f is None:#np.isinf(chi2norm):
                                    raise RuntimeError(result)
                                continue
                            thischi2norm = thischi2/thisndof
                            betterchi2 = thischi2norm<chi2norm if chi2norm>=1. else thischi2norm>chi2norm
                            if betterchi2:# or rtrigg_f is None:
                                chi2norm, pars, cov, chi2, ndof, index, rtrigg_f = thischi2, result, thiscov, thischi2, thisndof, i, rtrigg_distribution[i]
                        if has_mpl:
                            fit_shower(rtrigg_f, pars, cov, Rs, chi2, ndof, f_and_p=(f,p), fit_errors=None, savenpz=savenpz,
                                       savename=savename+'_%s-%s'%(thisLge, nextLge), verbose=verbose>1)
                        # ---*
                    if verbose>2:
                        raw_input(' press enter to continue...')
                    elif verbose>1:
                        pyp.pause(verbose)
                    pyp.close(f)
                    del f,p
                # area of this bin (dArea*dr*trig = 2pi*r*dr*trig => Aeff = 2pidr*sum(r*trigg); or Aeff = 2pi*int(r*rfit*dr) = 2pi*I
                # area uncert of this bin (DAeff = D(sum(dArea*dr*trig)) = 2pidr*sum(Dr*trig + r*Dtrig); or DAeff = D(2pi*r*rfit) = 2pi*int(drfit*r + dr*rfit) = 2pi*Ierr
                if fittings:
                    I,Ierr = quad(lambda r: r*rtrigg_f(r,*pars), 0, 5000)
                    Aeff += [pi2*I]
                    AeffError += [[pi2*Ierr]*2]
                else:
                    Aeff += [drpi2*np.sum(Rs*triggs*multi)]
                    AeffError += [[np.sum((Rerrors*triggs+Rs*triggErrors[:,0])*multi)*drpi2,np.sum((Rerrors*triggs+Rs*triggErrors[:,1])*multi)*drpi2]]
                thisLge = nextLge
                if len(lge[lge>=nextLge])==0:
                    break
    E, Etriggs, Ecounts = 10**np.array(E), np.swapaxes(Etriggs,0,1), np.swapaxes(Ecounts,0,1)
    if verbose:
        print "...plot the trigger efficiency"
    if r_selections is not None:
        labels = [('' if errorbar_kwd.get('label',None) is None else (errorbar_kwd.get('label')+': '))+'r$\leq$%d m'%int(round(float(rsel))) for rsel in r_selections]
        errorbar_kwd['label'] = labels
    plot_TriggerEfficiency(E, Etriggs, Ecounts, savename=savename+('_pe' if only_triggerEff else ''), savenpz=savenpz, xaxis='pe' if only_triggerEff else 'E',
                           xlabel=xlabel, **errorbar_kwd)
    return E, np.array(Aeff), np.swapaxes(AeffError,0,1)

def plot_effArea(E, Aeff, AeffError, flux, solid_angle, savename=None, factor=1.5, integral=True, show_uncertainties=False, label=None, savenpz=False, verbose=False):
    AEpopt, AEpcov, AEchi2, AEndof = ChiSquareFit(eval(effA_fit), np.log10(E), np.log10(Aeff), np.log10(AeffError), init_effA)#, bounds=bounds_rtrigg[i])
    #savename = None if outname is None else os.path.join(savedir,outname+'_effectiveArea')
    if savenpz:
        savedata(savename+'_effectiveArea', verbose=verbose, E=E, Aeff=Aeff, AeffError=AeffError)
    if AEchi2 is None:
        print AEpopt
    if has_mpl:
        AEfig, AEplot = EffectiveArea(E, Aeff, AeffError, savename=savename if AEchi2 is None else None, show=AEchi2 is None, label=label, savenpz=savenpz, verbose=verbose)
    uncert_propagation = dict(chi2=AEchi2, ndof=AEndof, fit_errors=effA_errScaled) if AEchi2 is not None and show_uncertainties else dict()
    scaled_f = "lambda e, *pars: 10**eval(effA_fit)(np.log10(e), *pars)"
    if AEchi2 is not None and has_mpl:
        fit_shower(scaled_f, AEpopt, AEpcov, E, f_and_p=(AEfig, AEplot), savename=savename+'_effectiveArea', savenpz=savenpz,
                   verbose=verbose, bbox_to_anchor=([0.0,1.0]), loc="upper left", **uncert_propagation)
    if not (None in [solid_angle,flux]):
        analytic = dict(rate=scaled_f, rate_err=uncert_propagation.get('fit_errors',None), pars=AEpopt, cov=AEpcov, chi2=None, ndof=None)
        if show_uncertainties:
            analytic.update(dict())
    return differentialTriggerRate(E, Aeff, AeffError, flux, solid_angle, factor=factor, integral=integral, analytic=analytic if AEchi2 is not None else None,
                                   savename=savename, savenpz=savenpz, verbose=verbose, label=label)

if __name__=='__main__':
    files = aspath(args['<input>'])
    savename = args['<output>']
    if savename is not None:
        savename = aspath(savename)
        savedir = dirname(savename)
        if savedir != '':
            mkdirs(savedir)
    verbose = args['--verbose']
    lge_width, r_width = float(args['--energy-width']), float(args['--radius-width']) # m
    integral, fittings, savenpz = args['--integral'], args['--fit'], args['--savenpz']
    if args.get('--gamma',False):
        flux, factor, solid_angle = crab_flux, 1., 1.
    elif args.get('--proton',False):
        flux, factor, solid_angle = proton_flux, 1., pi2*(1-np.cos(deg2rad*float(args['--coneview'])/2.))
    elif args.get('--cosmic-rays',False):
        flux, factor, solid_angle = proton_flux, 1.5, pi2*(1-np.cos(deg2rad*float(args['--coneview'])/2.))
    else:
        flux, factor, integral, solid_angle = (None,None,False,1.)
    # Get data
    lge, r, r_triggered, pe, ph = build_effArea_variables(files, savename if savenpz else None, verbose)
    if pe is not None:
        pesel = pe>0.
    if ph is not None:
        phsel = ph>0.
    infoprint = "{:>35}: {:<}"
    print "\n {:=^19} Data info {:=^19}".format('','')
    print infoprint.format("Events", len(lge))
    if ph is not None:
        print infoprint.format("Events with photons in camera", len(ph[phsel]))
    if pe is not None:
        print infoprint.format("...p.e.s in camera", len(pe[pesel]))
    if ph is not None and pe is not None:    
        print infoprint.format("...photons lost (PDE)", len(ph[(pe==0)&(phsel)]))
        print infoprint.format("...p.e.s not from photons (=0!)", len(ph[(pesel)&(ph==0)]))
    print infoprint.format("Triggered events", len(lge[r_triggered!=-1]))
    if pe is not None:
        print infoprint.format("...by event with p.e. on camera", len(pe[(pesel)&(r_triggered!=-1)]))
        print infoprint.format("...by NSB, dark count, etc.", len(pe[(pe==0)&(r_triggered!=-1)]))
    print " {:=^49}\n".format('')
    E, Aeff, AeffError = calc_effectiveArea(lge, r, r_triggered, lge_width=lge_width, r_width=r_width,
                                            r_selections=None if len(args['--distance-cut'])==0 else args['--distance-cut'],
                                            nevent=args['--nevents'], savename=savename, fittings=fittings, verbose=verbose, savenpz=savenpz,
                                            label=args['--label'])
    if pe is not None:
        calc_effectiveArea(np.log10(pe[pesel]),r[pesel], r_triggered[pesel], lge_width=lge_width, r_width=r_width,
                           r_selections=None if len(args['--distance-cut'])==0 else args['--distance-cut'], xlabel='p.e.', only_triggerEff=True,
                           nevent=args['--nevents'], savename=savename, fittings=fittings, verbose=verbose, savenpz=savenpz, label=args['--label'])
    if verbose:
        print "...plot the effective area"
    plot_effArea(E, Aeff, AeffError, flux, solid_angle, savename, factor=factor, integral=integral, verbose=verbose>1,
                 show_uncertainties=args['--fit-uncertainties'], label=args['--label'], savenpz=savenpz)
    if has_mpl:
        raw_input(" => Press enter to close... <=")
    
    
    

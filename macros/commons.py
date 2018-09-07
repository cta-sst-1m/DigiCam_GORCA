#!/usr/bin/env python

import sys, os, glob, time
import numpy as np
from scipy import interpolate
from scipy.special import erf, erfc, betaincinv
from scipy.stats import norm
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.optimize import curve_fit, minimize
import ROOT as R
try:
    import matplotlib as mpl
    has_mpl = True
    import matplotlib.pyplot as pyp 
    rc, rcParams = mpl.rc, mpl.rcParams
    rcParams['lines.linewidth'] = 2.0
    rcParams['lines.markeredgewidth'] = 1.0
    rcParams['font.family'] = 'Times New Roman'
    rcParams['axes.titlesize'] = 'xx-large'
    rcParams['axes.labelweight'] = 'bold'
    rcParams['axes.formatter.useoffset'] = False
    rcParams['axes.prop_cycle'] = mpl.cycler('color', 'brkgmcy')
    rcParams['xtick.major.size'] = 8
    rcParams['xtick.minor.size'] = 4
    rcParams['xtick.major.width'] = 1
    rcParams['xtick.labelsize'] = 'xx-large'
    rcParams['ytick.major.size'] = 8
    rcParams['ytick.minor.size'] = 4
    rcParams['ytick.major.width'] = 1
    rcParams['ytick.labelsize'] = 'xx-large'
    rcParams['legend.fancybox'] = True
    rcParams['legend.numpoints'] = 1
    rcParams['legend.fontsize'] = 'x-large'
    rcParams['legend.handletextpad'] = 0.5
    rcParams['legend.shadow'] = True
    rcParams['legend.frameon'] = False
    rcParams['legend.scatterpoints'] = 1
    rcParams['figure.titlesize'] = 'xx-large'
    rcParams['figure.titleweight'] = 'bold'
    rcParams['figure.facecolor'] = 'w'
    rcParams['errorbar.capsize'] = 4
    rcParams['savefig.format'] = 'eps'

except ImportError:
    has_mpl = False
from ProgressBar import *

pi2 = 2*np.pi
ln10 = np.log(10.)
deg2rad = np.pi/180.
cardinal = lambda n: ('%s'+cardinalext(int(n)%10 if int(n)//10!=1 else 0) if str(n).isdigit() else '%s-th')%n
cardinalext = lambda n: 'st' if n==1 else ('nd' if n==2 else ('rd' if n==3 else 'th'))
colors = ['b', 'r', 'k', 'g', 'm', 'c', 'y']
markers = ['o', 's', '^', 'v', '>', '<', '+']
lines = ['-', '--', ':', '-.', '-', '--',':']

# == is...
def isbool(boolean):
    return isinstance(boolean, bool)

def isfloat(value):
    return isinstance(value, float)

def isdict(dictionary):
    return isinstance(dictionary, dict)

def isint(value):
    return isinstance(value, int)

def isiterable(iterable):
    """True if <iterable> is a tuple, list or ndarray"""
    return isinstance(iterable, (tuple, list, np.ndarray))

def isreal(value):
    """True if <value> is an int or a float"""
    return isinstance(value, (int,float))

def isstr(string):
    return isinstance(string, str)

# shells
def mkdirs(*dirs):
    """Multiple <dirs> can be crreated if they don't exist"""
    for d in dirs:
        if not d=='' and not os.path.exists(d):
            os.makedirs(d)

def abspath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def aspath(*paths):
    """<paths> will be joined"""
    return os.path.normpath(abspath(os.path.join(*paths)))

def basename(path):
    return os.path.basename(path)

def dirname(path):
    return os.path.dirname(path)

# Functions
def clopper_pearson_binomial(passed, total, sigma=1, CL=None):
    """
    Estimate the exact binomial from the clopper-pearson method
    - paramters -
    passed: int or array
        counts of passed elements
    total: int or array
        total elements
    sigma: float [default: 1]
        to estimate the CL automatically from the normal distribution at <sigma> sigmas
    CL: None or float [default: None]
        to specify a confidence level for the clopper-pearson. If None, it will be automatically estimated by <sigma>

    - return -
    eff: float or array:
        efficiency (<passed>/<total>). If <passed> and <total> are arrays, <eff> is an array
    uncertainties: 1d or 2d array
        array of the uncentainties: <uncertainties>[0] is the lower boundary, <uncertainties>[1] the upper one.
        If <passed> and <total> are arrays, <uncertainties>[0] and <uncertainties>[1] are arrays
    """
    eff, notpassed = np.float_(passed)/total, total-passed
    if CL is None:
        ybeta_low, ybeta_up = 1-norm.cdf(sigma,0,1), norm.cdf(sigma,0,1)
    else:
        ybeta_low, ybeta_up = (1-CL)/2, (1+CL)/2
    el, eu = eff-betaincinv(passed,notpassed+1,ybeta_low),betaincinv(passed+1,notpassed,ybeta_up)-eff
    if isiterable(el):
        el[np.isnan(el)] = eu[np.isnan(eu)] = 0.
    else:
        if np.isnan(el):
            el = 0.
        if np.isnan(eu):
            eu = 0.
    return eff, np.atleast_2d([el,eu])

def has_hexagonal_cluster(triggered_list, neighbours_map):
    for triggered in triggered_list:
        neighbours = neighbours_map[triggered] # RIFARE LISTA GRUPPI PER ORDINARLI SEQUENZIALMENTE!!!
        selected_neighbs = [n for n in neighbours if n>triggered]
        # Need the 3 major neighbours triggered
        if len(selected_neighbs)>3 and all(neighb in triggered_list for neighb in selected_neighbs) and\
            all(neighb in triggered_list for neighb in neighbours_map[selected_neighbs[-1]]): # and the last one (the center of the cluster) must have all of its neigbours triggered
            return True
    return False
def events_with_hexagonal_cluster(triggered_events, neighbours_map):
    return np.array([has_hexagonal_cluster(event, neighbours_map) for event in triggered_events])

# Fitting! <f> must be as 'f(x,*pars)'
rtrigg_distribution = [#lambda v, middle, mid2ground, maximum: maximum*0.5*erfc((v-middle)*2./mid2ground),
                       "lambda v, *pars: pars[0]*(1.0+(v/pars[1])**(-pars[2]/pars[3]))**(-pars[3])"#,
                       #lambda v, *pars: pars[0]*(v**pars[1])*(1.0+(v/pars[2])**((pars[1]-pars[3])/pars[4]))**(-pars[4])
                       ]
#f_fit? = lambda x,*pars: np.where(func(x,*pars)>pars[0],pars[0]*np.ones_like(len(func(x,*pars))),np.where(func(x,*pars)<0.,np.zeros(len(func(x,*pars))),func(x,*pars)))

init_rtrigg = [#[500.0,300.0,.5],
               [1.0,350.0,-10.0,1.0]#,
               #[1.0,0.0,350.0,-10.0,1.0]
               ]
bounds_rtrigg = [#[(0.,1.),(None,None),(None,None)],
                 [(0.,1.),(1e-6,None),(None,-1e-6),(1e-6,None)]#,
                 #[(1e-8,None),(None,None),(1e-6,None),(None,-1e-6),(1e-6,None)]
                 ]
#bounds_rtrigg = [[(1e-8,1.),(None,None),(1e-6,None)],
#                 [(0.,1.),(-0.1,10.),(10,2000),(-50.,-1),(0.01,50.)],
#                 [(0.,1.),(10,2000),(-50.,-1),(0.01,50.)]]
effA_fit = "lambda x, logA,B,C,D: logA + B*x - np.log10(1+((10**x)/C)**D)"
init_effA = [5., 1., 2., .5]
effA_err = ["lambda x, logA,B,C,D: 1",
            "lambda x, logA,B,C,D: B",
            "lambda x, logA,B,C,D: D*10**(x*D)/(ln10*(C**(D+1))*(1+((10**x)/C)**D))",
            "lambda x, logA,B,C,D: -(((10**x)/C)**D)*np.log((10**x)/C)/(ln10*(1+((10**x)/C)**D))"]
effA_errScaled = ["lambda x, *pars: 10**eval(effA_err[0])(np.log10(x), *pars)",
                  "lambda x, *pars: 10**eval(effA_err[1])(np.log10(x), *pars)",
                  "lambda x, *pars: 10**eval(effA_err[2])(np.log10(x), *pars)",
                  "lambda x, *pars: 10**eval(effA_err[3])(np.log10(x), *pars)"]

def curve_fit_wrap(f, xs, ys, yerrs, init_pars, **curve_fit_kwd):
    xs, ys, yerrs = np.atleast_1d(xs), np.atleast_1d(ys), np.atleast_2d(yerrs)
    if yerrs[0].shape != ys.shape or yerrs[-1].shape != ys.shape:
        raise ValueError("shapes of ys and yerrs differ")
    elif len(yerrs)>1:
        def red_f(xdata, *pars):
            res_y = f(xdata, *pars)-ys
            return np.where(res_y<0., res_y/yerrs[0], res_y/yerrs[-1])
    else:
        def red_f(xdata, *pars):
            return (f(xdata, pars)-ys)/yerrs[0]
    popt, pcov = curve_fit(red_f, xs, np.zeros(len(ys)), p0=init_pars, sigma=np.ones_like(ys), absolute_sigma=True)
    chi2 = np.sum(red_f(xs, *popt)**2)
    return popt, pcov, chi2, len(ys)-len(popt)

def ChiSquareFit(f, xs, ys, yerrs, init_pars, zeroes=None, **minimize_kwd):
    xs, ys, yerrs = np.atleast_1d(xs), np.atleast_1d(ys), np.atleast_2d(yerrs)
    if yerrs[0].shape != ys.shape or yerrs[-1].shape != ys.shape:
        raise ValueError("shapes of ys and yerrs differ")
    elif len(yerrs)>1:
        def ChiSq(pars):
            res_y = f(xs, *pars)-ys
            lowers = res_y<0.#; print res_y
            uncerts = np.where(lowers, yerrs[0], yerrs[-1])
            if zeroes is None:
                non_zeroes = uncerts>0.
                res_y, uncerts = res_y[non_zeroes], uncerts[non_zeroes]
            else:
                uncerts[uncerts==0.] = zeroes
                res_y[uncerts==0.] = 1./zeroes
            nans = np.isnan(uncerts)
            res_y, uncerts = res_y[~nans], uncerts[~nans]
            return np.sum((res_y/uncerts)**2)
    else:
        def ChiSq(pars):
            non_zeroes = yerrs[0]>0.
            return np.sum(((f(xs, *pars)-ys)[non_zeroes]/yerrs[0][non_zeroes])**2)
    methods = [minimize_kwd['method']] if minimize_kwd.has_key('method') else [None,'Nelder-Mead']
    for m in methods:
        minimize_kwd['method'] = m
        optres = minimize(ChiSq, init_pars, options=dict(maxiter=1e20), **minimize_kwd)
        if optres.success:
            if m is not None and len(methods)>1:
                print("[minimizer] standard method failed, but '%s' is succesfull, instead!"%m)
            break
    if optres.success:
        popt = optres.x
        if minimize_kwd.get('method') in ['Nelder-Mead']:
            pcov = np.nan
        else:
            pcov = (optres.hess_inv if minimize_kwd.get('bounds',None) is None else optres.hess_inv.todense())*2.
        ndof = len(ys[np.where(f(xs, *popt)-ys<0., yerrs[0], yerrs[-1])>0.])-len(popt)
        return popt, pcov, optres.fun, ndof#len(ys)-len(popt)#
    else:
        return optres.message, None,None,None
    
def propagate_errorfunctions(x, pars, cov, funcs):
    ys = [f(x,*pars) for f in np.atleast_1d(funcs)]
    error2, n, cov = 0., len(ys), np.atleast_2d(np.float_(cov))
    for i,y in enumerate(ys):
        error2 += y**2*cov[i][i]
        for j in xrange(i+1,n):
            error2 += 2*y*ys[j]*cov[i][j]
    return np.sqrt(error2)

# Plots
def set_axis(plot, axis, title, lims=None, scale=None, scalar_values=False, **title_kwds):
    size = title_kwds.pop('size',None)
    size = title_kwds.pop('size','xx-large')
    weight = title_kwds.pop('weight','bold')
    if axis.lower() not in ['x','y']:
        raise ValueError("'axis' must 'x', or 'y'!")
    if size=='xx-large':
        if axis=='x':
            plot.xaxis.set_label_coords(.5, -0.06)
        else:
            plot.yaxis.set_label_coords(-0.105,.5)
    set_title = eval('plot.set_%slabel'%axis)
    set_title(title, fontsize=size, weight=weight, **title_kwds)
    if lims is not None:
        set_limits = eval('plot.set_%slim'%axis)
        if isinstance(lims,tuple):
            set_limits(*lims)
        elif isinstance(lims,dict):
            set_limits(**lims)
        else:
            set_limits(lims)
    if scale is not None:
        set_scale = eval('plot.set_%sscale'%axis)
        set_scale(scale)
        if scalar_values:
            ax = eval('plot.%saxis'%axis)
            if isstr(scalar_values):
                ax.set_major_formatter(mpl.ticker.FormatStrFormatter(scalar_values))
            else:
                ax.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: (r"$%s$" % (int(x) if x>1 else (float(x) if x>1e-3 else (r"10^{%i}" % np.log10(x))))).replace('-',u'\u2011')))

def saveplot(fig, savename, extraname='', ext='eps', set_title=False, extratitle='', verbose=True, **npz_kwds):
    if set_title:
        if isstr(set_title):
            fig.suptitle(set_title, size="xx-large", weight='bold')
        else:
            name = basename(savename).replace("_"," ") + extratitle
            fig.suptitle(r"%s"%name, size="xx-large", weight='bold')
    savename += extraname
    if len(npz_kwds)>0:
        savedata(savename, verbose=verbose, **npz_kwds)
    fig.savefig(savename+'.'+ext)
    if verbose:
        print("   ...%splot saved as %s"%('' if isbool(verbose) else "'"+verbose+"' ", savename+'.'+ext))

def savedata(savename, verbose=False, **npz_kwds):
    np.savez(savename+".npz", **npz_kwds)
    if verbose:
        print("   ...%sdata saved as %s"%('' if isbool(verbose) else "'"+verbose+"' ", savename+'.npz'))
    

def Refficiency(Rs, triggs, triggerrors, Rerrors, thisLge, nextLge, N=None, savename=None, show=True, savenpz=False, verbose=False):
    f,p = pyp.subplots()
    title = "%1.2f$\mathbf{\leq E/TeV<}$%1.2f"%(10**thisLge,10**nextLge)
    p.errorbar(Rs, triggs, triggerrors, xerr=Rerrors, mec='b',fmt='o',mfc='b')
    if N is not None:
        for i,n in enumerate(N):
            p.annotate(str(n), xy=(Rs[i], triggs[i]), xycoords='data', xytext=(Rs[i]+10, triggs[i]+.01), textcoords='data', size='small', va='bottom', ha='left')
    set_axis(p, 'x', r"distance from the core (r) [m]")
    set_axis(p, 'y', r"trigger efficiency")
    if savename is not None:
        saveplot(f, savename, extraname='_Refficiency_%s-%s'%(thisLge,nextLge), set_title=True, extratitle=': '+title, verbose=verbose)
    if show:
        f.show()
    return f, p

def EffectiveArea(E, Aeff, AeffError, savename=None, fig_and_plot=None, show=True, savenpz=False, verbose=False, **errorbar_kwd):
    f, p = pyp.subplots() if fig_and_plot is None else fig_and_plot
    marker, color = errorbar_kwd.pop('marker','o'), errorbar_kwd.pop('color','b')
    p.errorbar(E, Aeff, AeffError, marker=marker,fmt=marker,mec=color,mfc=color, **errorbar_kwd)#,label=label)
    set_axis(p, 'x', r"energy [TeV]", scale='log', scalar_values=True)
    set_axis(p, 'y', r"effective area [m$^2$]", scale='log')
    if errorbar_kwd.get('label',None) is not None and show:
        p.legend()
    if savename is not None:
        saveplot(f, savename, extraname='_effectiveArea', set_title=True, verbose=verbose, **(dict(E=E, Aeff=Aeff, AeffError=AeffError) if savenpz else {}))
    if show:
        f.show()
    return f, p

def differentialTriggerRate(E, Aeff, AeffError, flux, solid_angle, factor=1., integral=False, analytic=None, verbose=False, **EffectiveArea_kwd):
    rate = factor*flux(E)*solid_angle*Aeff
    rate_err = factor*flux(E)*solid_angle*AeffError
    imax = np.where(rate==np.max(rate))[0][0]
    ishort = E>=E[imax]
    Eshort = E[ishort]
    savename = EffectiveArea_kwd.pop('savename',None)
    savenpz = EffectiveArea_kwd.pop('savenpz',None)
    results = None
    if len(Eshort)<2 and analytic is None and integral:
        print("The energy array is less than 2 elements: no integral will be performed!")
        integral = False
        results = None, None
    fig, plot = EffectiveArea(E, rate, rate_err, savename=None, show=False, **EffectiveArea_kwd)
    plot.set_ylabel(r"differential trigger rate (dR/dE) [Hz TeV$^{-1}$]",fontsize='xx-large')
    has_rate_f = False
    if analytic is None:
        show = True
        e=E
    else:
        Aeff_f, Aeff_ferr, popt, pcov = analytic.get('rate',None), analytic.get('rate_err',None), analytic.get('pars',None), analytic.get('cov',None)
        chi2, ndof, show = analytic.get('chi2',None), analytic.get('ndof',None), False
        if None is Aeff_ferr or None is pcov:
            rate_ferr = pcov = None
        else:
            rate_ferr = [lambda E, *pars: factor*flux(E)*solid_angle*eval(Aeff_ferr[0])(E, *pars),
                         lambda E, *pars: factor*flux(E)*solid_angle*eval(Aeff_ferr[1])(E, *pars),
                         lambda E, *pars: factor*flux(E)*solid_angle*eval(Aeff_ferr[2])(E, *pars),
                         lambda E, *pars: factor*flux(E)*solid_angle*eval(Aeff_ferr[3])(E, *pars)]
        if None is Aeff_f or None is popt:
            print(" ...wrong analytic dict formatting: it will not be plotted!")
        else:
            rate_f = lambda E, *pars: factor*flux(E)*solid_angle*eval(Aeff_f)(E, *pars)
            has_rate_f = True
            e = np.linspace(np.min(E)*.9,np.max(E)*1.1,100000)
            ph = rate_f(e,*popt)
            plot.set_ylim(np.min(ph), np.max(ph)*10)
            imax = np.where(ph==np.max(ph))[0][0]
            fit_shower(rate_f, popt, pcov, E, chi2, ndof, f_and_p=(fig, plot), fit_errors=rate_ferr, show=False, verbose=verbose, bbox_to_anchor=([1.,1.]), loc="upper right")
            if imax>0 and imax<len(e)-1:
                print "\n >>> THE Eth IS  %1.3f +/- %1.3f TeV <<<" % (e[imax], (e[1]-e[0])/2.)
                plot.vlines(e[imax],np.min(ph), np.max(ph)*10, colors='g',linestyle='-.')
                plot.annotate(r"E$_{th}$ ="+"\n%1.3f TeV"%e[imax],xy=(e[imax],np.sqrt(np.min(ph)*ph[imax])), xytext=(e[imax]*10,np.sqrt(np.min(ph)*ph[imax])*10**(-1)),
                              size='x-large', ha="left", arrowprops=dict(arrowstyle="fancy",fc="g", ec="none",connectionstyle="angle3,angleA=90,angleB=0"))
    if integral:
        if has_rate_f:
            trigger_rate, trigger_rate_err = quad(lambda e: rate_f(e,*popt), e[0], e[-1])
            trigger_rate_short, trigger_rate_err_short = quad(lambda e: rate_f(e,*popt), e[imax], e[-1])
        else:
            #e = E
            dE = np.diff([(3*E[0]-E[1])/2.]+list((E[:-1]+E[1:])/2.)+[(3*E[-1]-E[-2])/2.])
            trigger_rate, trigger_rate_err = np.sum(rate*dE), np.sum(rate_err*dE)
            dEshort = np.diff([(3*Eshort[0]-Eshort[1])/2.]+list((Eshort[:-1]+Eshort[1:])/2.)+[(3*Eshort[-1]-Eshort[-2])/2.])
            trigger_rate_short, trigger_rate_err_short = np.sum(rate[ishort]*dEshort), np.sum((rate_err[:,ishort] if len(rate_err.shape)>1 else rate_err[ishort])*dEshort)
        print "\n >>> THE TRIGGER RATE IS  %1.3f +/- %1.3f Hz <<<" % (trigger_rate, trigger_rate_err)
        print "\n >>> THE TRIGGER RATE WITH E>Eth IS  %1.3f +/- %1.3f Hz <<<" % (trigger_rate_short, trigger_rate_err_short)
        results = (trigger_rate, trigger_rate_err) if imax==0 or imax==len(e)-1 else (trigger_rate_short, trigger_rate_err_short)
        if results[1]/results[0]>0.8:
            results = np.nan, np.nan
            print "Integral larger than it's erro! This result is set as NaN!"
    if savename is not None:
        saveplot(fig, savename, extraname='_differentialTriggerRate', set_title=True, verbose=verbose,
                 **(dict(E=E, rate=rate, rate_err=rate_err, solid_angle=solid_angle, factor=factor, Eth=e[imax]) if savenpz else {}))
    fig.show()
    return results

def plot_TriggerEfficiency(X, triggered, counts, savename=None, savenpz=False, verbose=False, xaxis='E', xlabel="energy [TeV]", **errorbar_kwd):
    labels = errorbar_kwd.pop('label',None)
    if len(triggered.shape)==1:
        counter = xrange(1)
        triggered, counts = [triggered], [counts]
    elif len(triggered.shape)==2:
        counter = xrange(min(len(triggered),len(colors)))
    else:
        raise ValueError("shape of 'triggered' must be at most a length 2 tuple (it is %s)"%triggered.shape)
    labels = labels if isiterable(labels) else [labels]*len(triggered)
    if triggered.shape!=counts.shape or len(labels)!=len(triggered):
        raise ValueError("shape of 'triggered' (%s), 'counts' (%s) and labels (%s) differs!"%(triggered.shape, counts.shape, len(label)))
    f, p = pyp.subplots()
    for i in counter:
        triggProb, triggError = clopper_pearson_binomial(triggered[i], counts[i])
        if has_mpl:
            marker, color = errorbar_kwd.pop('marker',markers[i]), errorbar_kwd.pop('color',colors[i])
            mec, mfc, fmt = errorbar_kwd.pop('mec',color), errorbar_kwd.pop('mfc',color),  errorbar_kwd.pop('fmt',lines[i])
            kwds = errorbar_kwd.copy()
            if labels[i] is not None:
                kwds['label'] = r"%s"%labels[i]
            p.errorbar(X, triggProb, triggError, marker=marker, fmt=fmt, mec=color, mfc=color, color=color, **kwds)
    if has_mpl:
        set_axis(p, 'x', r"%s"%xlabel, scale='log', scalar_values=True)
        set_axis(p, 'y', r"trigger efficiency")
        p.legend(bbox_to_anchor=(1,0), loc='lower right')
        if savename is not None:
            npzdict = dict(triggProb=triggProb, triggError=triggError)
            npzdict[xaxis] = X
            saveplot(f, savename, extraname='_triggerEff', set_title=True, verbose=verbose, **(npzdict if savenpz else {}))
        f.show()

def fit_shower(f, pars, cov, xs, chi2=None, ndof=None, par_names=None, f_and_p=None, fit_errors=None, savename=None, show=True, savenpz=False, verbose=False, **legend_kwd):
    show_summary, propagate_errors = (False,)*2
    fstring, f = (f, eval(f)) if isstr(f) else (None, f)
    if None not in [chi2,ndof]:
        chi2norm = chi2/ndof
        show_summary = True
        if cov is not None:
            cov *= chi2norm
            propagate_errors = True
    if show_summary:
        if par_names is None:
            par_names = ['par %s'%i for i in xrange(len(pars))]
        elif isiterable(par_names) and len(par_names)<len(pars):
            par_names += ['par %s'%i for i in xrange(len(par_names),len(pars))]
        elif not isiterable(par_names):
            par_names = [par_names]+['par %s'%i for i in xrange(1,len(pars))]
        print "Fit results: chi2/ndof = %3.2f/%d = %2.2f (used to weight the cov. matrix)"%(chi2,ndof,chi2norm)
        for i,p in enumerate(pars):
            print " -- '{name}' = {p:5.5} +/- {perr:5.5}".format(name=par_names[i], p=p, perr=np.sqrt(cov[i][i]))
    if f_and_p is not None:
        if isiterable(f_and_p):
            fig, ax = f_and_p if len(f_and_p)>1 else (pyp,f_and_p)
        else:
            fig, ax = pyp, pyp
        x_max, x_min = np.max(xs), np.min(xs)
        x_range = x_max-x_min
        x_sample = np.linspace(x_min*.9, x_max*1.1, 100001)
        y_sample = f(x_sample, *pars)
        label = 'analytical model' if None in [chi2, ndof] else ('$\chi^2/$ndof = %3.2f/%d = %2.2f'%(chi2,ndof,chi2norm))
        ax.plot(x_sample, y_sample, 'r-', label=label)
        if fit_errors is not None and cov is not None:
            fit_errors = np.atleast_1d(fit_errors)
            fit_errors = [eval(ferr) if isstr(ferr) else ferr for ferr in fit_errors]
            err_sample = propagate_errorfunctions(x_sample, pars, cov, fit_errors)
            ax.plot(x_sample, y_sample-err_sample, 'r--', label="uncert. %s"%('$\cdot$ $\chi^2/$ndof' if label != 'analytical model' else label), lw=1)
            ax.plot(x_sample, y_sample+err_sample, 'r--', lw=1)
        ax.legend(**legend_kwd)
        if savename is not None:
            saveplot(fig, savename, extraname='_fit', verbose=verbose,
                     **(dict(f=np.atleast_1d(fstring), pars=np.atleast_1d(pars), cov=np.atleast_1d(cov),
                             chi2=np.atleast_1d(chi2), ndof=np.atleast_1d(ndof)) if savenpz else {}))
        if show:
            fig.show()
        return f_and_p

# read bias curve log files (WARNING! Same thresholds scan must be common to these files!!!). NO .gz!
def getNSBRatesFromLogs(*logfiles):
    thresholds, rates, errors = [], [], []
    for log in logfiles:
        with open(log) as f:
            th,rt,er = [], [], []
            for l in f.readlines():
                if l.find("NSB Telescope Trigger rate at")>-1:
                    words = l.split()
                    th += [float(words[5])]
                    rter = words[-2].split('+-')
                    rt += [float(rter[0])]
                    er += [float(rter[1])]
        thresholds += [th]
        rates += [rt]
        errors += [er]
    return np.mean(thresholds,axis=0), np.mean(rates,axis=0), np.sqrt(np.sum(np.power(errors,2),axis=0)/float(len(errors)))

# Physics
proton_flux = lambda E: 1.49e7*(E*1e3+2.15*np.exp(-0.21*np.sqrt(E*1e3)))**(-2.74) # m**(-2) * s**(-1) * TeV**(-1) * sr**(-1); E is in TeV
crab_flux = lambda E: (3.76e-7) * E**(-2.39) * np.exp(-E/14.3) # ph * m**(-2) * s**(-1) * TeV**(-1); E is in TeV; uncert.: (3.76+-0.07)e-7, 2.39+-0.03, 14.3+-2.1

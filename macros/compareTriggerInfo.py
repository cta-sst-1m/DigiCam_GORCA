#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys

    version = "Version: 0.1 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Plot NSB trigger rate from the log files in a folder

Usage:
  %(prog)s <input>...  [ -f FLUX -n NAME -o OUTPUT ] [ -l LABEL... ] [ -r RSEL... ] [ -V... ]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <input>                    Input folder containing the npz files of the TriggerInfo.py results..
                             If multiple <input>, each of them is considered a different rate to be plotted. Max supported inputs:
                             14 (are you sure you want more than 14 inputs? The plot it will too crowded and hard to be read...), beyond
                             that, color and markers will repeat and the result might be confusing...

General options:
  -o --output OUTPUT         If defined, this is the folder where all the plots will be saved
  -V --verbose               Level of verbose output

Setting options:
  -l --label LABEL           String used as label. It is 1 to 1 with the <input> list. If not given 'input <N>' will be used, whene <N>
                             is the number of the positional <input> arguments to which the legend is referring (starting from 1).
                             Note: if only an <input> is given, no legend is plotted.
  -f --flux FLUX             Flux to be user to convert the fit of the effective area to a differential trigger rate [default: crab]
  -n --name NAME             Root name of the files (porper appendix will be automatically appended) [default: gammas]
  -r --rsel RSEL             Cut on radius ti be inspected. If not given, no trigger efficiencies are produced

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from commons import *

if __name__=='__main__':

    line  = lambda i: ['-','--',':','-.'][i%4]
    mark = lambda i: ['o','s','^','v','<','>','x'][i%7]
    color = lambda i: ['b','r','k','g','y','c','m'][i%7]
    fill = lambda i: color(i) if i//7==0 else 'w'

    effAappendix, diffTrigEffAppendix = '_effectiveArea%s.npz', '_differentialTriggerRate.npz'
    triggEffAppendix_E, triggEffAppendix_pe = '_triggerEff.npz', '_pe_triggerEff.npz'

    folders = []
    for path in args['<input>'] :
        folders += glob.glob(aspath(path))
    folders = sorted(folders)
    flux = eval(args['--flux']+'_flux')
    labels = args['--label']
    name = args['--name']
    rsel = args['--rsel']
    outname = args['--output']
    if outname is not None:
        outname = aspath(outname)
        mkdirs(outname)

    feffA, peffA = pyp.subplots()
    fdifT, pdifT = pyp.subplots()
    ftr_E, ptr_E = [], []
    ftr_p, ptr_p = [], []
    for r in rsel:
        f, p = pyp.subplots()
        ftr_E += [f]
        ptr_E += [p]
        f, p = pyp.subplots()
        ftr_p += [f]
        ptr_p += [p]
    
    for i,path in enumerate(folders):
        effAfile = os.path.join(path,name+effAappendix)
        try:
            effA = np.load(effAfile%'')
        except IOError:
            print "file '%s' not found... let's try the next!"%(effAfile%'')
            continue
        try:
            difT = np.load(os.path.join(path,name+diffTrigEffAppendix))
        except IOError:
            print "file '%s' not found... let's try the next!"%os.path.join(path,name+diffTrigEffAppendix)
            continue
        peffA.errorbar(effA['E'], effA['Aeff'], effA['AeffError'], fmt='o', mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                       label=r"%s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
        pdifT.errorbar(difT['E'], difT['rate'], difT['rate_err'], fmt='o', mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                       label=r"%s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
        fitfile = effAfile%'_fit'
        if os.path.exists(fitfile):
            fit = np.load(fitfile)
            e = 10**np.linspace(np.log10(effA['E'][0]), np.log10(effA['E'][-1]), 10000)
            effAfit = eval(fit['f'][0])(e,*fit['pars'])
            peffA.errorbar(e, effAfit, color=color(i), fmt=line(i), marker='', label=r"%s fit"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
            pdifT.errorbar(e, effAfit*flux(e)*difT['factor']*difT['solid_angle'], color=color(i), fmt=line(i), marker='',
                       label=r"%s fit"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
        pdifT.axvline(difT['Eth'], color=color(i), ls=line(i), lw=1)
        vlim = pdifT.get_ylim()
        #ymin, ymax = vlim[0] if ylim['ymin'] is None else np.min(ylim['ymin'], vlim[0]), vlim[1] if ylim['ymax'] is None else np.max(ylim['ymax'], vlim[1])
        ytext = np.min(difT['rate'])*10**(i*2)#ymin*10**(i*1.5)
        pdifT.annotate(r"E$_{th}$ ="+"\n%1.3f TeV"%difT['Eth'], xy=(difT['Eth'], ytext), xytext=(difT['Eth'], ytext),
                       size='x-large', ha="left", va='bottom', color=color(i))
        for j,r in enumerate(rsel):
            rname = os.path.join(path,name+'_%s'%r)
            try:
                triggE = np.load(rname+triggEffAppendix_E)
            except IOError:
                print "file '%s' not found... let's try the next!"%(rname+triggEffAppendix_E)
                continue
            try:
                triggPE = np.load(rname+triggEffAppendix_pe)
            except IOError:
                print "file '%s' not found... let's try the next!"%(rname+triggEffAppendix_pe)
                continue
            ptr_E[j].errorbar(triggE['E'], triggE['triggProb'], triggE['triggError'], fmt=line(i), mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                              label=r"%s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))
            ptr_p[j].errorbar(triggPE['pe'], triggPE['triggProb'], triggPE['triggError'], fmt=line(i), mec=color(i), mfc=fill(i), marker=mark(i), color=color(i),
                              label=r"%s"%(args['--label'][i] if len(args['--label'])>i else ('input %s'%(i+1))))     
    set_axis(peffA, 'x', r"energy [TeV]", scale='log', scalar_values=True)
    set_axis(peffA, 'y', r"effective area [m$^2$]", scale='log')
    set_axis(pdifT, 'x', r"energy [TeV]", scale='log', scalar_values=True)
    set_axis(pdifT, 'y', r"differential trigger rate (dR/dE) [Hz TeV$^{-1}$]", scale='log')
    lowrightLeg = dict(bbox_to_anchor=(1,0), loc='lower right')
    if len(folders)>1:
        peffA.legend(**lowrightLeg)
        pdifT.legend()
    if outname is not None:
        effAout, difTout = os.path.join(outname,'effArea.pdf'), os.path.join(outname,'diffTRate.pdf')
        feffA.savefig(effAout)
        fdifT.savefig(difTout)
        print "saved '%s'"%effAout
        print "saved '%s'"%difTout
    feffA.show()
    fdifT.show()
    for j,r in enumerate(rsel):
        set_axis(ptr_E[j], 'x', r"energy [TeV]", scale='log', scalar_values=True)
        set_axis(ptr_E[j], 'y', r"trigger efficiency @ r$\leq$%s m"%r)
        set_axis(ptr_p[j], 'x', r"p.e.", scale='log', scalar_values=True)
        set_axis(ptr_p[j], 'y', r"trigger efficiency @ r$\leq$%s m"%r)
        if len(folders)>1:
            ptr_E[j].legend(**lowrightLeg)
            ptr_p[j].legend(**lowrightLeg)
        if outname is not None:
            trEout, trPEout = os.path.join(outname,'TeffE_%s.pdf'%r), os.path.join(outname,'TeffPE_%s.pdf'%r)
            ftr_E[j].savefig(trEout)
            ftr_p[j].savefig(trPEout)
            print "saved '%s'"%trEout
            print "saved '%s'"%trPEout
        ftr_E[j].show()
        ftr_p[j].show()
    raw_input("\n == Press enter to close ==\n")


    

    

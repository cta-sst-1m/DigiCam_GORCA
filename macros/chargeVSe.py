import docopt,os,sys

if __name__=='__main__':
    version = "Version: 0.2 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Make the proton rate

Usage:
  %(prog)s <pe>... [-m MIN -o OUTNAME -R RESFILE... -r REFFILE... -t TITLE -x XMIN]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <pe>                       Npz files with the pe distribution per pixel ('bins' and 'hist' keys)
                             They must be named wtarting with a fixed energy value in TeV, such as
                             1.23TeV_<whatevername>.npz

General options:
  -m --minimum-pe MIN        Minimum number of PE to be considered [default: 0.]
  -o --out-name OUTNAME      Name to save the figure: it will be saved in the <pe> folder
                             [default: ExpectedChargeResol]
  -R --resolution RESFILE    Text file(s) with the charge resolution values ('pe' 'resolution' in the columns).
                             If this option is given multiple times, the name of each file will be used as labels as follow:
                             1) split with '_' and got the last element
                             2) replaced all the '-' as ' '
                             Example: "file-1.txt" will be converte in 'file 1', as well as "test_file-1.txt" 
  -r --reference REFFILE     Text file(s) to be used as reference (same workings of --resolution)
  -t --title TITLE           Title of the plot [default: ]
  -x --xaxis-min XMIN        Minimum value for the x-axis [default: 5.e-2]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)


from chargeRes import *

def do_dQvsE(E, qresTuple, qerrTuple, resolutions, title='', savename=None):
    marks = ['bo', 'rs', 'k^','gv', 'm+', 'c>', 'y<']
    files = [f for F in resolutions for f in glob.glob(F)]
    fig,p = pyp.subplots()
    for i,f in enumerate(files):
        name = os.path.splitext(basename(f))[0]
        label = name.split('_')[-1].replace('-',' ')
        kwd = dict(marker=marks[i][1:], fmt=marks[i][1:], mec=marks[i][0], mfc=marks[i][0], color=marks[i][0], label=r"%s"%label, lw=1.)
        p.errorbar(E,qresTuple[i],qerrTuple[i],**kwd)
    fig.suptitle(r"%s"%title, size="xx-large", weight='bold')#"Expected distribution of the charge resolution"
    set_axis(p, 'x', r'energy [TeV]', scale='log', scalar_values=True)
    set_axis(p, 'y', r'$\mathbf{\sigma_Q/Q}$', lims=(1e-1,.8))
    if len(files)>1:
        p.legend()#bbox_to_anchor=(0.,1.), loc="upper left")
    if savename is not None:
        fig.savefig(outname+".eps")
    fig.show()


if __name__=='__main__':
    pefiles = [f for F in args['<pe>'] for f in glob.glob(F)]
    resolutions = args['--resolution']
    references = args['--reference']
    min_pe = float(args['--minimum-pe'])
    outname = args['--out-name']
    xmin = float(args['--xaxis-min'])
    title = lambda txt: '%s'%txt if args['--title'] in [None, '', 'None', 'none'] else args['--title']+' (%s)'%txt
    energies, dQ_means, dQ_uncert = [], [], []
    for pef in pefiles:
        namefile = basename(pef)
        energy = namefile.split('_')[0]
        if energy[-3:]!='TeV':
            print " *** the file must have as tructure as '<energy>TeV_<whatevername>.npz' to be suited for this script!"
            sys.exit(130)
        energies += [ float(energy[:-3]) ]
        fig, averages, ave_errs = makePlot(pef, resolutions, references, min_pe, True,
                                      os.path.join(dirname(aspath(pef)), basename(outname)+"_%s.eps"%energy), xmin, title(energy))
        raw_input(" -- Press enter to continue --")
        pyp.close(fig)
        dQ_means += [averages]
        dQ_uncert += [ave_errs]
    energies, dQ_means, dQ_uncert = np.array(energies), np.swapaxes(dQ_means, 0,1), np.swapaxes(dQ_uncert,0,1)
    do_dQvsE(energies, dQ_means, dQ_uncert, resolutions,
             title="%s $\leq$ E [TeV] $\leq$ %s (pe $\geq$ %s)"%(int(energies.min()),int(energies.max()), int(min_pe)), savename=outname)
    raw_input(" => Press enter to close... <=")

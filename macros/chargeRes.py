import docopt,os,sys

if __name__=='__main__':
    version = "Version: 0.2 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Make the proton rate

Usage:
  %(prog)s <pe> <resolutions>... [-a -m MIN -o OUTNAME -r REFFILE... -t TITLE -x XMIN]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Arguments:
  <pe>                       Npz file with the pe distribution per pixel ('bins' and 'hist' keys)
  <resolutions>              Text files with the charge resolution values ('pe' 'resolution' in the columns).
                             If multiple files are given, the name of each file will be used as labels as follow:
                             1) split with '_' and got the last element
                             2) replaced all the '-' as ' '
                             Example: "file-1.txt" will be converte in 'file 1', as well as "test_file-1.txt" 

General options:
  -a --average               Calculate the weighted average of the charge resolution
  -m --minimum-pe MIN        Minimum number of PE to be considered [default: 0.]
  -o --out-name OUTNAME      Name to save the figure: it will be saved in the <pe> folder
                             [default: ExpectedChargeResol]
  -r --reference REFFILE     Text file(s) to be used as reference (same workings of <resolutions>)
  -t --title TITLE           Title of the plot [default: '']
  -x --xaxis-min XMIN        Minimum value for the x-axis [default: 5.e-2]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from root2npz import *

def makePlot(pefile, resolutions, references, min_pe=0., do_average=False, outname="ExpectedChargeResol", xmin=5.e-2, title=''):
    with np.load(aspath(pefile)) as l:
        bins, hist = l['bins'], l['hist']
    entries = float(hist[np.where(bins[:-1]>=min_pe)[0]].sum())
    files = [f for F in resolutions for f in glob.glob(F)]
    marks = ['bo', 'cs', 'rv','g^', 'm+', 'c>', 'y<']#['bo', 'rs', 'k^','gv', 'm+', 'c>', 'y<']
    averages, ave_errs = [], []
    fig,p = pyp.subplots()
    for i,f in enumerate(files):
        if i>=len(marks):
            print "No more resolutions files are accepted!"
            break
        print("resolution file: %s"%f)
        #photoelectrons, resolutions, pe_error, res_error
        data = np.loadtxt(aspath(f), unpack=True)
        len_data = len(data)
        if len_data==2:
            photoelectrons, resolutions, pe_error, res_error = data[0], data[1], None, None
        elif len_data==3:
            photoelectrons, resolutions, pe_error, res_error = data[0], data[1], None, data[2]
        elif len_data==4:
            photoelectrons, resolutions, pe_error, res_error = data[0], data[1], data[2], data[3]
        elif len_data==5:
            photoelectrons, resolutions, pe_error, res_error = data[0], data[1], data[2], data[3:5]
        elif len_data==6:
            photoelectrons, resolutions, pe_error, res_error = data[0], data[1], data[2:4], data[4:6]
        else:
            raise AssertionError("The file '%s' contains %s columns: at least 2 and no more than 6 are expected"%(f,len_data))
        count, qreslist, qerrlist, ibin = [], [], [], []
        for pe in photoelectrons[photoelectrons>=min_pe]:
            iselect = np.where(bins<pe)[0][-1]
            if iselect>=len(hist):
                break
            ibin += [iselect]
            count += [hist[iselect]]
            index = np.where(photoelectrons==pe)[0][0]
            qreslist += [resolutions[index]]
            qerrlist += [res_error[index]] if len_data<5 else [[res_error[0][index],res_error[0][index]]]
        isort = np.argsort(ibin)
        ibin, count, qreslist, qerrlist = np.array(ibin)[isort], np.array(count)[isort], np.array(qreslist)[isort], np.array(qerrlist)[isort]
        name = os.path.splitext(basename(f))[0]
        label = name.split('_')[-1].replace('-',' ')
        ibin, iunique = np.unique(ibin,True)
        count = count[iunique]
        iunique = np.append(iunique,len(qreslist))
        qres, qerr = [], []
        if len_data<5:
            qerr = []
        else:
            qerrl = []
            qerrh = []
        for ind in xrange(1,len(iunique)):
            qres += [qreslist[iunique[ind-1]:iunique[ind]].mean()]
            if len_data<5:
                qerr += [np.sqrt((qerrlist[iunique[ind-1]:iunique[ind]]**2).sum()/(iunique[ind]-iunique[ind-1]))]
            else:
                qerrl += [np.sqrt((qerrlist[0][iunique[ind-1]:iunique[ind]]**2).sum()/(iunique[ind]-iunique[ind-1]))]
                qerrh += [np.sqrt((qerrlist[1][iunique[ind-1]:iunique[ind]]**2).sum()/(iunique[ind]-iunique[ind-1]))]
        qres, qerr = np.array(qres), np.array(qerr if len_data<5 else [qerrl,qerrh])
        eventRatio = count*100./entries
        if do_average:
            w = eventRatio/np.sum(eventRatio)
            average = np.sum(w*qres)
            ave_err = np.sqrt(np.sum(w*qerr**2))
            averages += [average]
            ave_errs += [ave_err]
            print "\t>>> %s => average charge resolution: %1.2f +- %1.2f <<<" %(label, average, ave_err)
            label += ' (mean = %1.2f$\pm$%1.2f)' % (average,ave_err)
        kwd = dict(marker=marks[i][1:], fmt=marks[i][1:], mec=marks[i][0], mfc=marks[i][0], color=marks[i][0], label=r"%s"%label, lw=1.)
        p.errorbar(count*100./entries,qres,qerr,**kwd)
    refs = [r for R in references for r in glob.glob(R)]
    lines = ['r-', 'k--', 'b-.','g:', 'm-', 'co', 'ys']
    for i,r in enumerate(refs):
        if i>=len(lines):
            print "No more reference files are accepted!"
            break
        print "reference file:", f
        photoelectrons, resolutions = np.loadtxt(aspath(r), unpack=True)
        count, qreslist, ibin = [], [], []
        for pe in photoelectrons:#[photoelectrons>float(args['--minimum-pe'])]:
            iselect = np.where(bins<pe)[0][-1]
            if iselect>=len(hist):
                break
            ibin += [iselect]
            count += [hist[iselect]]
            qreslist += [resolutions[np.where(photoelectrons==pe)[0][0]]]
        isort = np.argsort(ibin)
        ibin, count, qreslist = np.array(ibin)[isort], np.array(count)[isort], np.array(qreslist)[isort]
        ibin, iunique = np.unique(ibin,True)
        count = count[iunique]
        iunique = np.append(iunique,len(qreslist))
        qres = []
        for ind in xrange(1,len(iunique)):
            qres += [qreslist[iunique[ind-1]:iunique[ind]].mean()]    
            name = os.path.splitext(basename(r))[0]
        label = name.split('_')[-1]
        kwd = dict(marker='', fmt=lines[i][1:], mec=lines[i][0], mfc=lines[i][0], color=lines[i][0], label=label.replace('-',' '))
        p.errorbar(count*100./entries,qres, **kwd)
    
    fig.suptitle(r"%s"%title, size="xx-large", weight='bold')#"Expected distribution of the charge resolution"
    set_axis(p, 'x', r'fraction of pixels (pe$\geq$'+str(int(min_pe))+r') [%]', lims=(xmin,100.), scale='log', scalar_values=True)
    set_axis(p, 'y', r'$\mathbf{\sigma_Q/Q}$', lims=(3e-2,4.), scale='log')
    if len(files)>1 or len(refs)>0:
        p.legend(bbox_to_anchor=(0.,1.), loc="upper left")
    print outname
    fig.savefig(outname)
    fig.show()
    if do_average:
        return fig, averages, ave_errs
    else:
        return fig, None, None

if __name__=='__main__':
    pefile = args['<pe>']
    resolutions = args['<resolutions>']
    references = args['--reference']
    min_pe = float(args['--minimum-pe'])
    do_average = args['--average']
    outname = args['--out-name']
    if outname is not None:
        outname = os.path.join(dirname(aspath(pefile)), outname+(".eps" if os.path.splitext(outname)[1]=='' else ''))
    xmin = float(args['--xaxis-min'])
    title = args['--title']
    fig, _, _ = makePlot(pefile, resolutions, references, min_pe, do_average, outname, xmin, title)
    raw_input(" => Press enter to close... <=")
    pyp.close(fig)


        

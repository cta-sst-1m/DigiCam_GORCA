#!/usr/bin/env python

if __name__=='__main__':
    import docopt,os,sys
    #import sys, os, glob
    #from root_numpy import tree2rec

    version = "Version: 0.2 (A. Porcelli, ale_led@yahoo.it)"
    doc = """
Make effective area

Usage:
  %(prog)s (effarea|photons) <input>... [-o OUTPUT -V...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Command:
  effarea                    It will create and save a data set for the effective area calculation: lge, r, r_triggered
  info                       It will create and save a data set for the ph/PE info: lge, tiggeredTel, iPhotons, pixelPE
Arguments:
  <input>                    Input file(s). Format compatible to TChain for multiple parsing.

General options:
  -o --output OUTPUT         If defined, this name will be used to save the npz file, otherwise the <input> is used
  -V --verbose               Level of verbose output

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

from commons import *

# Data formatting
def build_effArea_variables(datafile, savenpz=None, verbose=2):
    """
    It reads the root files and return the array formatted for CameraView
    
    Parameters
    ----------
    datafile: string or list of string
        root file name(s) in the TChain format
    verbose: int (default: 0)
        level of verbosity
    savez: string (default: None)
        name of the npz file to save the dataset
    
    Returns
    ----------
    lge: ndarray
        log10 of the energy values
    r: ndarray
        radius of the events
    r_triggered: ndarray
        as <r>, but the not triggered events have a -1 value
    """
    if not isiterable(datafile):
        datafile = [datafile]
    if datafile[0][-4:]=='.npz':
        if verbose>1:
            print("datafile (npz compressed): %s" % ', '.join(datafile))
        return get_npz_elements(datafile, 'lge', 'r', 'r_triggered', 'pe', 'ph', verbose=verbose)
    else:
        if verbose>1:
            print("datafile (ROOT::TChain format compatible): %s" % ', '.join(datafile))
        tree = R.TChain("Events/tSimulatedEvents")
        tTel = R.TChain("Events/T0")
        try:
            for df in datafile:
                tree.Add(df)
                tTel.Add(df)
        except:
            print("TChain failed: try with npz files")
            return get_npz_elements(datafile, 'lge', 'r', 'r_triggered', 'pe', 'ph', verbose=verbose)
        else:
            lge, r, r_triggered, pe, ph, i = [], [], [], [], [], 0
            if verbose:
                print("tSimulatedEvents: chained root files")
                loadstep = int(tree.GetEntries()/100)
                progress = AnimatedProgressBar(end=100, width=100)
                print("Mapping tSimulatedEvents tree:")
                progress.progressing(0)
            for t in tree:
                lge += [np.log10(t.energy)]
                r += [np.sqrt(t.xcore**2+t.ycore**2)]
                r_triggered += [r[-1]  if get_branch_value(t.arrayTriggerBit) else -1]
                tTel.GetEntry(i)
                pe += [np.sum(tTel.vPEInPixel)]
                ph += [tTel.iPhotonsInFocalPlane]
                i += 1
                if verbose and i%loadstep==0:
                    progress.progressing(1, ('%s event'%cardinal(i)) if verbose>1 else None)
            lge, r, r_triggered, pe, ph = np.array(lge), np.array(r), np.array(r_triggered), np.array(pe), np.array(ph)
            infs = np.isinf(lge)
            lge, r, r_triggered, pe, ph = lge[~infs], r[~infs], r_triggered[~infs], pe[~infs], ph[~infs]
            if verbose:
                progress.finish("%s simulated events mapped"%i if verbose>1 else None)
            if savenpz is not None:
                np.savez(savenpz, lge=lge, r=r, r_triggered=r_triggered, pe=pe, ph=ph)
                if verbose:
                    print " ...saved varabiales as '%s'"%savenpz
            return lge, r, r_triggered, pe, ph

def get_branch_value(value):
    """Convert the branch <value> into ndarray compatible data"""
    if isinstance(value, str):
        return True if value=='\x01' else (False if value=='\x00' else value)
    elif hasattr(value, '__len__'):
        return [get_branch_value(value[i]) for i in xrange(len(value))]
    elif repr(value).find('<bool>')>-1:
        return bool(value)
    else:
        return value

def get_npz_elements(datafile, *keys, **kwd):
    if len(keys)==0:
        raise StandardError('No data are specified to be retrived... Are you sure this is what you want???')
    verbose = kwd.get('verbose',2)
    filelist = sorted([f for df in datafile for f in glob.glob(aspath(df)) if f[-4:]=='.npz'])
    if len(filelist)==0:
        raise StandardError('No npz files are found')
    if verbose:
        steps = len(filelist)
        print("Recovering data from files:")
        if steps>1:
            progress = AnimatedProgressBar(end=steps, width=steps)
            progress.progressing(0)
        else:
            print " Single file: %s"%filelist[0],
    data = [[] for i in xrange(len(keys))]
    for iF,f in enumerate(filelist):
        with np.load(f) as l:
            for ik,k in enumerate(keys):
                try:
                    data[ik] += list(l[k])
                except KeyError:
                    continue
        if verbose:
            if steps>1:
                progress.progressing(1, ('%s file loaded (%s)'%(cardinal(iF+1),os.path.basename(f))) if verbose>1 else None)
            else:
                print("...loaded!")
    for i,d in enumerate(data):
        if len(d)==0:
            data[i] = None
        else:
            data[i] = np.array(d)
    if verbose and steps>1:
        progress.finish("%s npz files from simulated events loaded"%(iF+1) if verbose>1 else None)
    return data
    
def build_photonInfo(datafile, savenpz=None, verbose=2):
    """
    It reads the root files and return the array formatted for CameraView
    
    Parameters
    ----------
    datafile: string or list of string
        root file name(s) in the TChain format
    verbose: int (default: 0)
        level of verbosity
    savez: string (default: None)
        name of the npz file to save the dataset
    
    Returns
    ----------
    lge: ndarray
        log10 of the energy values
    tiggeredTel: 2d ndarray
        triggered telescopes: (<ievent>,<itel>)
    iPhotons: 2d ndarray
        photons reaching the camera in each telescope: (<ievent>,<itel>)
    pixelPE: 3d ndarray
        number of PE in each pixel in each telescope: (<ievent>,<itel>,<ipixel>)
    """
    if not isiterable(datafile):
        datafile = [datafile]
    if datafile[0][-4:]=='.npz':
        if verbose>1:
            print("datafile (npz compressed): %s" % ', '.join(datafile))
        return get_npz_elements(datafile, 'lge', 'tiggeredTel', 'iPhotons', 'pixelPE', verbose=verbose)
    else:
        if verbose>1:
            print("datafile (ROOT::TChain format compatible): %s" % ', '.join(datafile))
        tSim = R.TChain("Events/tSimulatedEvents")
        try:
            for df in datafile:
                tSim.Add(df)
        except:
            print("TChain failed: try with npz files")
            return get_npz_elements(datafile, 'lge', 'tiggeredTel', 'iPhotons', 'pixelPE', verbose=verbose)
        else:
            lge, tiggeredTel, ntels, ievent = [], [], [], 0
            if verbose:
                print("tSimulatedEvents: chained root files")
                loadstep = int(tSim.GetEntries()/100)
                progress = AnimatedProgressBar(end=100, width=100)
                print("Mapping tSimulatedEvents tree:")
                progress.progressing(0)
            for t in tSim:
                lge += [np.log10(t.energy)]
                tiggeredTel += [[bool(t.vTelescopeTriggerBits[i]) for i in xrange(len(t.vTelescopeTriggerBits))]]
                ntels += [t.uNumTelescopes]
                ievent += 1
                if verbose and ievent%loadstep==0:
                    progress.progressing(1, ('%s event'%cardinal(ievent)) if verbose>1 else None)
            lge, tiggeredTel = np.array(lge), np.array(tiggeredTel)
            if verbose:
                progress.finish("%s simulated events mapped"%ievent if verbose>1 else None)
            tels = np.array(ntels).max()
            iPhotons, pixelPE = [], []
            for tel in xrange(tels):
                tTel =  R.TChain("Events/T%s"%tel)
                iPhotons += [[]]
                pixelPE += [[]]
                itevent = 0
                for df in datafile:
                    tTel.Add(df)
                if verbose:
                    print("tSimulatedEvents: chained root files")
                    loadstep = int(tSim.GetEntries()/100)
                    progress = AnimatedProgressBar(end=100, width=100)
                    print("Mapping T%s tree:"%tel)
                    progress.progressing(0)
                for t in tTel:
                    iPhotons[-1] += [t.iPhotonsInFocalPlane]
                    pixelPE[-1] += [list(t.vPEInPixel)]
                    itevent += 1
                    if verbose and itevent%loadstep==0:
                        progress.progressing(1, ('%s event'%cardinal(itevent)) if verbose>1 else None)
                if verbose:
                    progress.finish("%s T%s events mapped"%(itevent,tel) if verbose>1 else None)
            iPhotons =  np.swapaxes(iPhotons, 0, 1)
            pixelPE = np.swapaxes(pixelPE, 0, 1)
            infs = np.isinf(lge)
            lge, tiggeredTel, iPhotons, pixelPE = lge[~infs], tiggeredTel[~infs], iPhotons[~infs], pixelPE[~infs]
            if savenpz is not None:
                np.savez(savenpz, lge=lge, tiggeredTel=tiggeredTel, iPhotons=iPhotons, pixelPE=pixelPE)
                if verbose:
                    print " ...saved varabiales as '%s'"%savenpz
            return lge, tiggeredTel, iPhotons, pixelPE

def from_npz(file_name):
    with np.load(file_name) as l:
        variables = {k:a for k,a in l.iteritems()}
    return variables


if __name__=='__main__':
    files = sorted([f for ifile in args['<input>'] for f in glob.glob(aspath(ifile)) if f[-5:]=='.root'])
    verbose = args['--verbose']
    if len(files)==0:
        print "No good root files are given"
    for i,f in enumerate(files):
        if args['--output'] is None:
            outname = os.path.splitext(f)[0]+'_%s.npz'%sys.argv[1]
        else:
            outname = args['<output>']+'_{}_{:06d}'.format(sys.argv[1],i)+'.npz'
        if args['effarea']:
            results = build_effArea_variables(f, savenpz=outname, verbose=verbose)
        elif args['photons']:
            results = build_photonInfo(f, savenpz=outname, verbose=verbose)
        else:
            raise SystemExit("Option not yet implemented") #Not yet implemented
    if verbose>2:
        for r in results:
            print r

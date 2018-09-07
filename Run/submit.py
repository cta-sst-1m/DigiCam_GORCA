#!/usr/bin/env python
from parser import *

start, skip, verbosity = get_value('--start', int), get_value('--skip', int), args["--verbose"]
localenv = args["--localset"]
os.environ['thislocalsource'] = localenv
if verbosity:
    print("\n  --> Set environments from %s\n" % localenv)
bash_source(localenv, verbose=verbosity>1)

groptics_extras = lambda: [get_value('--defaults'), get_value('--cfg-gro'), get_value('--telarray'), str(int(args['--use-no-shadowing']==False)), get_value('--atmosphere')]
file2outname = lambda fname: os.path.split(os.path.normpath(os.path.splitext(fname)[0]))[-1]
path2outname = lambda fname: os.path.split(os.path.normpath(fname))[-1]
outname = lambda fname: file2outname(fname) if os.path.isfile(fname) else path2outname(fname)
runext = lambda i: "Run{:>06d}".format(i)
filelist = lambda paths=args['<input>']: sorted([f for p in paths for f in glob.glob(os.path.expandvars(os.path.os.path.expanduser(p)))])[skip:]

batchlog = lambda name: os.path.join(os.environ["out_dir"],args["<output>"],"Logs","%s.batchlog" % name)
dryrun_id = lambda name: ("%s.dryrun" % name) if args["--dryrun"] else None

if len(args['<input>'])>0:
    flist = filelist()
    runs = len(flist) if args['--runs'] is None else int(args['--runs'])
    iterator = (fname for fname in flist)
else:
    runs = 1 if args['--runs'] is None else int(args['--runs'])
    iterator = (None for i in xrange(runs))

def care_extras():
    extra = []
    for e in args['--extras']:
        e = e.split(" ")
        extra += [e[0] + (' \\\\\\\"\"\"'+' '.join(e[1:])+'\\\\\\\"\"\"' if len(e)>1 else '')]
    extra = [get_value('--defaults'), get_value('--cfg-care'), get_value('--telCARE')]+extra
    return extra
def care_others(extras):
    other = ['--notraces' if args['--notraces'] else "", "--writepedestals 1" if args['--writepedestals'] else ""]
    other = " ".join([e for e in other if e!=""])
    if len(other)>0:
        extras += ["'%s'" % other]
    return extras

if args[care]:
    extra = care_extras()
    if len(args['<input>'])>0:
        code = exe_chosen = run_care
        exe = os.path.join(os.environ['runPath'],code+'.sh')
        extra = care_others(extra)
        cmd = lambda i,fname,oname: " ".join([exe,args['<output>'],str(i),
                                              os.path.join(fname,"photonLocation.root") if os.path.isdir(fname) else fname,
                                              oname]+extra)
    else:
        code = exe_chosen = run_biascurve
        exe = os.path.join(os.environ['runPath'],code+'.sh')
        cmd = lambda i,fname,oname: " ".join([exe,args['<output>'],str(i),oname]+extra)

elif args[groptics]:
    code = exe_chosen = run_groptics
    exe = os.path.join(os.environ['runPath'],code+'.sh')
    cmd = lambda i,fname,oname: " ".join([exe,args['<output>'],str(i),fname,oname]+groptics_extras())

elif args[both]:
    exe_chosen = run_groptics+"+"+run_care
    code = both
    exe1 = os.path.join(os.environ['runPath'],run_groptics+'.sh')
    exe2 = os.path.join(os.environ['runPath'],run_care+'.sh')
    extra = care_others(care_extras())
    grout = lambda oname: os.path.join(os.environ["out_dir"],args["<output>"],'GrOptics',oname,'photonLocation.root')
    cmd = lambda i,fname,oname: " ".join([exe1,args['<output>'],str(i),fname,oname]+groptics_extras())+"; "+\
      " ".join([exe2,args['<output>'],str(i),grout(oname),oname]+extra)

if verbosity:
    print("  --> Executable chosen: %s\n" % exe_chosen)


print(" ==> Submitting...\n")
run = 0
qsub_kwd = dict(queue=os.environ["queue"] if args['--queue'] is None else args['--queue'], cwd=os.environ['runPath'],
                verbose=verbosity, PBS=True, mem=os.environ["mem"], vmem=os.environ["vmem"], walltime=os.environ["walltime"])
for fname in iterator:
    if run == runs:
        if verbosity:
            print("   ...the the inputs are more than the number of job requested, so the submitting ends here\n")
        break
    if fname is not None and os.path.basename(fname[:-1] if fname[-1]=='/' else fname) in ignorefolders:
        continue
    runID = run+start
    oname = outname(runext(runID) if fname is None else fname)
    name = code+'_'+oname
    qsub(cmd(runID,fname,oname), name=name, log=batchlog(name), dryrun_id=dryrun_id(name), **qsub_kwd)
    run += 1

print(" ==> %s job%s been submitted\n" % (("All "+str(run), "s have") if run>1 else (run, " has")))

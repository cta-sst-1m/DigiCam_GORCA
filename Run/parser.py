from shell_utils import *

version = "Version: 0.5 (A. Porcelli, ale_led@yahoo.it)"
version = "Version: 0.6 (I. Al Samarai)"
prog = os.path.basename(sys.argv[0])
thisdir = os.path.abspath(os.path.dirname(sys.argv[0]))
care, groptics, both = 'care', 'groptics', 'gro+care'
run_care, run_groptics, run_biascurve = 'care', 'groptics', 'biascurve'
doc = """
Submit to the batch PBS the simulations of the biascurve, CARE, grOptics or the chain 'grOptics+CARE'

Usage:
  {prog} {groptics} <output> <input>... [-A ATMOSPHERE -C GR_CFG -d SOURCE -l LOCAL -q QUEUE -r RUNS -S SKIP -s START -T TELARRAY -DU] [-V...]
  {prog} {care} <output> [<input>...] [-c CARE_CFG -d SOURCE -l LOCAL -q QUEUE -r RUNS -S SKIP -s STRAT_RUN -t TELCARE -Dnw] [-V...] [-e EXTRAS...]
  {prog} {both} <output> <input>... [-A ATMOSPHERE -C GR_CFG -c CARE_CFG -d SOURCE -l LOCAL -q QUEUE -r RUNS -S SKIP -s START -T TELARRAY -t TELCARE -DnwU] [-V...] [-e EXTRAS...]

  {prog} (-h | -v)

Commands:
  {care:<30}  Run {run_care}.sh (camera simulation), or biascurve.sh if no <input_folder> is given
  {groptics:<30}  Run {run_groptics}.sh (photon ray-tracing in the telescope) form coriska input
  {both:<30}  Do '{groptics}', than '{care}' from the result just obtained (all options available!)

Info options:
  -h --help                       Print this help
  -v --version                    Print the version

Arguments:
  <output>                        Name of the output, used to build the folder(s) where the outputs will be saved. The file name
                                  will be automatically built with the <input> argument if given, or with the run number otherwise
  <input>                         CORISKA file(s) for '{groptics}' or '{both}' commands; folder(s) where to find
                                  'photonLocation.root', or directly the GrOptics file(s), for '{care}' command.
                                  If not <input> is given in case of '{care}' command, the bias curve
                                  ('{run_biascurve}.sh') will be simulated.

Common options:
  -D --dryrun                     Execute a dryrun (only printout, no real execution)
  -d --defaults SOURCE            File where default parameter settings are
  -l --localset LOCAL             Paths of the extra environament [default: {thisdir}/local.sh]
  -q --queue QUEUE                Queue name (default set by local.sh or defined LOCAL with -l option)
  -r --runs RUNS                  Number of runs. Must be and int! If no <input> is given (only '{care}' command!)
                                  the default is 1, else it is the lenght of the expanded <input> argument
  -S --skip SKIP                  Number of <input>s to skip before start the runs. Useless in no <input is given>.
                                  Must be and int! [default: 0]
  -s --start START                Number of starting run. Must be and int! [default: 0]
  -V --verbose                    Verbose stdout. Multiple uses increase the verbosity level

'{groptics}' options:
  -A --atmosphere ATMOSPHERE      Atmosphere config file for corsikaIOreader beyond the default
  -C --cfg-gro GR_CFG             Path of the configuration file for GrOptics  (no extensions!)
  -T --telarray TELARRAY          Path of the telescope array configuration file for GrOptics (it's a python file!)
  -U --use-no-shadowing           To use no shadowing

'{care}' options:
  -c --cfg-care CARE_CFG          Path of the configuration file for CARE
  -e --extras EXTRAS...           Extra options to change the configuration parameters interactively. Each '-e' is a
                                  parameter to change with its change enquoted: the same example in the CARE tutorial
                                  will be declared with '-e "NSBRATEPERPIXEL 0 120"'
  -n --notraces                   Do not save the traces in the output root file (see CARE tutorial)
  -t --telCARE TELCARE            Telescope setting for CARE beyond the default
  -w --writepedestals             Simulate also the pedestal (see CARE tutorial)

{version}
Enjoy!
""".format(prog=prog, thisdir=thisdir, version=version, care=care, groptics=groptics, both=both,
           run_care=run_care, run_groptics=run_groptics, run_biascurve=run_biascurve)

args = docopt.docopt(doc, help=True, version=version)

def get_value(arg, this_type=str, arg_dict=args):
    if this_type==str:
        return '\\\"\\\"' if arg_dict[arg] is None else arg_dict[arg]
    try:
        value = this_type(arg_dict[arg])
    except ValueError:
        print("\n*** '%s' is not %s ***\n" % (arg, typef))
        docopt.docopt(doc, argv=[''], help=True, version=version)
        raise SystemExit
    else:
        return value

ignorefolders = ['Bads','BiasCurve','CARE','corsikaIOreader','GrOptics','Logs']
    
#def get_none_default(value):
#    return '' if value is None else value


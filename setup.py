#!/usr/bin/env python

# import modules and tools
try:
    from Run.shell_utils import *
except ImportError:
    print("ERROR: no Run.shell_utils module is found. Check what is wrong...")
    sys.exit(1)

# stdin browsable
try:
    import readline
    has_browsing = True
except ImportError:
    has_browsing = False
    print("warning: no readline module is found. CLI browsing not available in stdin")

# Colored terminal
try:
    import colorama
    has_colorama = True
except ImportError:
    has_colorama = False
    #print("info: if you install the colorama module, the outout will be fancier!")
if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
    if has_colorama:
        BCYAN = colorama.Style.BRIGHT+colorama.Fore.CYAN
        BBLACK = colorama.Style.BRIGHT+colorama.Fore.BLACK
        BBLUE = colorama.Style.BRIGHT+colorama.Fore.BLUE
        BRED = colorama.Style.BRIGHT+colorama.Fore.RED
        RESET = colorama.Style.RESET_ALL if has_colorama else ''
    else:
        BCYAN, BBLACK, BBLUE, BRED, RESET = '\x1b[1;36m','\x1b[1;30m','\x1b[1;34m','\x1b[1;31m','\x1b[0m'
else:
    BCYAN, BBLACK, BBLUE, BRED, RESET = '','','','','' 

# globals
get_input = lambda text: raw_input(text) if sys.version_info.major==2 else input(text)

prog = os.path.basename(sys.argv[0])
maindir = abspath(os.path.dirname(sys.argv[0]))
browsing="o The input line has CLI browsing: it is possible to go back in the line, catch previous lines, etc.\n"
rules = """
 RULES:
o The inserted paths (for the parent folders for output and temporary data, and root location)
  can be local (this/exmaple is the same of ./this/example) or absolute (/My/Home/this/example)
o sh variables can be used but they will be not expanded in the generated files. Example: $HOME/$USER/this/example
{browsing}
 NOTE:
The main path will be automatically located as the same of {prog}
""".format(prog=prog, browsing=browsing if has_browsing else '')
version="v1.0, environmental python v{major}.{minor} (A. Porcelli, ale_led@yahoo.it)".format(major=sys.version_info.major,
                                                                                              minor=sys.version_info.minor)
usage = """
Usage:
  python {prog} [ help | version | install ]
  python {prog} [ -h | -v | -i ]
  python {prog} [ --help | --version | --install ]
or, after a 'chmod +x {prog}',
  ./{prog} [ h | v | i ]
  ...
""".format(prog=prog)
doc = """
Generate the local configuration file(s) starting from the .ini in the Run folder.
If 'install' command is given, also an automatic installation of the software occures.
{usage}
Possible options (case insensitive):
  help -h --help          Print this help and exit
  version -v --version    Print the version and exit

  install -i --install    Install all the softwares. If no commands are given, only the local setting
                          configurations are generated
{rules}
Now, configure!

Version: {version}
""".format(usage=usage,rules=rules,version=version)

# printout format
title = lambda msg: BCYAN + msg + RESET
info = lambda msg: BBLACK + msg + RESET
state = lambda msg: BBLUE + msg + RESET
error = lambda msg:  BRED + msg + RESET
stateloc = lambda loc,msg: state("[%s]"%loc)+" "+msg
errorloc = lambda loc,msg: error("[%s]"%loc)+" "+msg

# building function
def build(what, *execute_pars, **execute_kwd):
    print("\n"+stateloc(what,"Building...\n"))
    try:
        output = execute(*execute_pars, **execute_kwd)
    except sp.CalledProcessError as e:
        print(e)
        print("\n"+errorloc(what,"...FAILED! Install it manually or check what is going on and retry\n"))
        return 1
    else:
        if output is not None:
            print(output)
        print("\n"+stateloc(what,"...DONE!\n"))
        return 0

commands_dict = dict(grogit=('GrOptics git', 'git clone http://www.gtlib.gatech.edu/pub/IACT/GrOptics.git'),
                     caregit=('CARE git', 'git clone http://www.gtlib.gatech.edu/pub/IACT/CARE.git'),
                     iotools=('corsikaSimulationTools', 'tar -xzf corsikaSimulationTools.tar.gz; cd corsikaSimulationTools; make'),
                     groptics=('{groptics}', 'cd {groptics}; make'),
                     linking=('data link', 'ln -s %s/corsikaSimulationTools/data {groptics}/data'%maindir),
                     vbf=('VBF', 'cd %s/{care}; tar -xzvf  VBF-0.3.4.tar.gz; cd VBF-0.3.4; ./configure --prefix=$PATHTOVBFLIBDIRECTORY; makes; make install'%maindir),
                     care=('{care}', 'cd {care}; make'),
                     makefile=('alterative Makefile for {software}','cp Makefile{software}_alternative {software}/Makefile')
                     )
command_seq = ['grogit', 'caregit', 'iotools', 'groptics', 'linking', 'vbf', 'care']


if __name__=='__main__':
    install = False
    if len(sys.argv)>1:
        opt = sys.argv[1].lower().replace('-','')
        if opt in ['h','help']:
            print(doc)
            sys.exit(0)
        elif opt in ['v','version']:
            print(version)
            sys.exit(0)
        elif opt in [ 'install', 'i' ]:
            install = True
        else:
            print(info("\n'%s': unknown command"%opt))
            print(usage)
            sys.exit(1)
    print(title("\n ==> Welcome to the package setup! <==\n"))
    print(state("   == Start the local configuration ==\n")+rules)
    print(" --> Main path idetified as '%s'\n"%maindir)
    outpath = abspath(get_input("  * Parent location for the outputs? [%s] "%os.environ['PWD']))
    print(" --> Output path defined as '%s'\n"%outpath)
    tmppath = abspath(get_input("  * Parent location for the temporary files? [%s] "%os.environ['PWD']))
    print(" --> Temporary path defined as '%s'\n"%tmppath)
    rootdef = os.environ.get('ROOTSYS')
    if rootdef is None:
        print("WARNING: no ROOTSYS is found in the system. If root is not installed, do it and set it BEFORE compiling the softwares!")
        rootdef = os.environ['HOME']
    rootsys = get_input("  * Where is located your root? [%s] "%abspath(rootdef))
    rootsys = abspath(rootdef) if rootsys=='' else abspath(rootsys)
    print(" --> ROOTSYS denfined as '%s'\n"%rootsys)
    groptics = get_input("  * How is named your default GrOptics? [GrOptics] ")
    groptics = 'GrOptics' if groptics=='' else groptics
    print(" --> GrOptics installation is named '%s'\n"%groptics)
    care = get_input("  * How is named your default CARE? [CARE] ")
    care = 'CARE' if care=='' else care
    print(" --> CARE installation is named '%s'\n"%care)
    print(info(" = Make output and temporary folders (if they don't exists already)"))
    mkdirs(outpath,tmppath)
    print("")
    initfiles = glob.glob(os.path.join(maindir,"Run/*.sh.ini"))
    for ifile in initfiles:
        with open(ifile) as f:
            text = f.read()
        text = text.replace('@MAINDIR@',maindir).\
          replace('@OUTPUT@',outpath).\
          replace('@TMP@',tmppath).\
          replace('@GROPTICS@',groptics).\
          replace('@CARE@',care).\
          replace('@ROOTSYS@',rootsys)
        ofile = os.path.splitext(ifile)[0]
        with open(ofile,'w') as f:
            f.write(text)
        print(info(" = Generated '%s' (from '%s')"%(ofile,ifile)))
    print(state("\n   == ...the local configuration is over =="))
    if install:
        print(state("\n   == Start installation ==\n"))
        os.environ['PATHTOVBFLIBDIRECTORY'] = os.path.join(maindir,care,'VBF-0.3.4/lib')
        os.environ['LD_LIBRARY_PATH'] = ':'.join([os.environ['LD_LIBRARY_PATH'], os.path.join(os.environ['PATHTOVBFLIBDIRECTORY'],'lib')])
        os.environ['PATH'] = ':'.join([os.path.join(os.environ['PATHTOVBFLIBDIRECTORY'],'bin'), os.environ['PATH']])
        for cmd in command_seq:
            what, command = commands_dict[cmd]
            what = what.format(groptics=groptics, care=care)
            command = command.format(groptics=groptics, care=care)
            if (what=='GrOptics git' and (groptics!='GrOptics' or path_exists(groptics))) or\
              (what=='CARE git' and (care!='CARE' or path_exists(care))) or\
              (cmd=='linking' and path_exists(command.split()[-1])) or\
              (what==groptics and not path_exists(groptics)) or\
              (what==care and not path_exists(care)):
                if cmd=='linking':
                    print("\n"+errorloc(what,"Link '%s' alredy existing: link building skipped\n"%command.split()[-1]))
                continue
            ret = build(what, command, stdout=None)
            if ret:
                if what in [groptics,care]:
                    mwhat, mcmd = commands_dict['makefile']
                    mwhat, mcmd = mwhat.format(software=what), mcmd.format(software=what)
                    if path_exists(mcmd.split()[1]):
                        build(mwhat, mcmd)
                        build(what, command, stdout=None)
                    else:
                        print("\n"+errorloc(mwhat,"Alternative makefile '%s' doesn't exist: compiling attempts end here\n"%mcmd.split()[1]))
                elif what == 'VBF':
                    print(error('\n   *** FORCED EXIT ***'))
                    sys.exit(1)
    print(state("\n   == Instalation finished =="))
    print(title("\n ==> Have a nice day! <==\n"))

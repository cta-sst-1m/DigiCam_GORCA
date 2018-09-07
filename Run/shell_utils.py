import shutil as sh
import os, sys, glob, docopt
import gzip as gz
import tarfile
import subprocess as sp
#import pickle as pk

def execute(command, verbose=0, dryrun=None, **Popen_kwd):
    stdout, stderr = Popen_kwd.pop("stdout",sp.PIPE), Popen_kwd.pop("stderr",sp.STDOUT)
    shell, executable = Popen_kwd.pop("shell",True), Popen_kwd.pop("executable","/bin/bash")
    if verbose:
        print(" [%s] command line: %s" % (("execute" if dryrun is None else "dryrun"),command))
        if verbose>1:
            print(", ".join([" [subprocess] Popen kwd: stdout=%s, stderr=%s"]+(["executable=%s"%executable] if shell else [])+["%s=%s"%(k,Popen_kwd[k]) for k in Popen_kwd]))
    if dryrun is None:
        proc = sp.Popen(command, shell=shell, executable=executable, stdout=stdout, stderr=stderr, **Popen_kwd)
        if stdout==sp.PIPE:
            output, _ = proc.communicate()
        else:
            proc.wait()
        retcode = proc.poll()
        if retcode:
            if stdout==sp.PIPE:
                print(output)
            raise sp.CalledProcessError(retcode, command, output=output if stdout==sp.PIPE else None)
    else:
        output = dryrun+"\n"
    return output if stdout==sp.PIPE else None

def abspath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def mkdirs(*dirs):
    for d in dirs:
        d = abspath(d)
        if not d==abspath('') and not path_exists(d):
            os.makedirs(d)
        else:
            print(" ('%s' NOT created: the path already exists or is the current one)"%d)

def path_exists(path):
    return os.path.exists(abspath(path))#access(path, os.F_OK)

def process_success(process_output):
    return process_output.returncode==0

def pwd():
    return abspath(os.environ["PWD"])

def symb_link(what, to_what, verbose=False):
    out = execute("ln -s %s %s"%(what, to_what), verbose=verbose)
    if verbose:
        print(out)

def gzip_file(file_in,file_out=None):
    with gz.open((file_in+".gz") if file_out is None else file_out, "wb") as g:
        g.write(file(file_in).read())

def tar_folder(folder_in,folder_out=None):
    with tar.open((folder_in+".tar.gz") if folder_out is None else folder_out, "w:gz" ) as t:
        for name in os.listdir(folder_in):
            t.add(name)

def qsub(executable, cwd=None, name=None, log=None, queue=None, verbose=0, PBS=True, dryrun_id=None, **extra):
    popen_kwd = dict(verbose=verbose, dryrun=dryrun_id, env=extra.pop('env',os.environ))
    if cwd is not None:
        executable = "cd %s; %s" % (cwd, executable)
        popen_kwd['cwd'] = cwd
    cmd = 'echo "%s"' % executable
    qsub_args = ['-V','-j oe' if PBS else '-j y']
    if name is not None:
        qsub_args += ['-N %s' % name]
    if log is not None:
        if dryrun_id is None:
            mkdirs(os.path.dirname(log))
        qsub_args += ['-o %s' % log]
    if queue is not None:
        qsub_args += ['-q %s' % queue]
    qsub_args = ' '.join(qsub_args+['-l %s=%s'%(k,extra[k]) for k in extra])
    commandline = cmd+" | qsub "+qsub_args
    stdout = execute(commandline, **popen_kwd)
    print(" [qsub] job '%s' (batchlog: %s)" % (name, log))
    print("        submitted with job_id %s" % stdout)

def bash_source(source_file, skip_keys=None, verbose=False):
    proc = sp.Popen('source %s && env'%source_file, shell=True, executable="/bin/bash", stdout=sp.PIPE, stderr=sp.STDOUT)
    if verbose:
        print("...Current environments:")
    for line in proc.stdout:
        (key, _, value) = str(line).partition("=")
        if skip_keys is not None and key in skip_keys:
            if verbose:
                print(" [environment] '%s' is skipped..." % key)
            continue
        os.environ[key.strip()] = value.strip()
        if verbose:
            print(" [environment] '%s': %s" % (key.strip(), value.strip()))
    if verbose:
        print()
    proc.communicate()


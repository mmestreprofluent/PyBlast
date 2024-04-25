import os, tempfile, time, sys
from multiprocessing import Pool, cpu_count
from subprocess import check_output

from Bio import SeqIO

class TMPFname() :

    def __init__(self, delete=True, ext="", quiet=False, ** kwargs) :
        suffix = self.format_ext(ext)
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, ** kwargs)
        self.fname = tmp.name 
        self.delete = delete
        self.quiet = quiet

    def __str__(self) :
        return self.fname

    def format_ext(self, ext) :
        ext = str(ext)
        if not ext : return ext
        if not ext.startswith(".") : ext = "." + ext.strip()
        return ext

    def exist(self) :
        return os.path.isfile(self.fname)

    def remove(self) :
        if self.exist() and self.delete :
            if not self.quiet : print ("Delete temporary file at : %s" %(self.fname))
            os.remove(self.fname)

    def __del__(self) :
        self.remove()

class TMPFasta(TMPFname) :

    def __init__(self, data, delete=True, wodesc=False, ** kwargs) :
        super(TMPFasta, self).__init__(delete=delete, ext="fa", ** kwargs)
        if data : self.write(data, wodesc=wodesc)

    @staticmethod
    def remove_desc(feature) :
        feature.description = ""
        return feature

    def write(self, data, mode="w", wodesc=False) :
        if wodesc : data = (TMPFasta.remove_desc(feat) for feat in data)
        with open(self.fname, mode) as handle :
            SeqIO.write(data, handle, "fasta")


"""
PROCESSING
"""

# decorator to time functions
def ftiming(f):
    def wrap(*args, **kwargs):
        # Receive 'quiet' from kwargs and use it directly
        quiet = kwargs.get('quiet', False)  # This will not modify kwargs

        funcname = kwargs.pop("funcname", f.__name__)
        funcname = funcname or "Unknown funcname"

        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()

        # Only print if 'quiet' is False
        if not quiet:
            print('{:s} function took {:.3f} ms'.format(funcname, (time2 - time1) * 1000.0))

        return ret

    return wrap

def run_cmdline(cmdline, funcname=None, quiet=False):
    # Pass all original parameters to the function
    if not quiet:
        return run_cmdline_loud(cmdline, funcname=funcname, quiet=quiet)
    return check_output(cmdline, universal_newlines=True, shell=True)

@ftiming
def run_cmdline_loud(cmdline, funcname=None, quiet=False):
    return check_output(cmdline, universal_newlines=True, shell=True)

"""
MULTIPROCESSING
"""

class FunArgs() :
 
    def __init__(self, fun, * args, ** kwargs) :
        self.fun = fun
        self.args = args
        self.kwargs = kwargs
 
    def launch_fun(self) :
        return self.fun(* self.args, ** self.kwargs)
 
class Multiprocess() :
 
    @staticmethod
    def lambda_fun(farg) :
        return farg.launch_fun()
 
    def run(self, fargs, ncore=1, quiet=False) :
        if not quiet: print ("ncore : %i - available : %i" %(ncore, cpu_count()))
        if ncore == 1 : return [farg.launch_fun() for farg in fargs]
 
        pool = Pool(ncore)
        func = Multiprocess.lambda_fun
        fargs = ((farg, ) for farg in fargs)

        try :
            data = pool.starmap(func, fargs)
            pool.close()
            pool.join()
            return data
 
        except KeyboardInterrupt :
            print("Interrupt childs process")
            pool.terminate()
            sys.exit()
 
def mproc(fargs, ncore=1, quiet=False) :
    mp = Multiprocess()
    return mp.run(fargs, ncore=ncore, quiet=quiet)
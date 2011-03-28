
import math, sys, os
from lxml import etree

def which(program):
    """
    Mimics the behavior of the UNIX 'which' command.
    From: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


# since we cannot rely on numpy being available, we
# define a floating point range function
# From:
# http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments/
#  
def frange6(*args):
    """
    Range function that accepts floats (and integers).
    
    Usage:
    frange(-2, 2, 0.1)
    frange(10)
    frange(10, increment = 0.5)
    
    The returned value is an iterator.  Use list(frange) for a list.
    """
    start = 0.0
    step = 1.0

    l = len(args)
    if l == 1:
        end = args[0]
    elif l == 2:
        start, end = args
    elif l == 3:
        start, end, step = args
        if step == 0.0:
            raise ValueError, "step must not be zero"
    else:
        raise TypeError, "frange expects 1-3 arguments, got %d" % l

    v = start
    while True:
        if (step > 0 and v >= end) or (step < 0 and v <= end):
            raise StopIteration
        yield v
        v += step


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

        sys.stdout.write("\nSaving band structure to %s... " % outFile)
        sys.stdout.flush()
        plt.savefig(outFile)
        sys.stdout.write("done!\n\n")

def OpenVaspRunXml():
    filename = 'vasprun.xml'
    if len(sys.argv) > 0:
        filename = sys.argv[0]
    filename = 'vasprun.xml'
    if not os.path.isfile(filename):
        print "\nERROR: File '%s' was not found!\n" % (filename)
    fsize = float(os.path.getsize(filename))/1024**2
    sys.stdout.write("Parsing %s (%.1f MB) using libxml... " % (filename, fsize))
    if fsize > 100:
        sys.stdout.write("\nThis may take some time... ")
    sys.stdout.flush()
    doc = etree.fromstring(open(filename,'r').read())
    sys.stdout.write("done!\n")
    return doc, filename


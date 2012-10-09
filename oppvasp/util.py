
import sys, os
import numpy as np

try:
    from lxml import etree
except:
    print "Warning: lxml python package not found"

def get_pairs(n, indices = None, include_self = False):
    """
        Returns a list of all n*(n-1)/2 unordered pairs of numbers 0...n-1
        as a (n*(n-1)/2, 2) numpy array.
            include_self : bool (default False)
                Whether to check for bonds between an atom and its _own_ periodic image.
                This is only needed for very small unit cells, typically with a single atom.

        Example: n = 5. 
            # pairs: n(n-1)/2 = 10

            j | n-1-j | jsum  | indices | i,j pairs
            --|-------|-------|---------|----------------
            0 |   4   |  0    | 0:4     | 1,0 2,0 3,0 4,0
            1 |   3   |  3    | 4:7     | 2,1 3,1 4,1
            2 |   2   |  5    | 7:9     | 3,2 4,2
            3 |   1   |  6    | 9:10    | 4,3 

        Example: n = 5 with include_self = True
            # pairs: (n+1)n/2 = 15

            j |  n-j  | jsum  | indices | i,j pairs
            --|-------|-------|---------|----------------
            0 |   5   |   0   | 0:5     | 0,0 1,0 2,0 3,0 4,0
            1 |   4   |   4   | 5:9     | 1,1 2,1 3,1 4,1
            2 |   3   |   7   | 9:12    | 2,2 3,2 4,2
            3 |   2   |   9   | 12:14   | 3,3 4,3 
            4 |   1   |  10   | 14:15   | 4,4
    """
    if indices == None:
        indices = np.arange(n)
    else:
        indices = np.array(indices)
    if include_self:
        m = n
    else:
        m = n-1
    pairs = np.empty((m*(m+1)/2,2), dtype = int)
    for j in xrange(m):
        if include_self:
            # jsum = sum(n-j, n-1)
            jsum = j*(2*n-j-1)/2
        else:
            # jsum = sum(n-1-j, n-2)
            jsum = j*(2*n-j-3)/2
        #print jsum
        #print "j=%d" % j, " =>",j+isum,j+isum+(n-1-j)
        #print pairs[j+isum:j+isum+(n-1-j),1]
        pairs[jsum+j:jsum+m,1] = indices[j]
        if include_self:
            pairs[jsum+j:jsum+m,0] = indices[range(j,n)]
        else:
            pairs[jsum+j:jsum+m,0] = indices[range(j+1,n)]
    return pairs

def unitcell_components(cell):
    A = cell[0]
    B = cell[1]
    C = cell[2]
    a = np.sqrt(np.dot(A,A))
    b = np.sqrt(np.dot(B,B))
    c = np.sqrt(np.dot(C,C))
    alpha = np.arccos( np.dot(B,C) / (b*c) ) * 360 / (2.*np.pi)
    beta  = np.arccos( np.dot(C,A) / (c*a) ) * 360 / (2.*np.pi)
    gamma = np.arccos( np.dot(A,B) / (a*b) ) * 360 / (2.*np.pi)
    return a,b,c,alpha,beta,gamma

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

# DEPRECATED
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

class FileType(object):
    """
    Copied from argparse, but adding a warning before overwriting files!
    """

    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, string):
        # the special argument "-" means sys.std{in,out}
        if string == '-':
            if 'r' in self._mode:
                return _sys.stdin
            elif 'w' in self._mode:
                return _sys.stdout
            else:
                msg = _('argument "-" with mode %r') % self._mode
                raise ValueError(msg)

        # all other arguments are used as file names
        if self._mode == 'w' and os.path.exists(string):
            if query_yes_no("The file %s exists. Do you want to overwrite it?" % string) == "no":
                sys.exit(1)
        try:
            return open(string, self._mode, self._bufsize)
        except IOError as e:
            message = _("can't open '%s': %s")
            raise ArgumentTypeError(message % (string, e))

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '%s(%s)' % (type(self).__name__, args_str)

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


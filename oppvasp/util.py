
import sys, os
import numpy as np

try:
    from lxml import etree
except:
    print "Warning: lxml python package not found"

def direct_to_cartesian(positions, basis):
    """
    Converts positions in direct coordinates into cartesian coordinates using a given basis.

    Parameters
    ----------
    positions : num_steps x num_atoms x 3 numpy array
        Array containing the direct coordinates of num_atoms atoms for num_steps steps
    basis : num_steps x 3 x 3 numpy array
        Array containg the lattice vectors in Cartesian coordinates

    """
    if positions == None:
        return None
    #print "converting to cartesian basis..."
    t1 = time.clock()
    
    # Alternative 1
    #pos = np.array([np.dot(p,b) for p,b in zip(positions, basis)]) # there is surely some faster way!

    if basis.ndim == 3:
        # time-dependent basis
    
        # Alternative 2 is about 10 times faster than alternative 1
        pos = np.zeros(positions.shape)
        if pos.ndim == 2 or pos.ndim == 3:
            i = 0
            for p,b in zip(positions,basis):
                pos[i] = np.dot(p,b)
                i += 1
        #elif pos.ndim == 2:
            # single coordinate set
            # for some reason, this is really slow:
            #   pos = np.dot(positions, basis)

        else:
            raise StandardError("positions is of wrong dimensions")
    
    else:

        if np.allclose(basis[[0,0,1,1,2,2],[1,2,0,2,0,1]],np.zeros(6)):
            # orthorhombic
            if np.allclose(basis.diagonal(),np.tile(basis[0,0],3)):
                # cubic
                pos = positions * basis[0,0]
            else:
                print "not supported yet"
        else:
            print "not supported yet"

    
    # Alternative 3: if the cell is static (3 x 3 instead of nsteps x 3 x 3):
    # perhaps we can use tensordot anyway? I'm not sure.
    #pos = np.tensordot(positions, basis, axes=([2],[0]))
    
    tdiff = time.clock() - t1
    if tdiff > 1:
        print "Conversion to cartesian basis took",round(tdiff, 3),"seconds. This function should be optimized!"
    # We could make a Fortran module doing the conversion... 
    #for i in range(pos.shape[0]): # axis 0   : 4 
    #    pos3[i] = np.dot(pos[i],basis[i])
    #    for j in range(pos.shape[1]): # axis 1 : 5
    #        #pos2[i,j] = np.dot(pos[i,j],basis[i,j])
    #        for k in range(pos.shape[2]): # axis 2 : 3
    #            for l in range(pos.shape[2]): # axis 2 : 3
    #                pos2[i,j,k] += pos[i,j,l] * basis[i,l,k]
    return pos


def cartesian_to_direct(positions, basis):
    """
    Converts positions in cartesian coordinates into direct coordinates using a given basis.

    Parameters
    ----------
    positions : num_steps x num_atoms x 3 numpy array
        Array containing the cartesian coordinates of num_atoms atoms for num_steps steps
    basis : num_steps x 3 x 3 numpy array
        Array containg the lattice vectors in Cartesian coordinates

    """
    if positions == None:
        return None
    pos = np.zeros(positions.shape)
    i = 0
    inv_basis = np.linalg.inv(basis)
    
    if pos.ndim == 3:
        # trajectory
        for p,b in zip(positions, inv_basis):
            pos[i] = np.dot(p,b)
            i += 1
    elif pos.ndim == 2:
        # single coordinate set
        pos = np.dot(positions, inv_basis)
    else:
        raise StandardError("positions is of wrong dimensions")

    return pos

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


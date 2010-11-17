
import math
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


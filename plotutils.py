# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

# Tip: Use \showthe\columnwidth to get column width
import numpy as np

import os,copy,sys
from copy import copy

from oppvasp.vasp.parsers import IterativeVasprunParser, PoscarParser

# print matplotlib.__version__

def prepare_canvas(width = 350.0, height = 'auto', fontsize = 10, fontsize_small = 8, lw = 0.5): 
    """
    Prepares a figure with specified width <width> and height 
    calculated according to the golden mean ratio. 
    
    Arguments:
       - width : float specifying the width in points, or string holding a value and its unit separated by space.
                 Valid units are 'cm', 'pt', 'in'
    
    Example:
        from oppvasp import plotutils
        plotutils.prepare_canvas( width = '7.2 cm' )
        import matplotlib.pyplot as plt  # import after prepare_canvas

        fig = plt.figure()
        ax1 = fig.add_axes([ 0.15, 0.15, 0.8, 0.8 ])
        ... and so on
        plt.savefig('out.pdf')
        
    """
    from matplotlib import rc
    
    inches_per_pt = 1.0/72.27 # According to TeX 
    inches_per_cm = 1.0/2.54  # 
    
    try:
        width = width.split()
        width_value = float(width[0])
        width_units = width[1]
    except AttributeError:
        width_value = float(width)
        width_units = 'pt'  # assume points
    
    if width_units == 'pt':
        fig_width = width_value * inches_per_pt
    elif width_units == 'cm':
        fig_width = width_value * inches_per_cm
    elif width_units == 'in':
        fig_width = width_value

    golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
    if height == 'auto':
        fig_height = fig_width*golden_mean      # height in inches
    else:
        try:
            height = height.split()
            height_value = float(height[0])
            height_units = height[1]
        except AttributeError:
            height_value = float(height)
            height_units = 'pt'  # assume points
        
        if height_units == 'pt':
            fig_height = height_value * inches_per_pt
        elif height_units == 'cm':
            fig_height = height_value * inches_per_cm
        elif height_units == 'in':
            fig_height = height_value
    
    fig_size = [fig_width,fig_height]

    rc('figure', figsize=fig_size)
    rc('figure.subplot', **{
        'left'   : 0.15,
        'right'  : 0.95,
        'bottom' : 0.15,
        'top'    : 0.93
        })
    rc('lines', linewidth=lw)
    #rc('font', family='serif', serif=['Latin Modern Roman','Palatino'], size=fontsize)
    rc('font', size=fontsize, family='serif')
    rc('text', usetex=False)
    rc('legend', fontsize=fontsize)
    rc('axes', labelsize=fontsize)
    rc('xtick', labelsize=fontsize_small)
    rc('ytick', labelsize=fontsize_small)

    #mpl.rcParams['figure.subplot.wspace']  = 0.2    # the amount of width reserved for blank space between subplots
    #mpl.rcParams['figure.subplot.hspace']  = 0.2    # the amount of height reserved for white space between subplots
    

def get_minmax(sets):
    mi = -1
    ma = -1
    for s in sets:
        mi0= np.min(s)
        ma0 = np.max(s)
        if mi == -1 or mi0 < mi:
            mi = mi0
        if ma == -1 or ma0 > ma:
            ma = ma0
    return mi,ma

def symmetric_running_mean(data, n, edgehandling = 'mirror'):
    """
    Calculates a symmetric running mean over the first axis (typically the time axis) 
    in a dataset 'data', in order to smooth out short-term fluctuations and highlight 
    longer-term trends or cycles.

    For each point <i> at this axis, the average is computed from the subset [i-n, i+n].
    The subset size is therefore 2*n+1.
    
    No true running mean therefore exists for the first n steps and the last n steps. 
    In these ranges

    edgehandling : either 'symmetric' or 'asymmetric' or 'mirror'
        - symmetric: smoothing window decreases gradually to zero. Makes the initial and 
          final positions correct, but unsmoothed.
        - asymmetric: smoothing windows decrases only on one side, until reaching n/2. 
          Better smoothing, but initial and final positions will be off.
        - mirroring about the edge. The current implementation mirrors about x = 0, y = 0 and x = x[-1]
          This preserves the initial position, but not the final position. 
          A line to mirror about y = y[-1] has been commented out, since it made the smoothing very
          sensitive to what happened to be the final value

    >>> data = np.ones((100,3)) * 5 * np.random.rand(100,3)
    >>> data[:,0] += np.arange(100)*.5
    >>> data[:,1] += np.arange(100)*.2
    >>> smoothed = symmetric_running_mean(data,10) 
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(data)
    >>> plt.plot(smoothed)

    """
    running_avg = np.zeros(data.shape)
    current_avg = running_avg[0]
    for i in np.arange(data.shape[0]): 
        if i <= n: # No true symmetric average exists for the first n values..
            # increase sample size gradually until reaching the desired size
            if edgehandling == 'symmetric':
                # symmetric mean makes the initial position correct, but results in large initial fluctuations
                current_avg = np.mean(data[0:1+2*i],axis=0) 
            elif edgehandling == 'asymmetric':
                # asymmetric mean is more stable, but the initial position will be off
                current_avg = np.mean(data[0:i+n], axis=0)
            elif edgehandling == 'mirror':
                # mirrors the first i+n frames about x = 0, giving correct initial positions and smoothing
                q = np.concatenate([data[0]-(data[n-i:0:-1]-data[0]), data[0:i+n+1]])
                #print i, "sample size:",len(q)
                current_avg = np.mean(q, axis=0)
        elif i >= data.shape[0]-n: # no symmetric average exists for the last n values
            # decrease sample size gradually
            r = 2*(data.shape[0]-i)
            if edgehandling == 'symmetric':
                # symmetric mean, decrease size of subset gradually until reaching 1
                current_avg = np.mean(data[-r:],axis=0) 
            elif edgehandling == 'asymmetric':
                # asymmetric mean, decrease size of subset gradually from 2n+1 to n+1
                current_avg = np.mean(data[i-n:], axis=0) 
            elif edgehandling == 'mirror':
                # mirror about x = mean(data[-n:] to keep subsetsize at 2n+1
                L = data.shape[0]
                # edgehandling = 'extrapolate'?:
                #q = np.concatenate([data[i-n:], 2*data[-n:].mean()-data[-2:L-i-3-n:-1]])
                q = np.concatenate([data[i-n:], data[-2:L-i-3-n:-1]])
                #print i, "sample size:",len(q)
                current_avg = np.mean(q, axis = 0)

        else:
            current_avg = current_avg - data[i-1-n]/(2*n+1) + data[i+n]/(2*n+1) # update symmetric running mean
        running_avg[i] = current_avg
    return running_avg

def symmetric_running_median(data, n):
    """
    Calculates a symmetric running median over the first axis in a dataset 'data'. 
    his function is surprisingly enough not too slow, but it could clearly be optimized!

    From a statistical point of view, the moving average, when used to estimate the underlying trend in a time series,
    is susceptible to rare events such as rapid shocks or other anomalies. 
    A more robust estimate of the trend is the simple moving median.
    
    """
    running_avg = np.zeros(data.shape)
    current_avg = running_avg[0]
    for i in np.arange(data.shape[0]): 
        if i <= n: # No true average exists for the first n values..
            current_avg = np.median(data[0:1+2*i], axis=0) # assymetric median
        elif i >= data.shape[0]-n: # no symmetric median exists for the last n values
            r = 2*(data.shape[0]-i)
            current_avg = np.median(data[-r:], axis=0) # assymetric median
        else:
            current_avg = np.median(data[i-1-n:i+n], axis=0) # update symmetric running median 
        running_avg[i] = current_avg
    return running_avg


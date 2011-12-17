# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

# Tip: Use \showthe\columnwidth to get column width
import numpy as np

import os,copy,sys
from copy import copy
from scipy.optimize import leastsq

from oppvasp.vasp.parsers import IterativeVasprunParser, PoscarParser
from oppvasp.md import Trajectory 

from matplotlib import rc
import matplotlib.pyplot as plt
# print matplotlib.__version__

def prepare_canvas(width = 350.0, fontsize = 10, fontsize_small = 8, lw = 0.5): 
    """
    Prepares a figure with specified width <width> and height 
    calculated according to the golden mean ratio. 
    
    Arguments:
       - width : float specifying the width in points, or string holding a value and its unit separated by space.
                 Valid units are 'cm', 'pt', 'in'
    
    Example:
        prepare_canvas( width = '7.2 cm' )
    """
    try:
        width = width.split()
        width_value = float(width[0])
        width_units = width[1]
    except AttributeError:
        width_value = float(width)
        width_units = 'pt'  # assume points
    
    inches_per_pt = 1.0/72.27 # According to TeX 
    inches_per_cm = 1.0/2.54  # 
    if width_units == 'pt':
        fig_width = width_value * inches_per_pt
    elif width_units == 'cm':
        fig_width = width_value * inches_per_cm
    elif width_units == 'in':
        fig_width = width_value

    golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
    fig_height = fig_width*golden_mean      # height in inches
    fig_size = [fig_width,fig_height]

    rc('figure', figsize=fig_size)
    rc('figure.subplot', **{
        'left'   : 0.15,
        'right'  : 0.95,
        'bottom' : 0.15,
        'top'    : 0.93
        })
    rc('lines', linewidth=lw)
    rc('font', family='sans-serif', serif=['Latin Modern Roman','Palatino'], size=fontsize)
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

def symmetric_running_mean(data, n):
    """
    Calculates a symmetric running mean over the first axis (typically the time axis) 
    in a dataset 'data', in order to smooth out short-term fluctuations and highlight 
    longer-term trends or cycles.

    For each point <i> at this axis, the average is computed from the subset [i-n, i+n].
    The subset size is therefore 2*n+1.
    
    No true running mean therefore exists for the first n steps and the last n steps. 
    In these ranges

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
            current_avg = np.mean(data[0:1+2*i],axis=0) # assymetric mean, increase size of subset gradually until reaching the desired size
        elif i >= data.shape[0]-n: # no symmetric average exists for the last n values
            r = 2*(data.shape[0]-i)
            current_avg = np.mean(data[-r:],axis=0) # assymetric mean, decrease size of subset gradually until reaching zero 
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



class DisplacementPlot(object):
    """
    The DisplacementPlot class is a simple script for making 
    some quick plots from MD trajectories.
    """

    def __init__(self, trajectories):
        """
        Construct a DisplacementPlot object.

        Parameters
        ----------
        trajectories : Trajectory or array
            Either a single Trajectory object or a list of such objects.
        
        Examples
        ----------
        >>> from oppvasp.vasp.parsers import read_trajectory
        >>> traj = read_trajectory(unwrap_pbcs = True)
        >>> # To only plot the part from step 10000 to step 20000:
        >>> traj.set_selection(10000, 20000)
        >>> 
        >>> dp = DisplacementPlot(traj)
        >>> dp.add_plot( what = 'r2', smoothen = True, 
        >>>     style = { 'color' : 'black' } ) # avg r^2 for all atoms 
        >>> dp.add_plot( what = 'r2', atom_no = 0, smoothen = True, linear_fit = True)
        >>> dp.add_plot( what = 'r2', atom_no = 0, smoothen = False, 
        >>>     style = { 'zorder' : -1, 'alpha': 0.4, 'color': 'gray' } )

        """
        if not 'pop' in dir(trajectories):
            trajectories = [trajectories]
        for traj in trajectories:
            traj.r2 = traj.get_displacements_squared(coordinates = 'cartesian')  # displacement r^2
            traj.avg_r2 = np.sum(traj.r2, axis=1) / traj.num_atoms       # displacement r^2 averaged over all atoms
        self.trajs = trajectories
        
        prepare_canvas('10 cm') 
        p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        self.fig = plt.figure()
        self.ax1 = self.fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        self.ax1.grid(True, which='major', color = 'gray', alpha = 0.5)
        self.ax1.set_xlabel(r'Time [ps]')
        self.plotdata = []
    
    def get_diffusion_coeff(self, atoms):
        """ 
        Returns the diffusion coefficient for a set of atoms in cm^2/s. 

        Parameters
        ----------
        atoms : list
            The atoms to calculate the diffusion coefficient for. 
            The atoms should be of the same type.
        """
        for traj in self.trajs:
            # to convert from Å^2/fs to cm^2/s: Å^2/fs * (1e-8 cm/Å)^2 * 1e15 s/fs = 0.1 cm^2/s
            # divide by 6 in 3 dimensions:
            fac = 0.1 * 1./6
            D = np.zeros(traj.r2.shape[0])
            # note: at first it seems more intelligible to just use traj.r2[1:atoms] since 
            # NumPy allows such a splicing.
            # but for some reason this completely hangs the computer. Perhaps a bug in NumPy.
            for atno in atoms:
                D[1:] += traj.r2[1:,atno] / traj.time[1:] 
            traj.D = D * fac / len(atoms)
        return [traj.D for traj in self.trajs]


    def add_plot(self, traj_no = 0, atom_no = -1, what = 'r2', smoothen = False, style = { 'color': 'black' }, 
            linear_fit = False):
        """
        Adds a plot for the property 'what' for atom 'atom_no' or an average for all the atoms if 'atom_no = -1.
        'what' can take the following values:
             'r' : displacement plot
            'r2' : square displacement plot 
        """
        if not 'ax1' in dir(self):
            self.prepare_plot()
        traj = self.trajs[traj_no]
        x = traj.time[:]/1.e3
        if what == 'r2' or what == 'r':
            self.ax1.set_ylabel(ur'$\Delta r^2$ [Å$^2$]')
            if atom_no == -1:
                y = traj.avg_r2[:]
            else:
                y = traj.r2[:, atom_no]
            if what == 'r':
                y = np.sqrt(y)
                self.ax1.set_ylabel(ur'$\Delta r$ [Å]')
        elif what == 'x' or what == 'y' or what == 'z':
            r = traj.get_displacements( coordinates = 'cartesian' )  # displacement vectors
            if atom_no == -1:
                raise ValueError("averaging not implemented")
            else:
                if what == 'x':
                    y = r[:, atom_no, 0]
                    self.ax1.set_ylabel(ur'$\Delta x$ [Å]')
                elif what == 'y':
                    y = r[:, atom_no, 1]
                    self.ax1.set_ylabel(ur'$\Delta y$ [Å]')
                elif what == 'z':
                    y = r[:, atom_no, 2]
                    self.ax1.set_ylabel(ur'$\Delta z$ [Å]')
        elif what == 'D':
            D = self.get_diffusion_coeff([atom_no])
            y = D[0]
            self.ax1.set_ylabel(ur'$D$ [cm$^2$/s]')
        
        if not x.shape == y.shape:
            raise ValueError("Uh oh, shape of x,",x.shape,", does not match shape of y,",y.shape) 

        if smoothen:
            #y = symmetric_running_mean(y, 250)
            y = symmetric_running_mean(y, 500)
            #y = symmetric_running_median(y, 1000)
        self.ax1.plot(x,y, **style)
        self.ax1.set_xlim(x[0],x[-1])
        self.plotdata.append( { 'x': x, 'y': y } )

        if linear_fit:
            fitfunc = lambda p, x: p[0] + p[1]*x  # Target function 
            errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
            p0 = [0,np.max(y)/np.max(x)] # Initial guess for the parameters
            p1, success = leastsq(errfunc, p0[:], args=(x,y))
            print p1
            self.ax1.plot(x, fitfunc(p1, traj.time/1.e3), color = 'red' )

    def save_plot(self, filename = 'displacement.pdf'):
        sys.stdout.write("Writing %s... " % os.path.basename(filename))
        sys.stdout.flush()
        plt.savefig(filename)
        sys.stdout.write("done!\n")



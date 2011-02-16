# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

# Tip: Use \showthe\columnwidth to get column width
import numpy as np
from matplotlib import rc

import os,copy,sys
from copy import copy
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

from oppvasp import getAtomicNumberFromSymbol
from oppvasp.vasp.parsers import IterativeVasprunParser, PoscarParser
from oppvasp.md import Trajectory, pair_correlation_function


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
    rc('lines', linewidth=lw)
    rc('font', family='sans-serif', serif=['Latin Modern Roman','Palatino'], size=fontsize)
    rc('text', usetex=False)
    rc('legend', fontsize=fontsize)
    rc('axes', labelsize=fontsize)
    rc('xtick', labelsize=fontsize_small)
    rc('ytick', labelsize=fontsize_small)

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
    
    No true srunning mean therefore exists for the first n steps and the last n steps. 
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




class DisplacementPlot:

    def __init__(self, last_step = -1):
        self.trajs = []
        self.last_step = last_step

    def read_trajectory(self, dir = '', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', npz_file = 'trajectory_nopbc.npz', POSCAR_file = ''):
        """
        This function will read a trajectory from the vasprun.xml file given as xml_file and save 
        the trajectory as a NumPy npz cache file npz_pbc_file. The function will then unwrap the 
        periodic boundary conditions and save the trajectory for the unwrapped trajectory as a 
        new NumPy npz cache file npz_file.

        The function will check if a npz file exists, and if it does, read the npz cache instead of the 
        vasprun.xml file. This is *much* faster.

        If POSCAR_file is specified, the initial positions will be read from this file. This may be 
        useful for some visualization purposes. While coordinates in POSCAR (and CONTCAR) may be 
        negative, the coordinates in vasprun.xml and XDATCAR are always wrapped around to become positive.
        """
        
        if os.path.isfile(dir + npz_file):
            traj = Trajectory(filename = dir +npz_file)
        else:
            p = IterativeVasprunParser(dir + xml_file)
            traj = p.get_all_trajectories()
            traj.save(dir + npz_pbc_file)
            if POSCAR_file != '':
                poscar = PoscarParser(dir + POSCAR_file)
                pos = poscar.get_positions( coordinates = 'direct' )
                print "Unwrapping using given initial pos"
                traj.unwrap_pbc( init_pos = pos)
            else:
                traj.unwrap_pbc()
            traj.save(dir + npz_file)
        traj.r2 = traj.get_displacements_squared(coordinates = 'cartesian')  # displacement r^2
        traj.avg_r2 = np.sum(traj.r2, axis=1) / traj.num_atoms       # displacement r^2 averaged over all atoms
        self.trajs.append(traj)
        if self.last_step == -1 or traj.avg_r2.shape[0] < self.last_step:
            self.last_step = traj.avg_r2.shape[0]

    def get_diffusion_coeff(self, atom_no):
        """ Returns the diffusion coefficient for atom atom_no in cm^2/s. Assumes there is only atom of that kind. """
        for traj in self.trajs:
            # to convert from Å^2/fs to cm^2/s: Å^2/fs * (1e-8 cm/Å)^2 * 1e15 s/fs = 0.1 cm^2/s
            # divide by 6 in 3 dimensions:
            fac = 0.1 * 1./6
            traj.D = traj.r2[0:self.last_step,atom_no] / traj.time[0:self.last_step] * fac
        return [traj.D for traj in self.trajs]

    def prepare_plot(self):
        prepare_canvas('10 cm') 
        p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        self.fig = plt.figure()
        self.ax1 = self.fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        self.ax1.grid(True, which='major', color = 'gray', alpha = 0.5)
        self.ax1.set_xlabel(r'Time [ps]')
        self.ax1.set_ylabel(ur'$\Delta r^2 = [\vec{r}(t)-\vec{r}(0)]^2$ [Å$^2$]')
        self.plotdata = []

    def add_plot(self, traj_no = 0, atom_no = -1, what = 'r2', smoothen = False, style = { 'color': 'black' }, linear_fit = False, fit_startstep = 0):
        """
        Adds a plot for the property 'what' for atom 'atom_no' or an average for all the atoms if 'atom_no = -1.
        'what' can take the following values:
             'r' : displacement plot
            'r2' : square displacement plot 
        """
        if not 'ax1' in dir(self):
            self.prepare_plot()
        traj = self.trajs[traj_no]
        x = traj.time[0:self.last_step]/1.e3
        if what == 'r2' or what == 'r':
            if atom_no == -1:
                y = traj.avg_r2[0:self.last_step]
            else:
                y = traj.r2[0:self.last_step, atom_no]
            if what == 'r':
                y = np.sqrt(y)
        elif what == 'x' or what == 'y' or what == 'z':
            r = traj.get_displacements( coordinates = 'cartesian' )  # displacement vectors
            if atom_no == -1:
                raise ValueError("averaging not implemented")
            else:
                if what == 'x':
                    y = r[0:self.last_step, atom_no, 0]
                elif what == 'y':
                    y = r[0:self.last_step, atom_no, 1]
                elif what == 'z':
                    y = r[0:self.last_step, atom_no, 2]
        
        if not x.shape == y.shape:
            raise ValueError("Uh oh. x shape:",x.shape," y: ",y.shape) 

        if smoothen:
            #y = symmetric_running_mean(y, 250)
            y = symmetric_running_mean(y, 500)
            #y = symmetric_running_median(y, 1000)
        self.ax1.plot(x,y, **style)
        self.ax1.set_xlim(0,x[-1])
        self.plotdata.append( { 'x': x, 'y': y } )

        if linear_fit:
            fitfunc = lambda p, x: p[0] + p[1]*x  # Target function
            errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
            p0 = [0,np.max(y)/np.max(x)] # Initial guess for the parameters
            p1, success = leastsq(errfunc, p0[:], args=(x[fit_startstep:],y[fit_startstep:]))
            self.ax1.plot(x[fit_startstep:], fitfunc(p1, time[fit_startstep:]), **styles[0])

    def save_plot(self, filename = 'displacement.pdf'):
        sys.stdout.write("Writing %s... " % os.path.basename(filename))
        sys.stdout.flush()
        plt.savefig(filename)
        sys.stdout.write("done!\n")



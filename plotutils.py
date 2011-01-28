# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

# Tip: Use \showthe\columnwidth to get column width
import numpy as np
from matplotlib import rc

import os,copy,sys
from copy import copy
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

from oppvasp import getAtomicNumberFromSymbol
from oppvasp.vasp.parsers import IterativeVasprunParser 
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
    Each average is calculated as the average during an interval of <n> steps before and <n> steps 
    after the current point. The first and last <n> steps are therefore "undefined".
    """
    running_avg = np.zeros(data.shape[0])
    current_avg = 0.
    for i in np.arange(data.shape[0]): 
        if i <= n: # No true mean exists for the first n values..
            current_avg = np.mean(data[0:i+n+1]) # assymetric mean
        elif i >= data.shape[0]-n: # no symmetric mean exists for the last n values
            current_avg = np.mean(data[i-n:-1]) # assymetric mean
        else:
            current_avg = current_avg - data[i-1-n]/(2*n+1) + data[i+n]/(2*n+1) # update symmetric running mean
        running_avg[i] = current_avg
    return running_avg

def symmetric_running_median(data, n):
    """
    Symmetric running median. This function is surprisingly enough not too slow,
    but it is clearly in need for optimization :)
    """
    running_avg = np.zeros(data.shape[0])
    current_avg = np.median(data[0:n+1+n]) # symmetric median at n
    running_avg[0:n+1] = current_avg # No true median exists for the first n values! 
    for i in np.arange(n+1,data.shape[0]-n,1):
        current_avg = np.median(data[i-1-n:i+n]) # update symmetric running median 
        running_avg[i] = current_avg
    running_avg[-n:] = current_avg # No true median exists for the last n values! 
    return running_avg




class DisplacementPlot:

    def __init__(self):
        self.trajs = []

    def add_trajectory(self, xml_filename ='vasprun.xml', npz_pbc_filename = 'trajectory_pbc.npz', npz_filename = 'trajectory_nopbc.npz'):
        """
        This function will read a trajectory from the vasprun.xml file given as xml_filename and save the trajectory 
        as a NumPy npz cache file npz_pbc_filename. The function will then unwrap the periodic boundary conditions and 
        save the trajectory for the unwrapped trajectory as a new NumPy npz cache file npz_filename.

        The function will check if a npz file exists, and if it does, read the npz cache instead of the vasprun.xml file. 
        This is *much* faster.
        """
        
        if os.path.isfile(npz_filename):
            traj = Trajectory(filename = npz_filename)
        else:
            p = IterativeVasprunParser(xml_filename)
            traj = p.get_all_trajectories()
            traj.save(npz_pbc_filename)
            traj.unwrap_pbc()
            traj.save(npz_filename)
        traj.r2 = traj.get_displacements(coordinates = 'cartesian')  # displacement r^2
        traj.avg_r2 = np.sum(traj.r2, axis=0) / traj.num_atoms       # displacement r^2 averaged over all atoms
        self.trajs.append(traj)

    def get_diffusion_coeff(self, atom_no, time_stop = -1):
        """ Calculates diffusion coefficient for atom atom_no """
        for traj in self.trajs:
            traj.D = np.sum(traj.r2[atom_no,0:time_stop])/traj.time[time_stop] * 10 # to convert from Å^2/fs to cm^2/s
        return [traj.D for traj in self.trajs]

    def prepare_plot(self):
        prepare_canvas(10*72/2.54) # 6 cm in points
        p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        self.fig = plt.figure()
        self.ax1 = self.fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        self.ax1.grid(True, which='major', color = 'gray', alpha = 0.5)
        self.ax1.set_xlabel(r'Time [ps]')
        self.ax1.set_ylabel(ur'$\Delta r^2 = [\vec{r}(t)-\vec{r}(0)]^2$ [Å$^2$]')
        self.plotdata = []
    
    def add_plot(self, traj_no = 0, atom_no = -1, smoothen = False, style = { 'color': 'black' }, linear_fit = False, fit_startstep = 0):
        """
        Adds a square displacement plot for atom atom_no,
        or an average square displacement for all the atoms if atom_no = -1.
        """
        if not 'ax1' in dir(self):
            self.prepare_plot()
        traj = self.trajs[traj_no]
        x = traj.time/1.e3
        if atom_no == -1:
            y =  traj.avg_r2
        else:
            y = traj.r2[atom_no]
        if smoothen:
            y = symmetric_running_mean(y, 250)
            y = symmetric_running_mean(y, 250)
        self.ax1.plot(x,y, **style)
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


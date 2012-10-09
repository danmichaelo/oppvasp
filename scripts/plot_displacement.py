#!/usr/bin/env python
#encoding=utf8

#import matplotlib
#from matplotlib import rc
#matplotlib.use('pdf')
from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
#prepare_canvas('8 cm') 
#rc('figure.subplot', bottom = 0.20 ) # increase bottom margin

import math
import sys, os
import numpy as np
import argparse 
from oppvasp.util import get_pairs
from oppvasp import direct_to_cartesian
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress


class DisplacementPlot(object):

    def __init__(self, trajectory, atom = 0, step = 50):
        """
        DisplacementPlot class
        
        Parameters:
            - trajectory: a oppvasp.Trajectory object
            - atom: index of the atom to follow
            - step : step size
        """
        
        atpos = trajectory.get_positions( coords = 'cart', atom = atom )
        nsteps = atpos.shape[0]
        atpos = atpos - atpos[0]

        atpos = symmetric_running_mean(atpos,300) 
        #atpos = symmetric_running_mean(atpos,500) 

        print " atpos 0: ",atpos[0]
        print " atpos -1: ",atpos[-1]

        r = np.sqrt(np.sum(atpos**2,axis=1))

        self.x = trajectory.time[0:nsteps:step]/1000.
        self.y = r[0:nsteps:step]

    
    def plot(self, axes):
        """
            axes: the matplotlib.Axes object to plot to
        """
        axes.plot(self.x, self.y, color = 'blue')
        axes.set_xlabel(u'Time [ps]')
        #axes.set_ylabel(u'$r\'(t)-r\'(0)$ [Å]')
        axes.set_xlim(self.x[0],self.x[-1])
        axes.set_ylim(0.0,np.max([1.0, 1.1*np.max(self.y)]))

    def plot_to_file(self,outfile):
        import matplotlib
        matplotlib.use('pdf')
        from oppvasp.plotutils import prepare_canvas
        prepare_canvas('8 cm')
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fig.subplots_adjust(left=0.17, bottom=0.20, right=0.80, top=0.95)
        self.plot(ax)
        plt.savefig(outfile)


if __name__ == '__main__':
    from oppvasp import read_trajectory
    traj = read_trajectory(unwrap_pbcs = True)

    atom = 0
    pl = DisplacementPlot(traj, atom = atom)
    pl.plot_to_file('%s_%s.pdf' % (os.getcwd().split('/').pop(), 'displacement' ))


##!/usr/bin/env python
## -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
## vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

#import os, sys
#import numpy as np
#from optparse import OptionParser
#from oppvasp import read_trajectory
#from oppvasp.plotutils import DisplacementPlot

################################## OPTIONS ########################################

## Parse cmd line args
#parser = OptionParser()
#pdf_file = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
#parser.add_option('--pdf', dest = 'pdf_file', default = pdf_file, help = 'Destination PDF')
#parser.add_option('--xml', dest = 'xml_file', default = 'vasprun.xml', help = 'Input vasprun.xml')
#parser.add_option('--poscar', dest = 'poscar_file', default = 'POSCAR', help = 'Input POSCAR')
#parser.add_option('-f', '--firststep', metavar='STEP', dest = 'first_step', type='int', default = 1, help = 'First ion step to include')
#parser.add_option('-l', '--laststep', metavar='STEP', dest = 'last_step', type='int', default = -1, help = 'Last ion step to include')
#parser.add_option('-a', '--atom', metavar='ATOM', dest = 'atom_no', type='int', default = 0, help = 'Index of atom to watch (first index is 0)')
#(options, args) = parser.parse_args()

#print "Range: %i:%i. Atom: %i" % (options.first_step, options.last_step, options.atom_no)

############################### MAIN SCRIPT ######################################

## Read trajectory:
#traj = read_trajectory(xml_file = options.xml_file, unwrap_pbcs = True, poscar_file = options.poscar_file)
#traj.set_selection(options.first_step, options.last_step)

## Make plot
#dp = DisplacementPlot(traj)
#dp.add_plot( what = 'r2', smoothen = True, style = { 'linestyle' : '--', 'dashes' : (3,1), 'color' : 'black' } ) # avg r^2 for all atoms
##dp.add_plot( what = 'x', atom_no = 0, smoothen = True, linear_fit = True)
#dp.add_plot( what = 'r2', atom_no = 0, smoothen = True, linear_fit = False)
##dp.add_plot( what = 'r', atom_no = 0, smoothen = True, linear_fit = True)
##dp.add_plot( what = 'x', atom_no = 0, smoothen = False, style = { 'zorder' : -1, 'alpha': 0.4, 'color': 'gray' } )

#ymax = 0
#ymin = 0
#for p in dp.plotdata:
#    ymax_tmp = np.max(p['y'])
#    if ymax_tmp > ymax:
#        ymax = ymax_tmp
#    ymin_tmp = np.min(p['y'])
#    if ymin_tmp < ymin:
#        ymin = ymin_tmp
#ymax = (ymax)*1.1
#ymin = (ymin)*1.1
#dp.ax1.set_ylim(ymin,ymax)

#dp.save_plot(options.pdf_file)


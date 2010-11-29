# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
#
# @created: Nov 14, 2010 21:06:06 
# Last modified: Nov 29, 2010 12:52:02

import os,copy,sys
import matplotlib
matplotlib.use("pdf")

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from scitools.std import seq
import glob # for finding files using wildcards
from oppvasp.utils import query_yes_no
from oppvasp.vasp.parsers import vasprunParser
import oppvasp.plotting
__docformat__ = "restructuredtext en"

class ConvergenceTestPlot:
    """
    This class parses all vasprun.xml files in the current directory 
    with the filename pattern 'vasprun*.xml'. A numpy array is created
    that holds the convergence parameter and the total energy from all the xml files.
    The array is sorted on volume and then plotted usign matplotlib.

    Example usage:

    >>> plot = ConvergenceTestPlot(parameter = 'ENCUT')
    >>> plot.setXLabel('$E$')
    >>> plot.plot('out.pdf')
    >>> plot.exportCSV('out.csv')
    
    """

    def __init__(self, parameter = 'ENCUT', width = 350., normalfontsize = 10, smallfontsize = 8, linewidth = 0.5, plotmargins = [0.125, 0.15, 0.05, 0.05]):
        
        plotting.prepareCanvas(width, normalfontsize, smallfontsize, linewidth, plotmargins)
        self.plt = plt

        self.xlabel = 'Volume'
        self.ylabel = 'Energy'

        #
        # Find and parse vasprun xml-files:
        #
        xmlFiles = glob.glob('vasprun*.xml')
        print "Found %d files matching 'vasprun*.xml'" % (len(xmlFiles))
        self.xy = np.zeros((2,len(xmlFiles)))
        for i in range(len(xmlFiles)):
            p = vasprunParser(xmlFiles[i])
            self.xy[0,i] = p.getIncarProperty(parameter)
            self.xy[1,i] = p.getTotalEnergy()

        #
        # Sort self.xy by increasing parameter value 
        #
        idx = np.argsort(self.xy[0])
        self.xy = np.array([[self.xy[j,i] for i in idx] for j in [0,1]])
        
        #
        # Print to screen
        #
        print "------------------------------------------------"
        print "File\t\t"+parameter+"\t\tEnergy"
        print "------------------------------------------------"
        for i in range(self.xy.shape[1]):
            print "%s\t%.2f\t\t%.4f" % (xmlFiles[idx[i]],self.xy[0,i],self.xy[1,i])
        print 

    def setXlabel(self,lab):
        self.xlabel = lab

    def setYlabel(self,lab):
        self.ylabel = lab

    def plot(self, filename = 'convergencetestplot.pdf'):

        # Create cubic-spline interpolation between data points
        # http://www.tau.ac.il/~kineret/amit/scipy_tutorial/               
        #
        tck = interpolate.splrep(self.xy[0],self.xy[1],s=0)   # Find the B-spline representation of 1-D curve.
        cs_x = np.arange(self.xy[0,0],self.xy[0,-1],0.001)    # Create x axis
        cs_y = interpolate.splev(cs_x,tck,der=0)    # Evaluate a B-spline or its derivatives
        # yder = interpolate.splev(xnew,tck,der=1)  # Derivative of spline  

        plt.plot(self.xy[0], self.xy[1], 'bx',zorder=2)                           # Plot raw data points
        plt.plot(cs_x,cs_y,color='gray',linestyle='solid',zorder=1)     # Plot interpolated line

        #min_arg = np.argmin(self.xy[1])
        #min_pt = self.xy[:,min_arg]
        #cs_min_arg = np.argmin(cs_y)
        #cs_min_pt = [cs_x[cs_min_arg],cs_y[cs_min_arg]]

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)

        yMin = np.min(self.xy[1]); yMax = np.max(self.xy[1]); yRange = yMax-yMin
        xMin = self.xy[0,0]; xMax = self.xy[0,-1]; xRange = xMax-xMin

        plt.xlim(xMin,xMax)
        plt.ylim(yMin-yRange*.05, yMax+yRange*0.05)

        # Write to file:
        if os.path.exists(filename):
            print "\nWARNING: The file \"%s\" already exists." % filename 
            if query_yes_no("         Do you want to overwrite it?") == 'no':
                return
                #return min_pt, cs_min_pt

        sys.stdout.write("\nSaving plot to %s... " % filename)
        sys.stdout.flush()
        plt.savefig(filename)
        sys.stdout.write("done!\n\n")

        #return min_pt, cs_min_pt

    def exportCSV(self, filename = 'convergencetestplot.csv'):
        """
        Exports the internal x,y numpy array as a .csv file, 
        that can be imported into Excel or what have you.
        """
        sys.stdout.write("\nSaving table to %s... " % filename)
        sys.stdout.flush()
        f = open(filename,'w')
        for i in range(xy.shape[1]):
            f.write('%.6f,%.6f\n' % (xy[0,i],xy[1,i]))
        f.close()
        sys.stdout.write("done!\n\n")



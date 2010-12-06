# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
#
# @created: Nov 14, 2010 21:06:06 
# Last modified: Dec 06, 2010 12:33:57

import matplotlib
matplotlib.use("pdf")
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

import os,copy,sys
from scipy import interpolate
from scitools.std import seq
from oppvasp.utils import query_yes_no
__docformat__ = "restructuredtext en"

# Tip: Use \showthe\columnwidth to get column width
def prepareCanvas(fig_width_pt = 350.0, s1 = 10, s2 = 8, lw = 0.5): 
    #
    # Configure matplotlib environment
    # http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    #
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]

    rc('figure', figsize=fig_size)
    rc('lines', linewidth=lw)
    rc('font', family='serif', serif=['Latin Modern Roman','Palatino'], size=s1)
    rc('text', usetex=True)
    rc('legend', fontsize=s1)
    rc('axes', labelsize=s1)
    rc('xtick', labelsize=s2)
    rc('ytick', labelsize=s2)


class GenericPlot:

    def __init__(self, width = 350., normalfontsize = 10, smallfontsize = 8, linewidth = 0.5, plotmargins = [0.125, 0.15, 0.05, 0.05]):
        prepareCanvas(width, normalfontsize, smallfontsize, linewidth)
        self.data = []
        p = plotmargins
        self.fig = plt.figure()
        self.ax = self.fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        self.xlabel = 'Parameter'
        self.ylabel = 'Energy'

    def addData(self,xy, style = '',title=''):
        self.data.append({
            'xy': xy,
            'style': style,
            'title': title
        })

    def plotAndSave(self, filename = 'plot.pdf'):


    def plot(self, filename = 'plot.pdf'):

        # Create cubic-spline interpolation between data points
        # http://www.tau.ac.il/~kineret/amit/scipy_tutorial/               
        #
        yMin = 1000; yMax = 0; xMin = 1000; xMax = 0;
        titles = []
        for d in self.data:
            xy = d['xy']
            tck = interpolate.splrep(xy[0],xy[1],s=0)   # Find the B-spline representation of 1-D curve.
            cs_x = np.arange(xy[0,0],xy[0,-1],0.001)    # Create x axis
            cs_y = interpolate.splev(cs_x,tck,der=0)    # Evaluate a B-spline or its derivatives
            # yder = interpolate.splev(xnew,tck,der=1)  # Derivative of spline  
            
            self.ax.plot(xy[0], xy[1], **d['style'] )                           # Plot raw data points
            titles.append(d['title'])
            #plt.plot(cs_x,cs_y,color='gray',linestyle='solid',zorder=1)     # Plot interpolated line

            yMin0 = np.min(xy[1])
            if yMin0 < yMin:
                yMin = yMin0.copy()
            yMax0 = np.max(xy[1])
            if yMax0 > yMax:
                yMax = yMax0.copy()
            xMin0 = xy[0,0]
            if xMin0 < xMin:
                xMin = xMin0.copy()
            xMax0 = xy[0,-1]
            if xMax0 > xMax:
                xMax = xMax0.copy()
        if not titles[0] == '':
            self.ax.legend(tuple(titles),'upper right',frameon=False,numpoints=1)
        
        xRange = xMax-xMin
        yRange = yMax-yMin
        
        #min_arg = np.argmin(self.xy[1])
        #min_pt = self.xy[:,min_arg]
        #cs_min_arg = np.argmin(cs_y)
        #cs_min_pt = [cs_x[cs_min_arg],cs_y[cs_min_arg]]


        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)

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
    
    def setXlabel(self,lab):
        self.xlabel = lab

    def setYlabel(self,lab):
        self.ylabel = lab



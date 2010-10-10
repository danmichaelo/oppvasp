# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
#
# Created on 14. mars 2010
# @author: danmichael


import os,copy,sys
import matplotlib
matplotlib.use("pdf")

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from scitools.std import seq

from matplotlib import rc
rc('lines', linewidth=0.5)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('font', family='serif')
rc('text', usetex=True)


#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
#from matplotlib.font_manager import FontProperties
#import matplotlib


from query_yes_no import query_yes_no

class PlotBands():
    '''
    classdocs
    '''

    def __init__(self, specialpoints={}):
        '''
        Constructor
        '''
        
        # Initialise some arrays:
        self.bandstructs = []
        self.k_axis_specialPoints = []
        self.k_axis_specialPointLabels = []
        self.specialpoints = {}

        # Do a quick validation of the special points sent to us:
        if 'points' in specialpoints:
            sp = specialpoints['points']
            for i in range(len(sp)):
                try:
                    if len(sp[i]) != 2:
                        raise StandardError('k point %i: not on the format (label,vector)' % i)
                    if not isinstance(sp[i][0],str):
                        raise StandardError('k point %i: invalid label' % i)
                    if len(sp[i][1]) != 3:
                        raise StandardError('Special point %i is not a vector of length 3.' % i)
                except TypeError:
                    print '\n ERROR: Special points should be on the format: (LABEL, (x,y,z))\n' \
                      ' Point %i is not on this format:' % i 
                    print " ",sp[i],"\n"
                    raise
            self.specialpoints = specialpoints

        # Now, let's prepare arrays that can easily be plotted.
        # We will save the points for the abscissa axis in 'k_axis', while we
        # save the points for the ordinate axis in 'bands' (this will be an
        # array of arrays, since we can have several eigenvalues for each
        # k-point. We will save the k-vectors themselves in 'k_vectors' so we
        # can check for special points later on. And finally we will save the
        # unique step lengths dk between each k point to check for variation
        # in k-point density.
    
    def addBands(self, bandstructure, fermiLevel = 0.0, color='blue', label='', linestyle='solid'):
        if not bandstructure.isValid():
            return
        self.bandstructs.append({
            'bandstructure': bandstructure,
            'fermiLevel': fermiLevel,
            'color': color,
            'label': label,
            'linestyle': linestyle
        })
        
    def prepareBandstructure(self, bandstruct, basis):
        """
        This procedure will prepare a simple, scalar k-axis from the 
        vector k-s found in 'bandstruct'.
        """
        kpoints = bandstruct['bandstructure'].getBands(basis)
        kpointsCount = len(kpoints)
        bandsCount = len(kpoints[0].eigenvals)
        print "    Adding to plot: %i bands, %i kpoints" % (bandsCount, kpointsCount)
        
        k_vectors = []
        k_axis = []        
        bands = [[] for i in range(bandsCount)]
        dk_array = []
        k = 0.0
        
        if 'basis' in self.specialpoints: # If special points have been assigned
            if self.specialpoints['basis'] == 'reciprocal' and basis == 'cartesian':
                # Convert special points from reciprocal to cartesian basis:
                specialPoints = { 'basis': 'cartesian', 'points': [] }
                S = bandstruct['bandstructure'].getReciprocalLatticeVectors()
                for sp in self.specialpoints['points']:
                    specialPoints['points'].append((sp[0],np.dot(S, sp[1])))
            elif self.specialpoints['basis'] == 'cartesian' and basis == 'reciprocal':
                # Convert special points from cartesian to reciprocal basis:
                specialPoints = { 'basis': 'reciprocal', 'points': [] }
                S = bandstruct['bandstructure'].getReciprocalLatticeVectors()
                invS = np.linalg.inv(S)
                for sp in self.specialpoints['points']:
                    specialPoints['points'].append((sp[0],np.dot(invS, sp[1])))       
            else:
                # No conversion needed: 
                specialPoints = { 'basis': basis, 'points': [] }
                for sp in self.specialpoints['points']:
                    specialPoints['points'].append((sp[0],np.array(sp[1])))   
        
        for kpoint in kpoints:
            kvec = kpoint.getVector()
            k_vectors.append(kvec)
            if len(k_vectors) > 1:      # If not first k-vector:
                dkvec = kvec - k_vectors[len(k_vectors)-2]
                dk = np.sqrt(np.dot(dkvec,dkvec))
                if dk not in dk_array:
                    dk_array = np.append(dk_array,dk)
                k += dk
            k_axis.append(k)
            bandNo = 0
            for eigval in kpoint.eigenvals:
                bands[bandNo].append(eigval['energy'] - bandstruct['fermiLevel'])
                bandNo += 1
            
            # Check if we are at a special point:
            if 'points' in self.specialpoints: # If special points have been assigned
                for p in specialPoints['points']:
                    # WARNING:
                    # A difficulty arises in comparing k points, since 
                    # e.g. [0,-1,0] should be considered equivalent to [0,0,1].
                    # For simplicity, we just compare the length of the vectors.
                    # This seems to work fine, but it MAY give erroneous results!
                    if np.sqrt(np.dot(kvec,kvec)) == np.sqrt(np.dot(p[1],p[1])):
                        newPoint = True
                        for op in self.k_axis_specialPoints:
                            if abs(op-k) < 0.0001:
                                newPoint = False
                        if newPoint:
                            print "     - Found ",p[0],"point (", kvec,") at",k
                            self.k_axis_specialPoints.append(k)
                            self.k_axis_specialPointLabels.append(p[0])

        
        bandstruct['k_axis'] = k_axis
        bandstruct['E_axis'] = bands
    
    def plotAndSaveToFile(self, outFile, yRange = [-4, 4], yTicks = [-4,4,1], basis='reciprocal'): 

        print "[+] Starting plotting procedure"
        
        # Plot the bandstructures themselves:
        for bandstruct in self.bandstructs:
            
            self.prepareBandstructure(bandstruct, basis)
            
            kaxis = bandstruct['k_axis']
            xMax = kaxis[-1]              # temporary solution!!
            for i in range(len(bandstruct['E_axis'])):
                band = bandstruct['E_axis'][i]
                if len(kaxis) != len(np.unique(kaxis)):
                    raise Exception, "Uh oh, non-unique points on k-axis!"
                
                # Cubic-spline interpolation. We set s=0 for no smoothing.
                # This ensures that the line passes through all points. 
                # See http://www.tau.ac.il/~kineret/amit/scipy_tutorial/               
                tck = interpolate.splrep(kaxis,band,s=0)                
                xnew = np.arange(0,kaxis[-1],0.001)
                ynew = interpolate.splev(xnew,tck,der=0)
                if (bandstruct['label'] != "" and i == 0):
                    plt.plot(xnew,ynew,color=bandstruct['color'],linestyle=bandstruct['linestyle'],label=bandstruct['label'], zorder=2)
                else:
                    plt.plot(xnew,ynew,color=bandstruct['color'],linestyle=bandstruct['linestyle'], zorder=2)
                # Alternative: Plot data directly:
                # plt.plot(kaxis,band,'.')
                
                # Alternative 2:
                #xnew = np.arange(0,kaxis[-1],0.001)
                #intf = interp1d(xnew,band,2)
                #plt.plot(kaxis,intf(xnew),'b')

        # Plot horizontal and vertical lines behind the bandstructures:
        xticksX = []
        xticksLabels = []
        plt.axhline(y=0, color='grey', zorder=1)
        for i in range(len(self.k_axis_specialPoints)):
            #print " - Adding special point ",self.k_axis_specialPointLabels[i]
            plt.axvline(x=self.k_axis_specialPoints[i], color='black', zorder=1)
            xticksX.append(self.k_axis_specialPoints[i])
            xticksLabels.append(self.k_axis_specialPointLabels[i])
        # Set ranges, ticks and legend:
        plt.xlim(0, xMax)
        plt.ylim(yRange[0], yRange[1])
        plt.yticks(seq(yTicks[0],yTicks[1],yTicks[2]))
        if len(xticksX) != 0:
            plt.xticks(xticksX, xticksLabels)
            plt.legend()
        
        # Write to file:
        if os.path.exists(outFile):
            print "\nWARNING: The file \"%s\" already exists." % outFile
            if query_yes_no("         Do you want to overwrite it?") == 'no':
                return

        sys.stdout.write("\nSaving band structure to %s... " % outFile)
        sys.stdout.flush()
        plt.savefig(outFile)
        sys.stdout.write("done!\n\n")


if __name__ == '__main__':
    import abpytho.vasp.parsebandstructure as vpb
    import abpytho.espresso.parsebandstructure as epb
    
    specialKPoints = { 
        'basis': 'reciprocal', 
        'points': (
            ('L'      , (0.5, 0.5,  0.5)),
            ('\\Gamma', (0.0, 0.0, 0.0)),
            ('X'      , (0.5, 0.5, 0.0)),
            ('W'      , (0.5, 0.75, 0.25))
        )
    }    
    plot = PlotBands(specialKPoints)
    
    plot.addBands(vpb.ParseBandStructure('/Users/danmichael/Dropbox/code/abpytho/src/tests/PROCAR','/Users/danmichael/Dropbox/code/abpytho/src/tests/OUTCAR'), 
        fermiLevel = 6.54, color='blue', label = 'VASP', linestyle='dashed')
    plot.addBands(epb.ParseBandStructure('/Users/danmichael/Dropbox/code/abpytho/src/tests/si.band','/Users/danmichael/Dropbox/code/abpytho/src/tests/si.band.out'), 
        fermiLevel = 6.54, color='red', label = 'Quantum Espresso')
    plot.plotAndSaveToFile(outFile = 'bands.pdf', yMin = -13, yMax = 8.0, basis='cartesian') 

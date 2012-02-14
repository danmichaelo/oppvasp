# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
#
# Created on 14. mars 2010
# @author: danmichael


import os,copy,sys

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from scitools.std import seq

#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
#from matplotlib.font_manager import FontProperties
#import matplotlib


from util import query_yes_no

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
    
    def addBands(self, bandstructure, fermiLevel = 0.0, label='', style = {}):
        #color='blue', label='', linestyle='solid', marker='d', markersize=2.0):
        #if not bandstructure.isValid():
        #    return
        self.bandstructs.append({
            'bandstructure': bandstructure,
            'fermiLevel': fermiLevel,
            'label': label,
            'style': style
        })
        return self.bandstructs[-1]
        
    def prepareBandstructure(self, bandstruct, basis, collapse_limit = 0.3):
        """
        This procedure will prepare a simple, scalar k-axis from the 
        vector k-s found in 'bandstruct'.
        """
        eig = bandstruct['bandstructure'].get_eigenvalues()
        kpoints = bandstruct['bandstructure'].get_kpoints(basis)
        
        kpointsCount = eig.shape[0]
        bandsCount = eig.shape[1]
        print "    Adding to plot: %i bands, %i kpoints" % (bandsCount, kpointsCount)
        
        k_vectors = []
        k_axis = []        
        bands = [[] for i in range(bandsCount)]
        dk_array = []
        k = 0.0
        
        if 'basis' in self.specialpoints: # If special points have been assigned
            if self.specialpoints['basis'][0].lower() == 'r' and basis[0].lower() == 'c':
                # Convert special points from reciprocal to cartesian basis:
                print "[Debug] Converting special points from rec to cart"
                specialPoints = { 'basis': 'cartesian', 'points': [] }
                S = bandstruct['bandstructure'].getReciprocalLatticeVectors()
                for sp in self.specialpoints['points']:
                    specialPoints['points'].append((sp[0],np.dot(S, sp[1])))
            elif self.specialpoints['basis'][0].lower() == 'c' and basis[0].lower() == 'r':
                # Convert special points from cartesian to reciprocal basis:
                print "[Debug] Converting special points from cart to rec"
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
        
        for idx,kvec in enumerate(kpoints):
            #print "kvec:",kvec
            k_vectors.append(kvec)
            if len(k_vectors) > 1:      # If not first k-vector:
                dkvec = kvec - k_vectors[len(k_vectors)-2]
                dk = np.sqrt(np.dot(dkvec,dkvec))
                if dk not in dk_array:
                    dk_array = np.append(dk_array,dk)
                if dk > collapse_limit:
                    print "NOTE: Collapsing"
                    dk = 0.001

                k += dk
            k_axis.append(k)
            bandNo = 0
            for eigval in eig[idx]:
                bands[bandNo].append(eigval - bandstruct['fermiLevel'])
                bandNo += 1
            
            # Check if we are at a special point:
            if 'points' in self.specialpoints: # If special points have been assigned
                for p in specialPoints['points']:
                    # WARNING:
                    # A difficulty arises in comparing k points, since 
                    # e.g. [0,-1,0] should be considered equivalent to [0,0,1].
                    # For simplicity, we just compare the length of the vectors.
                    # This seems to work fine, but it MAY give erroneous results!
                    #print "  ",np.sqrt(np.dot(kvec,kvec)),' || ',np.sqrt(np.dot(p[1],p[1]))
                    if (abs(np.sqrt(np.dot(kvec,kvec)) - np.sqrt(np.dot(p[1],p[1]))) < 1.e-4):
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
    
    def plotAndSaveToFile(self, outFile = 'bandstructure.pdf', silent_overwrite = False, **kwargs):
            

        import matplotlib
        matplotlib.use("pdf")
        from oppvasp.plotutils import prepare_canvas
        from matplotlib import rc
        rc('lines', linewidth=0.5)
        rc('font',**{'family':'serif','serif':['Palatino']})
        rc('font', family='serif')
        rc('text', usetex=True)
        prepare_canvas('10 cm')
        fig = plt.figure()
        p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        ax0 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        self.plot(ax0, **kwargs)
        # Write to file:
        if not silent_overwrite and os.path.exists(outFile):
            print "\nWARNING: The file \"%s\" already exists." % outFile
            if query_yes_no("         Do you want to overwrite it?") == 'no':
                return

        sys.stdout.write("\nSaving band structure to %s... " % outFile)
        sys.stdout.flush()
        plt.savefig(outFile)
        sys.stdout.write("done!\n\n")

    def plot(self, ax, yRange = 'auto', yTicks ='auto', basis = 'reciprocal', interpolate = False):
        """ 
        Parameters:
            ax: a matplotlib axis 
            basis: reciprocal or cartesian
        """

        print "[+] Starting plotting procedure"
        
        # Plot the bandstructures themselves:
        ymin = 0.
        ymax = 0.
        for bandstruct in self.bandstructs:
            
            self.prepareBandstructure(bandstruct, basis)
            
            kaxis = bandstruct['k_axis']
            xMax = kaxis[-1]              # temporary solution!!
            for i in range(len(bandstruct['E_axis'])):
                band = bandstruct['E_axis'][i]
                if len(kaxis) != len(np.unique(kaxis)):
                    raise Exception, "Uh oh, non-unique points on k-axis!"
                
                pltoptions = copy.copy(bandstruct['style'])
                #if (bandstruct['label'] != "" and i == 0):
                #    pltoptions['label'] =  bandstruct['label'] 
                if interpolate:
                    # Cubic-spline interpolation. We set s=0 for no smoothing.
                    # This ensures that the line passes through all points. 
                    # See http://www.tau.ac.il/~kineret/amit/scipy_tutorial/               
                    tck = scipy.interpolate.splrep(kaxis, band, s=0)                
                    xnew = np.arange(0,kaxis[-1],0.001)
                    ynew = scipy.interpolate.splev(xnew, tck, der=0)
                    pltoptions['markersize'] = 0
                    ax.plot(xnew,ynew, zorder=2, **pltoptions)
                    if 'markersize' in bandstruct['style']:
                        pltoptions['markersize'] = bandstruct['style']['markersize']
                    pltoptions['linewidth'] = 0
                    ax.plot(kaxis, band, zorder=3, **pltoptions)
                else:
                    print len(kaxis),len(band)
                    ax.plot(kaxis, band, zorder=2, **pltoptions)
                # Alternative: Plot data directly:
                # plt.plot(kaxis,band,'.')
                
                # Alternative 2:
                #xnew = np.arange(0,kaxis[-1],0.001)
                #intf = interp1d(xnew,band,2)
                #plt.plot(kaxis,intf(xnew),'b')
                if np.min(band) < ymin:
                    ymin = np.min(band)
                if np.max(band) > ymax:
                    ymax = np.max(band)

        # Plot horizontal and vertical lines behind the bandstructures:
        xticksX = []
        xticksLabels = []
        ax.axhline(y=0, color='grey', zorder=1)
        for i in range(len(self.k_axis_specialPoints)):
            ax.axvline(x=self.k_axis_specialPoints[i], color='silver', zorder=1)
            xticksX.append(self.k_axis_specialPoints[i])
            xticksLabels.append(self.k_axis_specialPointLabels[i])
        # Set ranges, ticks and legend:
        ax.set_xlim(0, xMax)
        if yRange == 'auto':
            yRange = [ymin,ymax]
            yTicks = [ymin,ymax,(ymax-ymin)/10.]
        print "YRANGE",yRange
        ax.set_ylim(yRange[0], yRange[1])
        ax.set_yticks(seq(yTicks[0],yTicks[1],yTicks[2]))
        if len(xticksX) != 0:
            ax.set_xticks(xticksX)
            ax.set_xticklabels(xticksLabels)
            ax.legend()
        

if __name__ == '__main__':
    #import oppvasp.vasp.parsebandstructure as vpb
    import oppvasp.espresso.parsebandstructure as epb
    
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
    
    #plot.addBands(vpb.ParseBandStructure('/Users/danmichael/Dropbox/code/oppvasp/src/tests/PROCAR','/Users/danmichael/Dropbox/code/oppvasp/src/tests/OUTCAR'), 
    #    fermiLevel = 6.54, color='blue', label = 'VASP', linestyle='dashed')
    plot.addBands(epb.ParseBandStructure('../oppvasp-examples/si.band','/Users/danmichael/Dropbox/code/oppvasp/src/tests/si.band.out'), 
        fermiLevel = 6.54, color='red', label = 'Quantum Espresso')
    plot.plotAndSaveToFile(outFile = 'bands.pdf', yMin = -13, yMax = 8.0, basis='cartesian') 

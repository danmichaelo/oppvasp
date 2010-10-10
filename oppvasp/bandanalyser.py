# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
#
# Created on 14. mars 2010
# @author: danmichael


import os,copy,sys
import numpy as np

class BandAnalyser():

    def __init__(self, bandstructure, specialpoints={}, fermiLevel = 0.00):
        
        # Initialise some arrays:
        self.bandstructure = bandstructure
        self.fermiLevel = fermiLevel
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

    def analyse(self, basis='reciprocal'):
        """
        This function will look for band gaps and the corresponding valence band
        maxima and conduction band minima. Perhaps we could also automaticly 
        calculate the bands' effective mass.
        """
        kpoints = self.bandstructure.getBands(basis)
        kpointsCount = len(kpoints)
        bandsCount = len(kpoints[0].eigenvals)
        print "[+] Band analyser input: %i bands, %i kpoints" % (bandsCount, kpointsCount)
        
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
        
        
        # Calculate band ranges:
        bandrange = range(bandsCount)
        bands = [[] for i in bandrange]
        ks = []
        for kpoint in kpoints:
            for bandNo in bandrange:
                bands[bandNo].append(kpoint.eigenvals[bandNo]['energy'] - self.fermiLevel)
            ks.append(kpoint)
        bandranges = [[] for i in range(bandsCount)]
        #print "Band ranges:"
        for bandNo in range(len(bands)):
            bandranges[bandNo] = (min(bands[bandNo]),max(bands[bandNo]))
            #print " - Band",bandNo,": ",bandranges[bandNo]


        # Create sets of overlapping bands:
        bandsets = []
        for i in bandrange:
            for j in range(i):
                if bandranges[i][0] < bandranges[j][1] and bandranges[j][0] < bandranges[i][1]:
                    # Band i and j overlap (how cute). Which band set are they in?
                    #print "Band",j,"overlaps",i
                    bset=-1
                    for k in range(len(bandsets)):
                        if i in bandsets[k] and j in bandsets[k]:
                            bset = k
                        if i in bandsets[k] and not j in bandsets[k]:
                            bset = k
                            bandsets[k].append(j)
                            #print "  Add %i to set %i" % (j,k)
                        elif j in bandsets[k] and not i in bandsets[k]:
                            bset = k
                            bandsets[k].append(i)
                            #print "  Add %i to set %i" % (i,k)
                    if bset == -1:
                        # i and j does not seem to belong to any existing band sets.
                        # Let's create a new one.
                        bandsets.append([i,j])
                        #print "  create new set"


        # Calculate band set ranges:
        bandsetranges = []
        for bs in bandsets:
            bsr = [0.0,0.0]
            for j in bs:
                br = bandranges[j]
                if bsr[0] == 0.0 and bsr[1] == 0.0:
                    bsr[0] = br[0]
                    bsr[1] = br[1]
                else:
                    if br[1] > bsr[1]:
                        bsr[1] = br[1]
                    if br[0] < bsr[0]:
                        bsr[0] = br[0]
            bandsetranges.append(bsr)
        
        vbset = -1 # valence band set
        cbset = -1 # conduction band set
        if len(bandsets) < 2:
            print "\n ERROR: No band gap found!\n"
        for b in range(len(bandsetranges)):
            if bandsetranges[b][1] <= 0:
                if vbset == -1:
                    vbset = b
                elif bandsetranges[b][1] > bandsetranges[vbset][1]:
                    vbset = b
            else:
                if cbset == -1:
                    cbset = b
                elif bandsetranges[b][0] < bandsetranges[vbset][0]:
                    cbset = b

        #print "Valence band set:",vbset," Conduction band set:",cbset
        vb0 = bandsetranges[vbset][1]
        cb0 = bandsetranges[cbset][0]
        
        vbands = [] # Valence bands
        cbands = [] # Conduction bands
        for i in bandsets[vbset]:
            if bandranges[i][1] == vb0:
                vbands.append(i)
        for i in bandsets[cbset]:
            if bandranges[i][0] == cb0:
                cbands.append(i)

        print "  Valence band maximum at",vb0,"eV consists of"
        for vband in vbands:
            bandMax = max((x,i) for i,x in enumerate(bands[vband]))
            k = ks[bandMax[1]]
            print "     band %i at %s" % (vband,k)

        print "  Conduction band minimum at",cb0,"eV consists of"
        for cband in cbands:
            bandMin = min((x,i) for i,x in enumerate(bands[cband]))
            k = ks[bandMin[1]]
            print "     band %i at %s" % (cband,k)

        print "Energy gap: ",(cb0-vb0),"eV" 

        #for kpoint in kpoints:
            #kvec = kpoint.getVector()
            #k_vectors.append(kvec)
            #if len(k_vectors) > 1:      # If not first k-vector:
                #dkvec = kvec - k_vectors[len(k_vectors)-2]
                #dk = np.sqrt(np.dot(dkvec,dkvec))
                #if dk not in dk_array:
                    #dk_array = np.append(dk_array,dk)
                #k += dk
            #k_axis.append(k)
            #bandNo = 0
            #for eigval in kpoint.eigenvals:
                #bands[bandNo].append(eigval['energy'] - bandstruct['fermiLevel'])
                #bandNo += 1
            
            ## Check if we are at a special point:
            #if 'points' in self.specialpoints: # If special points have been assigned
                #for p in specialPoints['points']:
                    ## WARNING:
                    ## A difficulty arises in comparing k points, since 
                    ## e.g. [0,-1,0] should be considered equivalent to [0,0,1].
                    ## For simplicity, we just compare the length of the vectors.
                    ## This seems to work fine, but it MAY give erroneous results!
                    #if np.sqrt(np.dot(kvec,kvec)) == np.sqrt(np.dot(p[1],p[1])):
                        #newPoint = True
                        #for op in self.k_axis_specialPoints:
                            #if abs(op-k) < 0.0001:
                                #newPoint = False
                        #if newPoint:
                            #print "     - Found ",p[0],"point (", kvec,") at",k
                            #self.k_axis_specialPoints.append(k)
                            #self.k_axis_specialPointLabels.append(p[0])


if __name__ == '__main__':
    print "Hello"

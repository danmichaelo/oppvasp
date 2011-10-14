# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4
#
# Created on 14. mars 2010
# @author: danmichael

import re, os, copy
import numpy as np
from oppvasp.kpoint import KPoint

class ParseBandStructure():
    '''
    This class will parse the band structure produced from PWSCF/Quantum Espresso
    into a standard format for use with 'plotbands.py'.
    '''

    def __init__(self, bandFile, bandOutFile = ''):
        """
        The 'bandFile' is the file produced by the post-processing code
        'bands.x' included with Quantum Espresso. This file contains the 
        energy eigenvalues for the selected k points. The k points are given 
        in Cartesian coordinates in units of a_0 (the lattice parameter).

        The 'bandOutFile' is a file containing the stdout from a 'band' 
        calculation carried out with 'pw.x'. This file is just used to read
        out the reciprocal lattice vectors, for conversion between 'cartesian'
        and 'reciprocal' coordinates.
        """
        self.bandstructure = []
        self.__successful = False

        print "[+] Parsing espresso bands from: %s" % bandFile
        
        self.__successful = self.readBandOutFile(bandOutFile)
        if self.__successful == True: 
            self.__successful = self.readBandFile(bandFile)

    def getBands(self, basis='cartesian'):
        """
        Returns the bands with the k-points either in cartesian coordinates or 
        in reciprocal coordinates, using b1,b2,b3 as basis.
        Quantum Espresso gives the k-vectors in the Cartesian basis, 2π/a(x,y,z),
        so if reciprocal basis is request, we transform the components to the 
        reciprocal basis (b1,b2,b3):
        """
        if basis == 'cartesian':
            return self.bandstructure
        elif basis == 'reciprocal':
            # Define transformation matrix to convert vectors from a 
            # cartesian (x,y,z) basis to a reciprocal (b1,b2,b3) basis. 
            # If x is a vector in the (x,y,z) basis and x' a vector
            # in the (b1,b2,b3) basis, the components are related by
            # x' = S^-1 x            
            S = self.reciprocalAxes
            invS = np.linalg.inv(S)
            b = copy.deepcopy(self.bandstructure)
            for kpoint in b:                
                # dot(A,v) treats v as a column vector:
                kpoint.setVector(np.dot(invS,kpoint.getVector())) 
            return b
        else:
            raise StandardError("Unknown basis requested")

    def getDirectLatticeVectors(self):
        return self.crystalAxes
    
    def getReciprocalLatticeVectors(self):
        return self.reciprocalAxes
    
    def isValid(self):
        return self.__successful
        
    def readBandOutFile(self, bandOutFile):    
        """
        This function will parse the std out file from a pw.x band calculation run,
        in order to find the definitions of the direct and reciprocal lattice vectors. 

        Example of what we are looking for::
        
              crystal axes: (cart. coord. in units of a_0)
                        a(1) = (  0.000000  0.500000  0.500000 )  
                        a(2) = (  0.500000  0.000000  0.500000 )  
                        a(3) = (  0.500000  0.500000  0.000000 ) 

              reciprocal axes: (cart. coord. in units 2 pi/a_0)
                        b(1) = ( -1.000000  1.000000  1.000000 )  
                        b(2) = (  1.000000 -1.000000  1.000000 )  
                        b(3) = (  1.000000  1.000000 -1.000000 )
        """
        crystalAxisHeaderMatcher = re.compile("\s*crystal axes:")
        crystalAxisMatcher = re.compile("\s*a\(\d\)\s*=[(\s]+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)")
        reciprocalAxisHeaderMatcher = re.compile("\s*reciprocal axes:")
        reciprocalAxisMatcher = re.compile("\s*b\(\d\)\s*=[(\s]+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)")
                
        self.crystalAxes = np.zeros((3,3))
        self.reciprocalAxes = np.zeros((3,3))

        if not os.path.exists(bandOutFile):
            print """"
            FATAL ERROR: 
            The specified input file was not found:
            "%s"
            """ % (bandOutFile)
            return False
        
        crystalAxes = []
        reciprocalAxes = []
        f = open(bandOutFile,'r')        
        lines = f.readlines()
        for i in xrange(len(lines)):

            if len(crystalAxes) == 0 and crystalAxisHeaderMatcher.match(lines[i]):
                # Good, we found the crystal axes header. Now lets check
                # if the next three lines are the crystal axes coordinates:
                for j in range(i+1,i+4):
                    reg = crystalAxisMatcher.match(lines[j])
                    if (reg):
                        c = [reg.group(1),reg.group(2),reg.group(3)]
                        c = [float(ci) for ci in c]
                        crystalAxes.append(c)

            elif len(crystalAxes) != 0 and reciprocalAxisHeaderMatcher.match(lines[i]):
                # Good, we found the reciprocal axes header. Now lets check
                # if the next three lines are the reciprocal axes coordinates:
                for j in range(i+1,i+4):
                    reg = reciprocalAxisMatcher.match(lines[j])
                    if (reg):
                        c = [reg.group(1),reg.group(2),reg.group(3)]
                        c = [float(ci) for ci in c]
                        reciprocalAxes.append(c)
                break
            
        if len(crystalAxes) == 3 and len(reciprocalAxes) == 3:
            self.crystalAxes = np.transpose(np.array(crystalAxes))
            self.reciprocalAxes = np.transpose(np.array(reciprocalAxes))
            return True
        else:
            print """
            FATAL ERROR: 
            Complete definitions of the crystal and reciprocal axes 
            were not found in the specified pw.x output file: 
            "%s"
            
            Found %i of 3 crystal axes and %i of 3 reciprocal axes.
            
            Please check that the file is a valid output file. If it is, it comes 
            from a version of pw.x not yet supported.
            """ % (bandOutFile,len(crystalAxes),len(reciprocalAxes))
            return False


    def readBandFile(self, bandFile):
        
        if not os.path.exists(bandFile):
            print """"
            FATAL ERROR: The specified input file was not found:
            \"%s\"
            """ % (bandFile)
            return False
        
        f = open(bandFile,'r')        
        
        # Get number of bands and number of k points from first line:
        firstline = f.readline()
        firstline = re.match(".*nbnd=[\s]*(?P<nbnd>[\d]+),[\s]*nks=[\s]*(?P<nks>[\d]+)",firstline)
        if firstline:
            nbnd = int(firstline.group('nbnd'))
            nks = int(firstline.group('nks'))
        else:
            print """
            FATAL ERROR:
            The first line of the band file was not recognized as a valid header line.
            Please check that it is a valid band file.
            """
            return False
        print "    Header reports %i bands and %i k-points" % (nbnd,nks)
        
        # Read rest of the file
        kpoints = []
        isKPointLine = True
        bandNo = 0
        for line in f.readlines():
            if isKPointLine:
                kvec = np.array([float(i) for i in line.split()])
                kpoint = KPoint(kvec[0],kvec[1],kvec[2],[])
                isKPointLine = False
            else:
                for val in line.split():
                    bandNo += 1
                    kpoint.appendEigenval({ 'energy': float(val) })
                if bandNo == nbnd:
                    kpoints.append(kpoint)
                    isKPointLine = True
                    bandNo = 0
        f.close()
        
        print "    File parsed successfully. Found %i k-points" % (len(kpoints))
        
        # Check for variation in k-point density
        #kdmin = np.min(dk_array)
        #kdmax = np.max(dk_array)
        #kdtreshold = 0.01
        #if np.abs(kdmax-kdmin) > kdtreshold:
        #    print "NOTICE: k-point density varies from %.3f to %.3f" %(1/kdmax, 1/kdmin)


        # Plot it, baby!
        #pb = PlotBandStructure(k_axis, bands, k_axis_specialpoints)
        #pb.plotAndSaveToFile(outFile, yMin, yMax)
        
        self.bandstructure = kpoints        
        return True
                

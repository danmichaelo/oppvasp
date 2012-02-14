# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4
#
# Created on 14. mars 2010
# @author: danmichael

import re,os,copy
import numpy as np
from oppvasp.kpoint import KPoint
from oppvasp.bandstructure import BandStructure
from oppvasp.vasp import parsers

def parse_espresso_bands(procar = '', eigenval='EIGENVAL', outcar = 'OUTCAR', vasprun = ''):
    """
    Specify either vasprun.xml, eigenval+outcar or procar+outcar

    Returns:
        oppvasp.bandstructure.Bandstructure object
    """
    parser = VaspBandParser(procar,eigenval,outcar,vasprun)
    return parser.get_band_structure()

class VaspBandParser():

    def __init__(self, procar = '', eigenval='EIGENVAL', outcar = 'OUTCAR', vasprun = ''):
        self.bandstructure = []
        self.__successful = False

        print "[+] Parsing bandstructure from VASP"
        
        #self.__successful = self.readOUTCAR(outcar)
        #if self.__successful == True: 
        #    self.__successful = self.readEIGENVAL(eigenval)
        #if self.__successful == True and procar != '': 
        #    self.__successful = self.readPROCAR(procar)
        if vasprun != '':
            vp = parsers.VasprunParser(vasprun)
            self.k = vp.get_kpoints()
            self.eig = vp.get_eigenvalues()
    
    def get_band_structure(self):
        """
        VASP gives the k-vectors in the cartesian basis, in units of ??? 
        """
        #bands = []
        #for s in range(self.eig.shape[0]):
        #    sbands = []
        #    for idx,k in enumerate(self.k):
        #        eigs = self.eig[s,:,idx]
        #        sbands.append(KPoint(k, eigenvals = eigs))
        #    bands.append(sbands)

        return BandStructure(eigenvalues = self.eig, kpoints = self.k coords = 'cart', basis = self.crystalAxes)
    
    
    def isValid(self):
        """ DEPRECATED"""
        return self.__successful
    
    #def getBands(self, basis='reciprocal'):
    #    """
    #    Returns the bands with the k-points either in cartesian coordinates or 
    #    in reciprocal coordinates, using b1,b2,b3 as basis.
    #    VASP gives the k-vectors in the reciprocal basis, (b1,b2,b3).
    #    """
    #    if basis[0].lower() == 'r':
    #        return self.bandstructure
    #    elif basis[0].lower() == 'c':
    #        # Define transformation matrix to convert vectors from a 
    #        # reciprocal basis to a cartesian basis. 
    #        # If x is a vector in the (x,y,z) basis and x' a vector
    #        # in the (b1,b2,b3) basis, the components are related by
    #        # x = S x'
    #        # We round off after four digits to avoid numerical noise
    #        # preventing us from finding special points
    #        S = self.reciprocalAxes
    #        b = copy.deepcopy(self.bandstructure)
    #        for kpoint in b:                
    #            # dot(A,v) treats v as a column vector:
    #            kpoint.setVector(np.dot(S,kpoint.getVector())) 
    #        return b
    #    else:
    #        raise StandardError("Unknown basis requested")

    
    def readKPointFromPROCAR(self, stream, bandsCount, ionsCount):
        
        kpointHeaderMatcher = re.compile("\s*k-point[^:]+:\s*(-?[0-9.]{10})\s*(-?[0-9.]{10})\s*(-?[^\s]{10})")
        # Example line to match:
        #  k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.05000000
        #  k-point    8 :   -0.37500000-0.37500000 0.12500000
        #  Note that above there is no space between two of the entries. Therefore we
        #  have to assume 10 digits is always used for each coordinate.

        bandHeaderMatcher = re.compile(".+energy\s+([^#]+)# occ\.\s(.+)$")
        # Example line to match:
        # band   2 # energy    6.52193506 # occ.  2.00316155
        
        orbitalHeaderMatcher = re.compile("ion\s+s.+tot$")
        # Example line to match:
        # ion      s      p      d    tot
        
        line = stream.readline()
        while not kpointHeaderMatcher.match(line):
            line = stream.readline()
        kpointHeaderLine = line
        m = kpointHeaderMatcher.match(kpointHeaderLine)
        
        print '   ',m.group(0);
        
        kpoint_x = float(m.group(1))
        kpoint_y = float(m.group(2))
        kpoint_z = float(m.group(3))


        bands = []
        for i in range(0,bandsCount):
            line = stream.readline()
            while not bandHeaderMatcher.match(line):
                line = stream.readline()
            bandHeaderLine = line
            energy = float(bandHeaderMatcher.match(line).group(1))

            line = stream.readline()
            while not orbitalHeaderMatcher.match(line):
                line = stream.readline()
            orbitalHeaderLine = line
            self.orbitalNames = orbitalHeaderLine.split()[1:-1]

            orbitals = []
            for j in range(0,ionsCount):
                line = stream.readline()
                orbitalAssignmentText = line.split()
                orbitals.append(orbitalAssignmentText)

            bands.append({ 'energy' : energy, 'orbitals' : orbitals })

        return KPoint(kpoint_x, kpoint_y, kpoint_z, bands)
    
    def readKPointFromEIGENVAL(self, lines):
        
        kpointHeader = lines[0].split()
        # Example line:
        # 0.5000000E+00  0.5000000E+00  0.5000000E+00  0.2000000E-01
        
        kpoint_x = float(kpointHeader[0])
        kpoint_y = float(kpointHeader[1])
        kpoint_z = float(kpointHeader[2])

        bands = []
        for line in lines[1:]:
            lineWords = line.split()
            energy = float(lineWords[1])
            bands.append({ 'energy' : energy, 'orbitals' : [] })

        return KPoint(kpoint_x, kpoint_y, kpoint_z, bands)

    def getMetadataFromPROCAR(self,line):
        metaMatcher = re.compile("#[^:]+:([^#]+)#[^:]+:([^#]+)#[^:]+:(.+)$")
        m = metaMatcher.match(line)
        kpointsCount = int(m.group(1))
        bandsCount = int(m.group(2))
        ionsCount = int(m.group(3))
        return (kpointsCount, bandsCount, ionsCount)
    
    def getMetadataFromEIGENVAL(self,line):
        metaMatcher = re.compile("^\s*(\d*)\s*(\d*)\s*(\d*)")
        m = metaMatcher.match(line)
        # What is m.group(1) ?
        kpointsCount = int(m.group(2))
        bandsCount = int(m.group(3))
        return (kpointsCount, bandsCount)

    def readPROCAR(self,procar):   

        if not os.path.exists(procar):
            print """
            FATAL ERROR: The specified input file was not found:
            "%s"
            """ % (procar)
            return False

        f = open(procar, "r")
        f.readline() # throwaway
        metaLine = f.readline() # metadata
        kpointsCount, bandsCount, ionsCount = self.getMetadataFromPROCAR(metaLine)
        print "    PROCAR header reports %i bands, %i kpoints and %i ions." % (bandsCount, kpointsCount, ionsCount)

        kpoints = []
        for i in range(0,kpointsCount):
            kpoint = self.readKPointFromPROCAR(f, bandsCount, ionsCount)
            if len(kpoints) != 0 and kpoint == kpoints[-1]:
                print "     - Duplicate k point found at ",kpoint
            else:
                kpoints.append(kpoint)
        print "    File parsed successfully. %i inequivalent kpoints found." % (len(kpoints))

        self.bandstructure = kpoints
        return True

    def readEIGENVAL(self,eigenval):   

        if not os.path.exists(eigenval):
            print """
            FATAL ERROR: The specified input file was not found:
            "%s"
            """ % (eigenval)
            return False

        f = open(eigenval, "r")
        lines = f.readlines()

        kpointsCount, bandsCount = self.getMetadataFromEIGENVAL(lines[5])
        print "    EIGENVAL header reports %i bands and %i kpoints." % (bandsCount, kpointsCount)

        kpoints = []
        for i in range(0,kpointsCount):
            kpoint = self.readKPointFromEIGENVAL(lines[(7+i*(bandsCount+2)):(6+(i+1)*(bandsCount+2))])
            if len(kpoints) != 0 and kpoint == kpoints[-1]:
                print "     - Duplicate k point found at ",kpoint
            else:
                kpoints.append(kpoint)
        print "    File parsed successfully. %i inequivalent kpoints found." % (len(kpoints))

        self.bandstructure = kpoints
        return True

    def readOUTCAR(self, outcar):    
        """
        This function will parse the outcar file from a band calculation run,
        in order to find the definitions of the direct and reciprocal lattice vectors. 

        Example of what we are looking for:

        LATTYP: Found a face centered cubic cell.
        ALAT       =     5.3893600000
         
        [...]

        direct lattice vectors                 reciprocal lattice vectors
        0.000000000  2.694680000  2.694680000    -0.185550789  0.185550789  0.185550789
        2.694680000  0.000000000  2.694680000     0.185550789 -0.185550789  0.185550789
        2.694680000  2.694680000  0.000000000     0.185550789  0.185550789 -0.185550789
        """
        alatMatcher = re.compile("\s*ALAT\s*=\s*([.\d]+)")
        axesHeaderMatcher = re.compile("\s*direct lattice vectors\s+reciprocal lattice vectors")
        axisMatcher = re.compile("\s*([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)")
                
        self.crystalAxes = np.zeros((3,3))
        self.reciprocalAxes = np.zeros((3,3))

        if not os.path.exists(outcar):
            print """
            FATAL ERROR: 
            The specified input file was not found:
            "%s"
            """ % (outcar)
            return False
        
        latticeParameters = {}
        crystalAxes = []
        reciprocalAxes = []
        f = open(outcar,'r')        
        lines = f.readlines()
        for i in xrange(len(lines)):

            if len(latticeParameters) == 0:
                m = alatMatcher.match(lines[i])
                if m:
                    # Good, we found the lattice parameter
                    # We save it in a flexible format in case we would
                    # like to include other lattice parameters later
                    latticeParameters = {'a': float(m.group(1)) }
            elif axesHeaderMatcher.match(lines[i]):
                # Good, we found the axes header. Now lets check
                # if the next three lines are the axes coordinates:
                for j in range(i+1,i+4):
                    m = axisMatcher.match(lines[j])
                    if m:
                        d = [float(m.group(i)) for i in range(1,4)]
                        r = [float(m.group(i)) for i in range(4,7)]
                        crystalAxes.append(d)
                        reciprocalAxes.append(r)
                break
            
        if len(crystalAxes) == 3 and len(reciprocalAxes) == 3:
            self.latticeParameters = latticeParameters

            # We round off after four digits to avoid numerical noise
            # preventing us from finding special points:
            self.crystalAxes = np.round(np.transpose(np.array(crystalAxes)/self.latticeParameters['a']),4)
            self.reciprocalAxes = np.round(np.transpose(np.array(reciprocalAxes)*self.latticeParameters['a']),4)
            
            return True
        else:
            print """
            FATAL ERROR: 
            Complete definitions of the crystal and reciprocal axes 
            were not found in the specified OUTCAR file: 
            "%s"
            
            Found %i of 3 crystal axes and %i of 3 reciprocal axes.
            
            Please check that the file is a valid OUTCAR file. If it is, it comes 
            from a version of VASP not yet supported.
            """ % (outcar,len(crystalAxes),len(reciprocalAxes))
            return False

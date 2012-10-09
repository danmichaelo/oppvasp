# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4
#
# Created on 14. mars 2010
# @author: danmichael

import re, os, copy
import numpy as np
from oppvasp.kpoint import KPoint
from oppvasp.bandstructure import BandStructure


def parse_espresso_bands(band_file, out_file):
    """
    Parses band structure(s) from PWSCF/Quantum Espresso into a
    standard representation.

    Parameters:

        band_file : a single band file
            produced by bands.x (by setting filband). should 
            contain the energy eigenvalues for the selected k points,
            the k points given in Cartesian coordinates in units of a_0 
            (the lattice parameter).

        out_file: A file containg the stdout from a pw.x 'band' calculation. 
            This is just for getting the lattice vectors. 
            A more elegant solution may be implemented at some time...

    Returns:
        
        oppvasp.bandstructure.BandStructure object
    """

    print "[+] Parsing espresso bands from: %s" % band_file
    kpoints,eig = read_band_file(band_file)
    crystalAxes, reciprocalAxes = read_out_file(out_file)
    # Quantum Espresso gives the k-vectors in the cartesian basis, in units of a_0
    return BandStructure( kpoints = kpoints, eigenvalues = eig, coords = 'cart', basis = crystalAxes)

def read_out_file(filename):    
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
            
    crystalAxes = np.zeros((3,3))
    reciprocalAxes = np.zeros((3,3))

    if not os.path.exists(filename):
        print """"
        FATAL ERROR: 
        The specified input file was not found:
        "%s"
        """ % (filename)
        return False
    
    crystalAxes = []
    reciprocalAxes = []
    f = open(filename,'r')        
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
        crystalAxes = np.transpose(np.array(crystalAxes))
        reciprocalAxes = np.transpose(np.array(reciprocalAxes))
        return crystalAxes, reciprocalAxes
    else:
        print """
        FATAL ERROR: 
        Complete definitions of the crystal and reciprocal axes 
        were not found in the specified pw.x output file: 
        "%s"
        
        Found %i of 3 crystal axes and %i of 3 reciprocal axes.
        
        Please check that the file is a valid output file. If it is, it comes 
        from a version of pw.x not yet supported.
        """ % (filename,len(crystalAxes),len(reciprocalAxes))
        return False


def read_band_file(filename):
    
    if not os.path.exists(filename):
        print """"
        FATAL ERROR: The specified input file was not found:
        \"%s\"
        """ % (filename)
        return False
    
    f = open(filename,'r')        
    
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
    eigenvalues = []
    isKPointLine = True
    bandNo = 0
    for line in f.readlines():
        if isKPointLine:
            kvec = np.array([float(i) for i in line.split()])
            kpoints.append(kvec)
            eigenvalues.append([])
            #kpoint = KPoint(kvec[0],kvec[1],kvec[2],[])
            isKPointLine = False
        else:
            for val in line.split():
                bandNo += 1
                eigenvalues[-1].append(float(val))
            if bandNo == nbnd:
                isKPointLine = True
                bandNo = 0
    f.close()
    
    print "    File parsed successfully. Found %i k-points" % (len(kpoints))
    
    return np.array(kpoints), np.array(eigenvalues)


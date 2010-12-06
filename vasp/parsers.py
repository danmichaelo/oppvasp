# This file contains code originally written by olem

import sys,re, math

#from xml.dom.ext.reader import Sax2
#from xml import xpath

from lxml import etree
from StringIO import StringIO

class vasprunParser:
    """
    Parser for vasprun.xml files, making use of libxml for relatively fast parsing.
    """
    
    def __init__(self, filename = 'vasprun.xml'):

        docstr = open(filename).read()

        # zap control characters that invalidates the xml
        docstr = re.sub('[\x00-\x09\x0B-\x1F]','',docstr)

        parser = etree.XMLParser()
        try:
            self.doc = etree.parse(StringIO(docstr), parser)
        except etree.XMLSyntaxError:
            print "Failed to parse xml file: ",filename
            error = parser.error_log[0]
            print(error.message)
            sys.exit(2)
    
    def getIncarProperty(self, propname):
        """ 
        Returns the value of a given INCAR property as a string,
        or throws a LookupError if the property was not found.
        Example: 
        >>> getIncarProperty('ENCUT')
        """
        results = self.doc.xpath( "/modeling/incar/i[@name='"+propname+"']")
        if results:
            return results[0].text
        else:
            raise LookupError('Value not found')
    
    def getTotalEnergy(self):
        """
        Returns the total energy in electronvolt
        """
        results = self.doc.xpath( "/modeling/calculation/energy/i[@name='e_fr_energy']")
        if results:
            return float(results[0].text)
        else:
            raise LookupError('Value not found')
    
    def getFinalVolume(self):
        """
        Returns the final volume in units Angstrom^3
        """
        results = self.doc.xpath( "/modeling/structure[@name='finalpos']/crystal/i[@name='volume']")
        if results:
            return float(results[0].text)
        else:
            raise LookupError('Value not found')

    def getSCsteps(self):
        """
        Returns array of electronic self-consistent steps
        """
        results = self.doc.xpath( "/modeling/calculation/scstep")
        if results:
            return results
        else:
            raise LookupError('Value not found')



class poscarParser:
    """
    Parser for POSCAR files
    """
    
    def __init__(self, poscarname='POSCAR'):
        self.selective = 0 
        self.readPOSCAR(poscarname)

    def readPOSCAR(self,poscarname):
        poscarfile = open( poscarname, 'r')  # r for reading
        commentline = poscarfile.readline()
        scaleline = poscarfile.readline()
        vec1line = poscarfile.readline()
        vec2line = poscarfile.readline()
        vec3line = poscarfile.readline()
        sixthline = poscarfile.readline()  # Test for vasp5 syntax
        try:
            dummy = int(sixthline.split()[0])
            atomnumberline = sixthline
        except:
            atomnumberline = poscarfile.readline()
        atomnumbers = map(int,atomnumberline.split())
        natoms = 0
        for i in range(len(atomnumbers)):
            natoms = natoms + atomnumbers[i]
        seventhline = poscarfile.readline()

        if seventhline[0] == 'S' or seventhline[0] == 's':
            self.selective = 1
            seventhline = poscarfile.readline()
        else:
            self.selective = 0

        if self.selective:
           x = []; y = []; z = []
           for j in range(natoms):
            line = poscarfile.readline()  # read a line
            if not line: break
            xyz = line.split()
            x.append(xyz[3])
            y.append(xyz[4])
            z.append(xyz[5])
        
        poscarfile.close()

class outcarParser:
    """
    Parser for OUTCAR files
    """

    def getIncarProperty(self, propname):
        outfile = open(self.filename, 'r')
        lines = outfile.readlines()
        s = re.compile('[\t ]*'+propname+'[\t ]*=[\t ]*([0-9.]*)')
        for l in lines:
            res = s.match(l)
            if res:
                return res.group(1)

        print "Failed to lookup INCAR property "+propname+" in "+self.filename
        sys.exit(1)


    def __init__(self, outcarname = 'OUTCAR', selective = 0):

        self.filename = outcarname

        self.kpoints = 0
        self.dist = 0
        self.toten = 0
        self.planewaves = 0
        self.cpu = 0
        self.memory = 0
        self.maxforce = 0
        self.maxdrift = 0
        self.maxpressure = 0

        outfile = open(outcarname, 'r')
        while 1:
            line = outfile.readline()
            kmatch = re.search('(\d+) +irreducible', line)
            ematch = re.search('free  energy   TOTEN  = +(-*\d+.\d+)', line)
            cpumatch = re.search('Total CPU time used \(sec\): +(\d+.\d+)', line)
            distmatch = re.search('Following cartesian coordinates:', line)

            #k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=     96
            #number of dos      NEDOS =    301   number of ions     NIONS =      8
            #non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8
            #total plane-waves  NPLWV =  32768

            planewavematch = re.search('NPLWV[ \t]*=[ \t]*([0-9])', line)
            nbandsmatch = re.search('NPLWV[ \t]*=[ \t]*([0-9])', line)
            if planewavematch:
                self.planewaves = int(planewavematch.group(1))
            if kmatch:
                self.kpoints = int(kmatch.group(1))
            elif ematch:
                self.toten= float(ematch.group(1))
            elif cpumatch:
                self.cpu = float(cpumatch.group(1))
            elif distmatch:
                if self.kpoints > 1:
                    tmpline = outfile.readline()
                    firstline = outfile.readline()
                    secondline = outfile.readline()
                    k1x,k1y,k1z,dummy = map(float, firstline.split())
                    k2x,k2y,k2z,dummy = map(float, secondline.split())
                    self.dist = math.sqrt((k2x-k1x)*(k2x-k1x)+(k2y-k1y)*(k2y-k1y)+(k2z-k1z)*(k2z-k1z))
                else:
                    self.dist = 0
            elif re.search(r'external pressure', line): 
                tmp,tmp,tmp,pressure,tmp,tmp,tmp,tmp,tmp,tmp = line.split()
                self.maxpressure = float(pressure)
            elif re.search(r'TOTAL\-FORCE', line):
                i=0
                line = outfile.readline()
                maxdrift = 0.0
                maxforce = 0.0
                while 1:
                    line = outfile.readline()
                    if re.search(r'----', line):
                        line = outfile.readline()
                        a,b,driftx,drifty,driftz = line.split()
                        if abs(float(driftx)) > maxdrift:
                            maxdrift = abs(float(driftx))
                        if abs(float(drifty)) > maxdrift:
                            maxdrift = abs(float(drifty))
                        if abs(float(driftz)) > maxdrift:
                            maxdrift = abs(float(driftz))
                        break
                    posx,posy,posz,forx,fory,forz = map(float, line.split())
                    if selective:
                        if (abs(forx) > maxforce) and (x[i] == 'T' or x[i] == 't'):
                            maxforce = abs(forx)
                            maxi = i
                        if (abs(fory) > maxforce) and (y[i] == 'T' or y[i] == 't'):
                            maxforce = abs(fory)
                            maxi = i
                        if (abs(forz) > maxforce) and (z[i] == 'T' or z[i] == 't'):
                            maxforce = abs(forz)
                            maxi = i
                    else:
                        if abs(forx) > maxforce:
                            maxforce = abs(forx)
                            maxi = i
                        if abs(fory) > maxforce:
                            maxforce = abs(fory)
                            maxi = i
                        if abs(forz) > maxforce:
                            maxforce = abs(forz)
                            maxi = i
                        i = i+1
                self.maxforce = maxforce
                self.maxdrift = maxdrift
            if not line:
                break
        outfile.close()

    def getTotalEnergy(self):
        return self.toten

    def getCPUTime(self):
        return self.cpu

    #def read_stress(self):
    #    for line in open('OUTCAR'):
    #        if line.find(' in kB  ') != -1:
    #            stress = -np.array([float(a) for a in line.split()[2:]]) \
    #                     [[0, 1, 2, 4, 5, 3]] \
    #                     * 1e-1 * ase.units.GPa
    #    return stress
        


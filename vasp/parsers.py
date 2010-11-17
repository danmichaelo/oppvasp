# This file contains code originally written by olem

from xml.dom.ext.reader import Sax2
from xml import xpath

class vasprunParser:
        
    
    def __init__(self,vasprun='vasprun.xml'):

        self.readXml(vasprun)

    def readXml(self,filename):

        doc = Sax2.FromXmlFile(filename).documentElement
        results = xpath.Evaluate( "/modeling/calculation/energy/i[@name='e_fr_energy']", doc)
        if results:
            self.toten = float(results[0].firstChild.nodeValue)
        
        print '%.6f' % (self.toten)


class poscarParser:
    
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

    def __init__(self, outcarname = 'OUTCAR', selective = 0):

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


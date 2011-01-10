"""
This file contains classes to easily extract information from VASP input
and output files. Functions are added as needed, and only a small 
subset of the total information available in the input and output files can 
currently be extracted using the classes in this file.

When time allows for it, functions to extract the same information from
both vasprun.xml and OUTCAR type files are implemented. Usually, extracting
data from vasprun.xml files are easier due to the excellent Python xml
parsers available, but it is believed that extracting data from OUTCAR files 
may be faster for very large files. 

Note that this file contains code originally written by olem
"""

import sys,re,math,os

# Useful notes on libxml: http://codespeak.net/lxml/parsing.html
from lxml import etree
from StringIO import StringIO
import numpy as np

class VasprunParser:
    """
    Parser for vasprun.xml files, making use of libxml for relatively fast parsing.
    """
    
    def __init__(self, filename = 'vasprun.xml', verbose = False):
        
        if verbose:
            print "Parsing %s (%.2f MB)... " % (filename,os.path.getsize(filename)/1024.**2)

        self.filename = filename
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
    
    def get_incar_property(self, propname):
        """ 
        Returns the value of a given INCAR property as a string,
        or throws a LookupError if the property was not found.
        Example: 
        >>> get_incar_property('ENCUT')
        """
        results = self.doc.xpath( "/modeling/incar/i[@name='"+propname+"']")
        if results:
            return results[0].text
        else:
            raise LookupError('Value not found')

    def get_num_kpoints(self):
        """Returns the number of k-points"""
        results = self.doc.xpath( "/modeling/kpoints/varray[@name='"+kpointlist+"']")
        if results:
            return results[0].text
        else:
            raise LookupError('Value not found')    
    
    def get_total_energy(self):
        """
        Returns the total energy in electronvolt
        """
        results = self.doc.xpath( "/modeling/calculation/energy/i[@name='e_fr_energy']")
        if results:
            return float(results[0].text)
        else:
            raise LookupError('Value not found')
    

    def get_sc_steps(self):
        """
        Returns array of electronic self-consistent steps
        """
        results = self.doc.xpath( "/modeling/calculation/scstep")
        if results:
            return results
        else:
            raise LookupError('Value not found')

    def get_force_on_atom(self,atom_no):
        """
        Returns the force on atom <atom_no> as a numpy vector (x,y,z),
        where 1 is the first atom.
        """
        forces = self.doc.xpath( "/modeling/calculation/varray[@name='forces']/v")
        force = np.array([float(f) for f in forces[atom_no-1].text.split()])
        return force

    def get_initial_positions(self):
        """
        Returns the final position of all atoms as a (n,3) numpy array, where n is the number of atoms
        """
        all_pos = self.doc.xpath( "/modeling/structure[@name='initialpos']/varray[@name='positions']/v")
        num_atoms = len(all_pos)
        pos_array = np.zeros((num_atoms,3))
        for i in range(num_atoms):
            pos_array[i] = [float(f) for f in all_pos[i].text.split()]
        return pos_array

    def get_final_positions(self):
        """
        Returns the final position of all atoms as a (n,3) numpy array, where n is the number of atoms
        """
        all_pos = self.doc.xpath( "/modeling/structure[@name='finalpos']/varray[@name='positions']/v")
        num_atoms = len(all_pos)
        pos_array = np.zeros((num_atoms,3))
        for i in range(num_atoms):
            pos_array[i] = [float(f) for f in all_pos[i].text.split()]
        return pos_array
    
    def get_initial_velocities (self):
        """
        Returns the initial velocities of all atoms as a (n,3) numpy array, where n is the number of atoms
        """
        all_vel = self.doc.xpath( "/modeling/structure[@name='initialpos']/varray[@name='velocities']/v")
        if not all_vel:
            raise LookupError('Velocities not found. Is this file from a MD run?')

        num_atoms = len(all_vel)
        vel_array = np.zeros((num_atoms,3))
        for i in range(num_atoms):
            vel_array[i] = [float(f) for f in all_vel[i].text.split()]
        return vel_array

    def get_final_velocities(self):
        """
        Returns the final velocities of all atoms as a (n,3) numpy array, where n is the number of atoms
        """
        all_vel = self.doc.xpath( "/modeling/structure[@name='finalpos']/varray[@name='velocities']/v")
        if not all_vel:
            raise LookupError('Velocities not found. Is this file from a MD run?')

        num_atoms = len(all_vel)
        vel_array = np.zeros((num_atoms,3))
        for i in range(num_atoms):
            vel_array[i] = [float(f) for f in all_vel[i].text.split()]
        return vel_array

    def get_final_atom_position(self,atom_no):
        """
        Returns the final position of atom <atom_no> as a numpy vector (x,y,z)
        """
        pos = self.doc.xpath( "/modeling/structure[@name='finalpos']/varray[@name='positions']/v")
        atpos = np.array([float(f) for f in pos[atom_no-1].text.split()])
        return atpos

    def get_final_volume(self):
        """
        Returns the final volume in units Angstrom^3
        """
        results = self.doc.xpath( "/modeling/structure[@name='finalpos']/crystal/i[@name='volume']")
        if results:
            return float(results[0].text)
        else:
            raise LookupError('Value not found')

    def get_max_force(self):
        """
        Returns the max force acting on any atom
        """
        forces = self.doc.xpath( "/modeling/calculation/varray[@name='forces']/v")
        max_force = 0.
        for f in forces:
            force = np.array(f.text.split())
            force_norm = np.sqrt(np.dot(force,force))
            if force_norm > max_force:
                max_force = force_norm
        return max_force

    def get_cpu_time(self):
        """
        Returns the CPU time spent. The value returned corresponds to the the value found in
        OUTCAR files on a line like this: 
        
        >>> LOOP+:  cpu time  482.23: real time  482.58
        
        This will value is always somewhat lower than the value found in OUTCAR file on a 
        line like this:
        
        >>>         Total CPU time used (sec):      490.877
        
        but it appears like this value is not available from vasprun.xml files.
        """
        time = self.doc.xpath( "/modeling/calculation/time[@name='totalsc']")
        if time:
            time = time[0].text.split()
            return float(time[0])
        else:
            raise LookupError('Value not found in %s' % (self.filename))
        



class FileIterator:
    """
    Abstract iterator for reading files
    """

    def __init__(self, filename, cache = True):
        """
        If <cache> is set to True, the whole file is read in at once. 
        This is usually faster, but may not be possible for very large files.
        """
        self.filename = filename
        self.file = open(filename,'r')
        self.cached = cache
        # Caching is preferred unless the file is too big too keep in memory 
        if cache:
            self.contents = self.file.readlines()
            self.numlines = len(self.contents)
        self.reset()

    def __iter__(self):
        return self

    def next(self):
        self.lineno += 1
        if self.cached and self.lineno >= self.numlines:
            raise StopIteration
        elif not self.cached:
            try:
                return "1 "+self.file.readline()
            except GeneralError:
                raise StopIteration
        else:
            return "0 "+self.contents[self.lineno]

    def reset(self):
        self.lineno = -1
        if not self.cached:
            self.file.seek(0)


class OutcarParser:
    """
    Parser for OUTCAR files
    """

    def __init__(self, outcarname = 'OUTCAR', selective = 0, verbose = False):
        
        if verbose:
            print "Parsing %s (%.1f MB)... " % (outcarname,os.path.getsize(outcarname)/1024.**2)

        self.filename = outcarname
        self.file = fileIterator(self.filename)
        self.selective = selective
        
        # Read the first lines to find the following parameters:
        config = { 'IBRION': 0, 'NSW': 0, 'POTIM': 0., 'TEIN': 0., 'TEBEG': 0., 'TEEND': 0., 'SMASS': 0. }
        config_patterns = {}
        keys_found = {}
        for key in config.keys():
            keys_found[key] = False
            config_patterns[key] = re.compile(key+'[ \t]*=[ \t]*([0-9.\-]+)')

        for line in self.file:
            allkeys_found = True
            for key in config.keys():
                m = config_patterns[key].search(line)
                if m:
                    config[key] = float(m.group(1))
                    keys_found[key] = True
                if not keys_found[key]:
                    allkeys_found = False
            if allkeys_found:
                break
        if not allkeys_found:
            print "WARNING! Not all config keys were found! Perhaps the OUTCAR format has changed?"
            print "Keys not found:"
            for key in config.keys():
                if not keys_found[key]:
                    print key
        self.config = config 

        #print (self.file.lineno+1),"lines read"
        self.file.reset()

    def get_ionic_steps(self):
        """
        returns a dictionary with entries for time and energy.
        
        >>> d = get_ionic_steps()
        >>> plt.plot(d['time'],d['energy']['total'])
        
        This function could be optimized, since all the lines we are interested
        in for each step occurs after each other, we don't really have to search
        for each line
        """
        self.file.reset()
        numsteps = self.config['NSW']
        print " extracting ionic step data for %.0f steps..." % (numsteps)
        a = {
            'time': np.arange(1,numsteps+1),
            'energy': {
                'ion_electron': np.zeros((numsteps)),
                'ion_kinetic': np.zeros((numsteps)),
                'total': np.zeros((numsteps))
            },
            'forces': {
                'max_atom': np.zeros((numsteps)),
                'rms': np.zeros((numsteps))
            }
        }
        if self.config['IBRION'] == 0.:  # If MD
            a['time'] *= self.config['POTIM']

        for line in self.file:
            
            m = re.search('Iteration[ \t]*([0-9]+)', line)
            if m: 
                stepno = int(m.group(1))
                
            m = re.search('FORCES: max atom, RMS[ \t]*([0-9.\-]+)[ \t]*([0-9.\-]+)', line)
            if m: 
                a['forces']['max_atom'][stepno-1] = float(m.group(1))
                a['forces']['rms'][stepno-1] = float(m.group(2))

            m = re.search('% ion-electron   TOTEN[ \t]*=[ \t]*([0-9.\-]+)', line)
            if m: 
                a['energy']['ion_electron'][stepno-1] = float(m.group(1))

            m = re.search('kinetic Energy EKIN[ \t]*=[ \t]*([0-9.\-]+)', line)
            if m: 
                a['energy']['ion_kinetic'][stepno-1] = float(m.group(1))

            m = re.search('maximum distance moved by ions[ \t]*:[ \t]*([0-9.\-E]+)', line)
            if m: 
                a['energy']['total'][stepno-1] = float(m.group(1))
        return a


    def readItAll(self):
        outfile = open(self.filename)
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
                    if self.selective:
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

    def get_incar_property(self, propname):
        outfile = open(self.filename, 'r')
        lines = outfile.readlines()
        s = re.compile('[\t ]*'+propname+'[\t ]*=[\t ]*([0-9.]*)')
        for l in lines:
            res = s.match(l)
            if res:
                return res.group(1)

        print "Failed to lookup INCAR property "+propname+" in "+self.filename
        sys.exit(1)
    
    def get_num_kpoints(self):
        return self.kpoints
    #def read_stress(self):
    #    for line in open('OUTCAR'):
    #        if line.find(' in kB  ') != -1:
    #            stress = -np.array([float(a) for a in line.split()[2:]]) \
    #                     [[0, 1, 2, 4, 5, 3]] \
    #                     * 1e-1 * ase.units.GPa
    #    return stress
        

class PoscarParser:
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


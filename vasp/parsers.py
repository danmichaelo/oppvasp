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

from copy import copy
import sys,re,math,os
from StringIO import StringIO
import numpy as np

import oppvasp
from oppvasp import get_atomic_number_from_symbol
from oppvasp.md import Trajectory
from oppvasp.structure import Structure

# Optional:
imported = { 'progressbar' : False, 'psutil' : False, 'lxml' : False }

try:
    from lxml import etree
    imported['lxml'] = True
    # Useful reading: 
    # - http://codespeak.net/lxml/parsing.html
    # - http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
except ImportError:
    print "Info: Module 'lxml' is not available"

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
        RotatingMarker, ReverseBar, SimpleProgress
    imported['progressbar'] = True
except ImportError:
    print "Info: Module 'progressbar' is not available"

try:
    import psutil 
    imported['psutil'] = True
except ImportError:
    print "Info: Module 'psutil' is not available"


def read_trajectory(dir = '', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', 
        npz_file = 'trajectory_nopbc.npz', poscar_file = '', unwrap_pbcs = False):
    """
    This function will read a trajectory from the vasprun.xml file given as xml_file and save 
    the trajectory as a NumPy npz cache file npz_pbc_file. The function will then unwrap the 
    periodic boundary conditions and save the trajectory for the unwrapped trajectory as a 
    new NumPy npz cache file npz_file.

    The function will check if a npz file exists, and if it does, read the npz cache instead of the 
    vasprun.xml file. This is *much* faster.

    If POSCAR_file is specified, the initial positions will be read from this file. This may be 
    useful for some visualization purposes. While coordinates in POSCAR (and CONTCAR) may be 
    negative, the coordinates in vasprun.xml and XDATCAR are always wrapped around to become positive.

    Parameters
    ----------
    dir : str
        Directory to read files from
        Default is current directory
    xml_file : str
        Filename of vasprun.xml file
        Default is 'vasprun.xml'
    npz_pbc_file : str
        Filename of the npz cache file to write. 
        Set to blank to disable writing of this file.
        Default is 'trajectory_pbc.npz'
    npz_file : str
        Filename of the npz cache file to write with periodic boundary conditions (PBCs) unwrapped. 
        Set to blank to disable the writing of this file.
        Default is 'trajectory_nopbc.npz'
    poscar_file : str
        Filename of poscar file to get ideal lattice sites from. If given, this will be used upon 
        unwrapping the periodic boundary conditions (PBCs).
        Default is '' 
    unwrap_pbcs : bool
        Set to True to unwrap the periodic boundary conditions (PBCs) or False to leave them.
        Default is False 
    """

    if unwrap_pbcs and os.path.isfile(dir + npz_file):
        traj = Trajectory(filename = dir +npz_file)
    elif not unwrap_pbcs and os.path.isfile(dir + npz_pbc_file):
        traj = Trajectory(filename = dir +npz_pbc_file)
    else:
        p = IterativeVasprunParser(dir + xml_file)
        traj = p.get_all_trajectories()
        if  npz_pbc_file != '':
            traj.save(dir + npz_pbc_file)
        if unwrap_pbcs:
            if poscar_file != '':
                poscar = PoscarParser(dir + poscar_file)
                pos = poscar.get_positions( coords = 'direct' )
                print "Unwrapping using given initial pos"
                traj.unwrap_pbc( init_pos = pos)
            else:
                traj.unwrap_pbc()
            if npz_file != '':
                traj.save(dir + npz_file)
    return traj


def ZapControlCharacters(filename):
    print "Document contains control characters that break the parser!"
    print "Trying to zap them... (this may take some time)"
    doc = open(filename).read()
    # zap control characters that invalidates the xml
    doc = re.sub('[\x00-\x09\x0B-\x1F]','',doc)
    f = file(filename,'w')
    f.write(doc)
    f.close()


class myFile(object):
    def __init__(self, filename):
        self.f = open(filename)

    def read(self, size=None):
        # zap control characters that invalidates the xml
        #return self.f.next().replace('\x1e', '').replace('some other bad character...' ,'')
        return re.sub('[\x00-\x09\x0B-\x1F]','',self.f.next())



def print_memory_usage():
    if imported['psutil']:
        p = psutil.Process(os.getpid())
        rss,vms = p.get_memory_info()
        print "[current physical memory usage: %.1f MB]" % (rss/1024.**2)

class IterativeVasprunParser(object):
    """
    Parser for very large vasprun.xml files, based on iterative xml parsing.
    The functionality of this parser is limited compared to VasprunParser.
    """

    def __str__(self):
        return "iterative vasprun parser"
    
    def __init__(self, filename = 'vasprun.xml', verbose = False):
        
        if not imported['lxml']:
            print "Error: The module 'lxml' is needed!"
            sys.exit(1)
        
        self.filename = filename
        self.verbose = verbose

        if not os.path.isfile(self.filename):
            print "Fatal error: The file '%s' was not found or is not a file." % (self.filename)
            sys.exit(1)

        print_memory_usage()
        
        # read beginning of file to find number of ionic steps (NSW) and timestep (POTIM)
        self.params = self._find_first_instance('parameters', self._params_tag_found)
        self.nsw = int(self.params.xpath("separator[@name='ionic']/i[@name='NSW']")[0].text)
        if self.nsw == 0:
            print "Note: This file contains no ionic motion (NSW=0)."
            self.nsw = 1 # to read the static structure

        # should make a try clause
        self.potim = float(self.params.xpath("separator[@name='ionic']/i[@name='POTIM']")[0].text)

        self.atoms = self._find_first_instance('atominfo',self._get_atoms)
        self.natoms = len(self.atoms)

        try:
            self.nsw
            #print "Number of ionic steps: %d" % (self.nsw) 
        except AttributeError:
            print "Could not find incar:NSW in vasprun.xml"
            sys.exit(1)

    def _params_tag_found(self, elem):
        return copy(elem)
    
    def _get_atoms(self, elem):
        atoms = []
        for rc in elem.xpath("array[@name='atoms']/set/rc"):
            atoms.append(get_atomic_number_from_symbol(rc[0].text))
        return np.array(atoms, dtype=int)

    def _fast_iter(self, context, func):
        for event, elem in context:
            func(elem)
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
        del context

    def _find_first_instance(self, tag, func):
        parser = etree.XMLParser()
        context = etree.iterparse(self.filename, tag=tag)
        ret = None
        try:
            for event, elem in context:
                ret = func(elem)
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
                break
        except etree.XMLSyntaxError:
            print "XML parsing failed:"
            type, message, traceback = sys.exc_info()
            print message
            for e in parser.error_log:
                print "XML Error: "+e.message
            if str(message).split()[0] == 'Char':
                ZapControlCharacters(self.filename)
                print
                print "You may now try to re-run the script."
                print

            sys.exit(1)
        del context
        return ret

    def get_num_ionic_steps(self):
        """ Returns the number of ionic steps """
        return self.nsw

    def get_num_atoms(self):
        """ Returns the number of atoms """
        return self.natoms

    def get_atoms(self):
        """ Returns an array with the types of the atoms """
        return self.atoms

    def _get_initial_structure(self,elem):
        basis= elem.xpath("crystal/varray[@name='basis']/v")
        basis = [[float(x) for x in p.text.split()] for p in basis]

        pos = elem.xpath("varray[@name='positions']/v")
        pos = [[float(x) for x in p.text.split()] for p in pos]

        vel = elem.xpath("varray[@name='velocities']/v")
        vel = [[float(x) for x in p.text.split()] for p in vel]

        return { 'basis': basis, 'positions': pos, 'velocities': vel }

    def get_initial_structure(self):
        """
        Returns a dictionary containing 'basis', 'positions' and 'velocities' 
        """
        return self._find_first_instance('structure',self._get_initial_structure) 

    def _calculation_tag_found(self, elem):

        bas = elem.xpath("structure/crystal/varray[@name='basis']/v")
        self.trajectory.set_basis(self.step_no, np.array([[float(x) for x in b.text.split()] for b in bas]))

        if self.trajectory.num_atoms == 1:
            pos = elem.xpath("structure/varray[@name='positions']/v[%d]" % (self.atom_no+1))
            forces = elem.xpath("structure/varray[@name='forces']/v[%d]" % (self.atom_no+1))
        else:
            pos = elem.xpath("structure/varray[@name='positions']/v")
            forces = elem.xpath("varray[@name='forces']/v")
        pos = [[float(x) for x in ap.text.split()] for ap in pos]
        forces = [[float(x) for x in ap.text.split()] for ap in forces]
        self.trajectory.set_positions(self.step_no, pos)
        self.trajectory.set_forces(self.step_no, forces)
        
        e_kin = elem.xpath("energy/i[@name='kinetic']")
        if e_kin:
            self.trajectory.set_e_kinetic(self.step_no, float(e_kin[0].text))
        
        e_pot = elem.xpath("energy/i[@name='e_fr_energy']")
        self.trajectory.set_e_total(self.step_no, float(e_pot[0].text))

        self.step_no += 1
        if imported['progressbar']:
            self.pbar.update(self.step_no)
        #print pos

    def _get_trajectories(self):
        atoms = self.get_atoms()
        self.trajectory = Trajectory(num_steps = self.nsw, timestep = self.potim, atoms = atoms)
        self.step_no = 0
        status_text = "Parsing %.2f MB... " % (os.path.getsize(self.filename)/1024.**2)
        if imported['progressbar']:
            self.pbar = ProgressBar(widgets=[status_text,Percentage()], maxval = self.nsw+1).start()
        
        parser = etree.XMLParser()
        context = etree.iterparse(self.filename, tag='calculation')
        try:
            self._fast_iter(context, self._calculation_tag_found)
        except etree.XMLSyntaxError:
            type, message, traceback = sys.exc_info()
            print "XML parsing halted:",message
            for e in parser.error_log:
                print "Warning: "+e.message

        if imported['progressbar']:
            self.pbar.finish()
        print "Found %d out of %d steps" % (self.step_no,self.nsw)
        self.trajectory.update_length(self.step_no)
        print_memory_usage()

    def get_all_trajectories(self):
        """
        get trajectories of all atoms
        """
        self.atom_no = -1
        self._get_trajectories()
        return self.trajectory

    def get_single_trajectory(self, atom_no):
        """
        Returns a single trajectory as (stepno, 3) array
        The index of the first atom is 0.
        """
        self.atom_no = atom_no
        self._get_trajectories()
        self.traj['positions'] = self.traj['atoms'][0]['positions']
        return self.trajectory


class IonicStep(object):
    def __init__(self, node, atoms):
        self._node = node
        self._atoms = atoms
    
    def get_total_energy(self):
        """Returns the total energy in electronvolt"""
        results = self._node.xpath( "energy/i[@name='e_fr_energy']")
        if results:
            return float(results[0].text)
        else:
            raise LookupError('Value not found')
    
    def get_stress(self):
        """Returns stress matrix in kB"""
        results = self._node.xpath( "varray[@name='stress']/v")
        if results:
            return np.array([[float(x) for x in p.text.split()] for p in results])
        else:
            raise LookupError('Value not found')

    def get_pressure(self):
        stress = self.get_stress()
        return (stress[0,0] + stress[1,1] + stress[2,2])/3.0
        
    def get_structure(self):
        struc = self._node.xpath("structure")[0]

        basis = struc.xpath("crystal/varray[@name='basis']/v")
        basis = [[float(x) for x in p.text.split()] for p in basis]
        
        rec_basis = struc.xpath("crystal/varray[@name='rec_basis']/v")
        rec_basis = [[float(x) for x in p.text.split()] for p in rec_basis]

        pos = struc.xpath("varray[@name='positions']/v")
        pos = [[float(x) for x in p.text.split()] for p in pos]

        vel = struc.xpath("varray[@name='velocities']/v")
        if vel:
            vel = [[float(x) for x in p.text.split()] for p in vel]
        else:
            # Found no velocities
            vel = None

        forces = self._node.xpath("varray[@name='forces']/v")
        if forces:
            forces = [[float(x) for x in p.text.split()] for p in forces]
        else:
            # Found no forces
            forces = None

        return Structure( cell = basis, atomtypes = self._atoms, positions = pos, velocities = vel, forces = forces, coords = 'direct' )

class VasprunParser(object):
    """
    Parser for vasprun.xml files, making use of libxml for relatively fast parsing.
    """
    
    def __str__(self):
        return "vasprun parser"
    
    def __init__(self, filename = 'vasprun.xml', verbose = False):
        
        if not imported['lxml']:
            print "Error: The module 'lxml' is needed!"
            sys.exit(1)
        
        #print_memory_usage()
        if verbose:
            print "Reading %s (%.2f MB)... " % (filename,os.path.getsize(filename)/1024.**2)

        self.filename = filename
        docstr = open(filename).read()

        # zap control characters that invalidates the xml
        docstr = re.sub('[\x00-\x09\x0B-\x1F]','',docstr)

        if verbose:
            print "Parsing... " 

        parser = etree.XMLParser() #recovers from bad characters.
        try:
            self.doc = etree.parse(StringIO(docstr), parser)
            #self.doc = etree.parse(self.filename, parser)
            for e in parser.error_log:
                print "Warning: "+e.message
        except etree.XMLSyntaxError:
            print "Failed to parse xml file: ",filename
            for e in parser.error_log:
                print "Warning: "+e.message
            sys.exit(2)
        #print_memory_usage()

    @property
    def ionic_steps(self):

        atoms = []
        for rc in self.doc.xpath("/modeling/atominfo/array[@name='atoms']/set/rc"):
            atoms.append(get_atomic_number_from_symbol(rc[0].text.strip()))
        #atoms = [rc[0].text.strip() for rc in atoms]
        atoms = np.array(atoms, dtype=int)

        results = self.doc.xpath( "/modeling/calculation")
        if results:
            return [IonicStep(r, atoms) for r in results]
        else:
            raise LookupError('Value not found')

    def get_final_energy(self):
        return self.ionic_steps[-1].get_total_energy()
    
    def get_final_pressure(self):
        return self.ionic_steps[-1].get_pressure()

    def get_total_energy(self):
        print "VasprunParser: Method get_total_energy is deprecated! "+\
                "Please use get_final_energy instead"
        return self.get_final_energy()

    def get_cpu_time(self):
        return self.get_time_spent()[0]
    
    def get_clock_time(self):
        return self.get_time_spent()[1]

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
    
    def get_kpoints(self):
        """Returns the number of k-points"""
        results = self.doc.xpath( "/modeling/kpoints/varray[@name='kpointlist']/v")
        if results:
            return np.array([[float(x) for x in k.text.split()] for k in results])
        else:
            raise LookupError('Value not found')    

    def get_num_kpoints(self):
        """Returns the number of k-points"""
        return self.get_kpoints().shape[0]

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

    def get_named_structure(self, name):
        """ Returns a named structure as a oppvasp.Structure object. """

        struc = self.doc.xpath("/modeling/structure[@name='%s']" % (name))[0]

        basis = struc.xpath("crystal/varray[@name='basis']/v")
        basis = [[float(x) for x in p.text.split()] for p in basis]
        
        rec_basis = struc.xpath("crystal/varray[@name='rec_basis']/v")
        rec_basis = [[float(x) for x in p.text.split()] for p in rec_basis]

        vol = struc.xpath("crystal/i[@name='volume']")[0]
        vol = float(vol.text)

        pos = struc.xpath("varray[@name='positions']/v")
        pos = [[float(x) for x in p.text.split()] for p in pos]

        vel = struc.xpath("varray[@name='velocities']/v")
        if vel:
            vel = [[float(x) for x in p.text.split()] for p in vel]
        else:
            # Found no velocities
            vel = None

        forces = struc.xpath("varray[@name='forces']/v")
        if forces:
            forces = [[float(x) for x in p.text.split()] for p in forces]
        else:
            # Found no forces
            forces = None

        atoms = []
        for rc in self.doc.xpath("/modeling/atominfo/array[@name='atoms']/set/rc"):
            atoms.append(get_atomic_number_from_symbol(rc[0].text.strip()))
        #atoms = [rc[0].text.strip() for rc in atoms]
        atoms = np.array(atoms, dtype=int)

        return Structure( cell = basis, atomtypes = atoms, positions = pos, velocities = vel, forces = forces, coords = 'direct' )

    def get_initial_structure(self):
        """ Returns the initial structure as a oppvasp.Structure object. """
        return self.get_named_structure('initialpos')

    def get_final_structure(self):
        """ Returns the final structure as a oppvasp.Structure object. """
        s = self.get_named_structure('finalpos')
        if not s.has_forces():
            forces = self.doc.xpath( "/modeling/calculation/varray[@name='forces']/v")
            tmp = []
            for force in forces:
                tmp.append(([float(f) for f in force.text.split()]))
            s.set_forces(np.array(tmp), 'direct')
        return s
    

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
    
    def get_atom_trajectory(self,atom_no):
        """
        Returns a (n,3) numpy array with the position of atom <atom_no> for n timesteps.
        The index of the first atom is 0.
        """
        steps = self.doc.xpath( "/modeling/calculation" )
        num_steps = len(steps)
        traj = np.zeros((num_steps,3))
        i = 0
        for step in steps:
            pos = step.xpath( "structure/varray[@name='positions']/v[%d]" % (atom_no+1) )[0].text.split()
            traj[i] = [float(p) for p in pos]
            i += 1
        
        print "Found %d steps" % (i)

        return traj 

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
        forces = self.doc.xpath( "/modeling/calculation[last()]/varray[@name='forces']/v")
        max_force = 0.0
        max_force_atm = -1
        for i,f in enumerate(forces):
            force = np.array([float(n) for n in f.text.split()])
            force_norm = np.sqrt(np.dot(force,force))
            if force_norm > max_force:
                max_force = force_norm
                max_force_atm = i
        return (max_force,max_force_atm)

    def get_time_spent(self):
        """
        Returns a tupple (CPU time, real time) of seconds spent. 
        Note that value is slightly lower than the value printed at
        the end of OUTCAR. The value corresponds to the sum of the
        times given on lines like this:
          LOOP+:  cpu time   35.95: real time   36.12
        This sum is usually slightly lower than the values reported
        on the last lines of the OUTCAR file
          Total CPU time used (sec):       37.031

        vasprun.xml sum: 157937, 162259
        OUTCAR: 157943, 162261
        """
        calc = self.doc.xpath( "/modeling/calculation" )
        cputime = 0.0
        realtime = 0.0
        for c in calc:
            stepcpu, stepreal = [float(t) for t in c.xpath("time[@name='totalsc']")[0].text.split()]
            #print stepcpu
            cputime += stepcpu
            realtime += stepreal
        return cputime, realtime



class FileIterator(object):
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


class OutcarParser(object):
    """
    Parser for OUTCAR files
    """
    
    def __str(self):
        return "OUTCAR parser"

    def __init__(self, outcarname = 'OUTCAR', selective_dynamics = 0, verbose = False):
        
        if verbose:
            print "Parsing %s (%.1f MB)... " % (outcarname,os.path.getsize(outcarname)/1024.**2)

        self.filename = outcarname
        self.file = FileIterator(self.filename)
        self.selective_dynamics = selective_dynamics 
        
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
        print " extracting ionic step data for %d steps..." % (numsteps)
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
            
            m = re.search('total energy   ETOTAL[ \t]*=[ \t]*([0-9.\-E]+)', line)
            if m: 
                a['energy']['total'][stepno-1] = float(m.group(1))

            #m = re.search('maximum distance moved by ions[ \t]*:[ \t]*([0-9.\-E]+)', line)
            #if m: 
            #    a['energy']['total'][stepno-1] = float(m.group(1))

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
                    if self.selective_dynamics:
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

    def get_max_drift(self):
        return self.maxdrift

    def get_max_pressure(self):
        return self.maxpressure

    def get_max_force(self):
        return self.maxforce

    def get_total_energy(self):
        return self.toten

    def get_cpu_time(self):
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
        

class PoscarParser(object):
    """
    Parser for POSCAR files
    """
    
    def __init__(self, poscar='POSCAR', silent=False):
        self.selective_dynamics = False
        if isinstance(poscar, file):
            self._file = poscar
            self.filename = poscar.name
        else:
            self.filename = poscar 
            self._file = open( self.filename, 'r')
        self._parse(silent)

    def strip_comments(self, line):
        # Comments start with "!"
        return re.sub('!.*','',line)

    def _parse(self, silent=False):

        poscarfile = self._file
        # Read comment line
        commentline = poscarfile.readline()

        # Read three unit cell lines
        self.scale_factor = float(poscarfile.readline()) # lattice constant
        vec1line = self.strip_comments(poscarfile.readline())
        vec2line = self.strip_comments(poscarfile.readline())
        vec3line = self.strip_comments(poscarfile.readline())
        self.basis = np.zeros((3,3))
        self.basis[0] = map(float,vec1line.split())
        self.basis[1] = map(float,vec2line.split())
        self.basis[2] = map(float,vec3line.split())
        self.basis *= self.scale_factor
        
        # Read the atom counts line
        # and check for optional line containing atom types (VASP >= 5.1)
        self.atomtypes = []
        sixthline = self.strip_comments(poscarfile.readline())
        try:
            dummy = int(sixthline.split()[0])
            atomcountline = sixthline
        except:
            self.atomtypes = np.array([get_atomic_number_from_symbol(i) for i in sixthline.split()])
            atomcountline = self.strip_comments(poscarfile.readline())
        self.atomcounts = map(int, atomcountline.split())
        self.natoms = sum(self.atomcounts)
        if len(self.atomtypes) == 0:
            if not silent:
                print "Note: Atom types are not specified in '%s'" % self.filename
            # assign negastive numbers, so we can still differentiate between the different types
            # (this will probably cause som nasty bugs at some time...)
            self.atomtypes = np.arange(-1,-len(self.atomcounts)-1,-1) 

        # Check for optional line specifying selective dynamics
        self.selective_dynamics = False
        line = self.strip_comments(poscarfile.readline())
        if line[0].upper() == 'S':
            self.selective_dynamics = True
            line = self.strip_comments(poscarfile.readline())
        
        # Read coordinate type line (direct or cartesian)
        self.direct_coords = False
        if line[0].upper() == 'D':
            self.direct_coords = True

        # Atomic positions
        self.positions = np.zeros((self.natoms,3))
        for j in range(self.natoms):
            line = self.strip_comments(poscarfile.readline())
            self.positions[j] = [float(x) for x in (line.split()[0:3])]
            if len(self.positions[j]) != 3:
                raise IOError('POSCAR file %s ended before all coordinates were read in' % self.filename)
        
        poscarfile.close()
        
    def get_atomtypes(self):
        if len(self.atomtypes) != len(self.atomcounts):
            return None
        atomarray = []
        for atcount, attype in zip(self.atomcounts, self.atomtypes):
            for j in range(atcount):
                atomarray.append(attype)
        return atomarray

    
    def get_positions(self, coords = 'direct'):
        if self.direct_coords == True and coords[0].lower() == 'c':
            return oppvasp.direct_to_cartesian(self.positions, self.basis)
        elif self.direct_coords == False and coords[0].lower() == 'd':
            return oppvasp.cartesian_to_direct(self.positions, self.basis)
        else:
            return self.positions.copy()

    def get_basis(self):
        return self.basis

    def get_scale_factor(self):
        """ lattice constant"""
        return self.scale_factor

    def get_structure(self):
        return Structure( cell = self.basis.copy(), positions = self.get_positions('direct'), atomtypes = self.get_atomtypes(), coords = 'direct' )
            


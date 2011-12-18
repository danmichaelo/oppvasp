# encoding=utf-8

import numpy as np
from oppvasp import get_atomic_number_from_symbol, elements, direct_to_cartesian, cartesian_to_direct
from oppvasp.util import get_pairs

def unique_list(seq, idfun=None):
    # Example: unique_list(a, lambda x: x.lower())
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

class Structure(object):
    """
    A Structure object defines a periodic atomic structure, keeping 
    information about its unit cell together with the types, positions 
    and possibly forces of its atoms. Internally, coordinates are
    stored in direct coordinates, to avoid more conversion than necessary

    Properties:

        - cell
        - atomtypes
        - positions
        - forces 
        - velocities
    """
    
    def __init__(self, cell = [], positions = [], atomtypes = None, forces = None, velocities = None, coords = 'direct'):
        self.set_cell(cell)
        self.set_atomtypes(atomtypes)
        self.set_positions(positions, coords)
        self.set_forces(forces, coords)
        self.set_velocities(velocities, coords)
   
    def get_num_atoms(self):
        """ Returns the number of atoms """
        return self._positions.shape[0]

    def validate(self):
        """ 
        Validates the Structure object is internally consistent;
        that arrays holding positions, atomtypes, .. share the
        same dimensions
        """
        if self._positions.shape[0] != self._atomtypes.shape[0]:
            print "WARNING: pos length:",self._positions.shape[0],"Atomtypes length:",self._atomtypes.shape[0]
            return False

        # Fill in more tests here

        return True

    
    #------------------- Properties -------------------------------

    def get_cell(self):
        """ 3x3 unit cell matrix """
        return self._cell

    def set_cell(self, cell):
        """ 
        Define the unit cell using a (3,3) numpy array
        """
        if type(cell).__name__ == 'list':
            self._cell = np.array(cell)
        elif type(cell).__name__ == 'ndarray':
            self._cell = cell.copy()
        else:
            print "Error: cell must be a list or numpy array"

    def get_atomtypes(self):
        """ 
        Array of atom numbers, having the same length as the number of 
        atoms in the structure. 

        Returns: List of strings or None
        """
        #if np.any(self._atomtypes < 0):
        #    print "Warning: Atom type unknown for one or more atoms!"
        return self._atomtypes

    def set_atomtypes(self, atoms):
        """
        If the array is set using a list of atomic symbols,
        it will be converted into an array of atomic numbers, so 
        the following two lines of code are equivalent:

        >>> set_atomtypes([14,16,16])
        >>> set_atomtypes(['Si','O','O'])
        """
        if type(atoms).__name__ == 'ndarray':
            self._atomtypes = np.array(atoms, dtype = 'int')
        elif type(atoms).__name__ == 'list':
            self._atomtypes = np.zeros(len(atoms), dtype = 'int')
            for i, atom in enumerate(atoms):
                if type(atom).__name__ == 'str':
                    self._atomtypes[i] = get_atomic_number_from_symbol(atom)
                else:
                    self._atomtypes[i] = atom
        elif type(atoms).__name__ == 'NoneType':
            self._atomtypes = None
        else:
            print "Error: atoms must be a list or array"
    
    def get_positions(self, coords):
        if coords[0].lower() == 'd':
            return self._positions.copy()
        elif coords[0].lower() == 'c':
            return direct_to_cartesian(self._positions, self._cell)

    def set_positions(self, pos, coords):
        """ 
        Define atom positions using a (n,3) numpy array as the first argument.
        The function assumes cartesian coordinates by default, but
        a second argument can be added to indicate if the coordinates
        are direct ('d') or cartesian ('c'). 
        """
        if type(pos).__name__ == 'ndarray':
            self._positions = pos.copy()
        elif type(pos).__name__ == 'list':
            self._positions = np.array(pos)
        else:
            print "Error: positions must be a list or numpy array"
        if coords[0].lower() == 'c':
            self._positions = cartesian_to_direct(self._positions, self._cell)
    
    
    def get_forces(self, coords):
        """ (n,3) array """
        if coords[0].lower() == 'd':
            return self._forces
        elif coords[0].lower() == 'c':
            return direct_to_cartesian(self._forces, self._cell)

    def set_forces(self, vel, coords):
        if type(vel).__name__ == 'ndarray':
            self._forces = vel.copy()
        elif type(vel).__name__ == 'list':
            self._forces = np.array(vel)
        elif type(vel).__name__ == 'NoneType':
            self._forces = None
        else:
            print "Error: forces must be a list or numpy array"
            self._forces = None
        if self._forces != None and coords[0] == 'c':
            self._forces = cartesian_to_direct(self._forces, self._cell)

    def get_velocities(self, coords):
        """ (n,3) array """
        if coords[0].lower() == 'd':
            return self._velocities
        elif coords[0].lower() == 'c':
            return direct_to_cartesian(self._velocities, self._cell)

    def set_velocities(self, vel, coords):
        if type(vel).__name__ == 'ndarray':
            self._velocities = vel.copy()
        elif type(vel).__name__ == 'list':
            self._velocities = np.array(vel)
        elif type(vel).__name__ == 'NoneType':
            self._velocities = None
        else:
            print "Error: velocities must be a list or numpy array"
            self._velocities = None
        if self._velocities != None and coords[0] == 'c':
            self._velocities = cartesian_to_direct(self._velocities, self._cell)

    #------------------- Methods -------------------------------
    
    def save(self, filename, format = 'POSCAR', direct_coords = True):
        """ Save as POSCAR """

        if not self.validate():
            print "ERROR: The current structure failed internal self-validation"
            return

        if isinstance(filename, file):
            f = filename
        else:
            f = open(filename, 'w')
        f.write('Saved ---\n')
        f.write('1.0\n')
        for cellvec in self.get_cell():
            f.write( "   % .16f   % .16f   % .16f\n" % tuple(cellvec))
        atomline1 = ""
        atomline2 = ""
        
        atom_numbers = self.get_atomtypes()
        atms_list = unique_list(atom_numbers)
        atms_mult = [(at == atom_numbers).sum() for at in atms_list]
        
        typesknown = True
        for mult,atomtyp in zip(atms_mult,atms_list):
            if atomtyp < 0:
                typesknown = False
            else:
                atomline1 += "  %s" % (elements[atomtyp]['symb']) 
            atomline2 += "  %d" % (mult) 
        if typesknown:
            f.write(atomline1+'\n')
        f.write(atomline2+'\n')
        f.write('Direct\n' if direct_coords else 'Cartesian\n')
        pos = self.get_positions('direct' if direct_coords else 'cart')
        for at in atms_list:
            for atomvec in pos[(at == atom_numbers)]:
                f.write( "   % .16f   % .16f   % .16f\n" % tuple(atomvec))
        f.close()

    def get_ase_atoms_object(self):
        """ Returns an ase.Atoms object """
        from ase import Atoms
        return Atoms(numbers = self.get_atomtypes(), positions = self.get_positions('cart'), cell = self.get_cell())

    def set_from_ase_atoms_object(self, atoms_obj):
        """ Imports an ase.Atoms object """
        from ase import Atoms
        self.set_cell(atoms_obj.get_cell())
        self.set_atomtypes(atoms_obj.get_atomic_numbers())
        self.set_positions(atoms_obj.get_positions())
        self.set_forces(atoms_obj.get_forces())
        self.set_velocities(atoms_obj.get_velocities())

    #-------------- Methods for analyzing the structure -----------

    def get_supercell_positions(self, sx, sy, sz, coords= 'cart'):
        """ Returns a (sx,sy,sz) supercell in the form of a (sx*sy*sz,3) numpy array """
        nat = self.get_num_atoms()
        sup = np.zeros((nat*sx*sy*sz,3))
        c = np.diag((1,1,1)) # since we use direct coordinates
        for i in range(sx):
            for j in range(sy):
                for k in range(sz):
                    m = i*sy*sz + j*sz + k
                    sup[m*nat:(m+1)*nat] = self._positions + np.dot(c,[i, j, k])
        if coords[0].lower() == 'c':
            return direct_to_cartesian(sup, self._cell)
        else:
            return sup

    def get_supercell_structure(self,sx,sy,sz):
        """ Returns a (sx,sy,sz) supercell in the form of a Structure object """
        nat = self.get_num_atoms()
        c = self._cell * [sx,sy,sz]
        p = self.get_supercell_positions(sx,sy,sz,coords='cart')
        t = np.zeros(p.shape[0])
        t0 = self.get_atomtypes()
        for i in range(sx*sy*sz):
            print i*nat,i*nat+nat
            t[i*nat:i*nat+nat] = t0
        return Structure( cell = c, positions = p, atomtypes = t, coords='cart')

    def get_nearest_neighbours(self, atom_no, tolerance = 0.1):
        """ Returns a list of the nearest neighbours to atom atom_no """
        nat = self.get_num_atoms()
        
        # Make a 3x3x3 cell to ensure neighbours across 
        # the unit cell boundary are counted as well.
        sup = self.get_supercell_positions(3,3,3)
        
        # Positions of center cell:
        cent = self._positions[atom_no] + np.dot(self.get_cell(),[1,1,1])

        # Calculate all distances:
        dist = np.sqrt(np.sum((sup-cent)**2,axis=1))

        order = np.argsort(dist)
        nearest_neighbours = [order[1]] 
        d0 = dist[order[1]]
        for n in order[2:]:
            if (dist[n] - d0) > tolerance:
                break
            nearest_neighbours.append(n)

        for n in nearest_neighbours:
            # use mod to get atomic number in unit cell, not supercell
            print n%nat,dist[n]

    def get_volume(self):
        """ Returns the volume in Angstrom^3 """
        return np.dot(self._cell[0],np.cross(self._cell[1],self._cell[2]))

    def wrap_pbc(self):
        """
        Wraps coordinates outside the unit cell into it
        """
        r = self.get_positions('direct')
        # Minimum image convention (may not be appropriate for extreme cell geometries!!)
        rc = r - (2*r-1).astype(np.int)
        self.set_positions(rc, 'direct')
    

    def unwrap_pbc(self, init_pos = None, silent = False):
        """ 
        Unwraps periodic boundary conditions (PBCs).
        Note that this procedure produces coordinates outside the unit cell. 
        
        An initial set of coordinates must be specified as a Structure object, 
        typically originating from a POSCAR file. 

        Example:
        ----------
        
        >>> pp = PoscarParser('POSCAR')
        >>> init_pos = pp.get_structure()
        >>> vp = VasprunParser('vasprun.xml')
        >>> fin_pos = vp.ionic_steps[-1].get_structure()
        >>> fin_pos = fin_pos.unwrap_pbc(init_pos)
        >>> print fin_pos - init_pos
        
        """
        #if self.shape[0] != r0.shape[0]:
        #    print "ERROR: Mismatch between number of atoms in initial and final structure! Aborting"
        #    return

        # Convert initial coordinates to final basis:
        r0 = cartesian_to_direct(init_pos.get_positions('cart'), self._cell)
        r = self.get_positions('direct')
        
        n = r.shape[0]
        dr = r-r0
        c = np.abs(np.array((dr-1, dr, dr+1)))
        diff = np.array([[np.argmin(c[:,i,j]) for j in range(3)] for i in range(n)]) - 1  # (-1,0,1)
        dr += diff

        if not silent and np.max(dr) > 0.1:
            print u"Warning: largest displacement of %.1f (direct units) is rather large..." % (np.max(dr))

        self.set_positions(r0+dr, 'direct')
    
    def get_shortest_bond(self):
        """Returns a tupple containing the shortest bond length in Angstrom,
        followed by the indices of the two atoms making up the bond"""

        pos = self._positions
        natoms = pos.shape[0]

        # Build array of pair indices:
        pairs = get_pairs(natoms)
        npairs = pairs.shape[0]

        # Find displacement vectors for all atoms:
        x = pos[pairs[:,0]] - pos[pairs[:,1]]

        # Use minimum image convention to threat bonds over PBCs
        # Note: This will not work with *very* tilted unit cells
        x = x - (2*x).astype('int')

        X = direct_to_cartesian(x, self._cell)

        r2 = (X**2).sum(axis=1)
        r = np.sqrt(r2)

        meanr2 = np.mean(r2,axis=0)
        meanr = np.mean(r,axis=0)

        minidx = np.argmin(r)
        minpair = pairs[minidx]
        minr = r[minidx]
    
        return (minr,minpair[0],minpair[1])


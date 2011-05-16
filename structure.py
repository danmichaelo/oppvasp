import numpy as np
from oppvasp import get_atomic_number_from_symbol, direct_to_cartesian

class Structure(object):
    """
    A Structure object defines a periodic atomic structure, keeping 
    information about its unit cell together with the types, positions 
    and possibly velocities of its atoms.

    Properties:

        - cell
        - atomtypes
        - positions
        - velocities
    """
    
    def __init__(self, cell = [], atomtypes = [], positions = [], velocities = []):
        self.set_cell(cell)
        self.set_atomtypes(atomtypes)
        self.set_positions(positions)
        self.set_velocities(velocities)
   
    #------------------- Properties -------------------------------

    def get_cell(self):
        """ 3x3 unit cell matrix """
        return self._cell

    def set_cell(self, cell):
        """ Define the unit cell using a (3,3) numpy array. """
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
        """
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
            for i in range(len(atoms)):
                if type(atoms[i]).__name__ == 'str':
                    self._atomtypes[i] = get_atomic_number_from_symbol(atoms[i])
                else:
                    self._atomtypes[i] = atoms[i]
        else:
            print "Error: atoms must be a list or array"
    
    def set_positions(self, pos):
        """ Define atom positions using a (n,3) numpy array. """
        if type(pos).__name__ == 'ndarray':
            self._positions = pos.copy()
        elif type(pos).__name__ == 'list':
            self._positions = np.array(pos)
        else:
            print "Error: positions must be a list or numpy array"
    
    def get_positions(self, coordinates = 'direct'):
        if coordinates == 'direct':
            return self._positions
        elif coordinates == 'cartesian':
            return direct_to_cartesian(self._positions)
    
    def get_velocities(self):
        """ Atom velocities array """
        return self._velocities

    def set_velocities(self, vel):
        if type(vel).__name__ == 'ndarray':
            self._velocities = vel.copy()
        elif type(vel).__name__ == 'list':
            self._velocities = np.array(vel)
        else:
            print "Error: velocities must be a list or numpy array"

    #------------------- Methods -------------------------------

    def get_ase_atoms_object(self):
        """ Returns an ase.Atoms object """
        from ase import Atoms
        return Atoms(numbers = self._atomtypes, positions = self._positions, cell = self._cell)

    def set_from_ase_atoms_object(self, atoms_obj):
        """ Imports an ase.Atoms object """
        from ase import Atoms
        self.set_cell(atoms_obj.get_cell())
        self.set_atomtypes(atoms_obj.get_atomic_numbers())
        self.set_positions(atoms_obj.get_positions())
        self.set_velocities(atoms_obj.get_velocities())


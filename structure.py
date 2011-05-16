import numpy as np
from oppvasp import get_atomic_number_from_symbol

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
        self.cell = cell
        self.atomtypes = atomtypes
        self.positions = positions 
        self.velocities = velocities 
   
    #------------------- Properties -------------------------------

    @property
    def cell(self):
        """ 3x3 unit cell matrix """
        return self._cell

    @cell.setter
    def cell(self, cell):
        """ Define the unit cell using a (3,3) numpy array. """
        if type(cell).__name__ == 'list':
            self._cell = np.array(cell)
        elif type(cell).__name__ == 'ndarray':
            self._cell = cell.copy()
        else:
            print "Error: cell must be a list or numpy array"
    
    @property
    def atomtypes(self):
        """ 
        Array of atom numbers, having the same length as the number of 
        atoms in the structure. If the array is set using a list of atomic symbols,
        it will be converted into an array of atomic numbers, so 
        the following two lines of code are equivalent:

        >>> structure.atomtypes = [14,16,16]
        >>> structure.atomtypes = ['Si','O','O']
        """
        return self._types

    @atomtypes.setter
    def atomtypes(self, atoms):
        if type(atoms).__name__ == 'ndarray':
            self._types = np.array(atoms, dtype = 'int')
        elif type(atoms).__name__ == 'list':
            self._types = np.zeros(len(atoms), dtype = 'int')
            for i in range(len(atoms)):
                if type(atoms[i]).__name__ == 'str':
                    self._types[i] = get_atomic_number_from_symbol(atoms[i])
                else:
                    self._types[i] = atoms[i]
        else:
            print "Error: atoms must be a list or array"
    
    @property
    def positions(self):
        """ Atom positions array """
        return self._positions

    @positions.setter
    def positions(self, pos):
        """ Define atom positions using a (n,3) numpy array. """
        if type(pos).__name__ == 'ndarray':
            self._positions = pos.copy()
        elif type(pos).__name__ == 'list':
            self._positions = np.array(pos)
        else:
            print "Error: positions must be a list or numpy array"
    
    @property
    def velocities(self):
        """ Atom velocities array """
        return self._velocities

    @velocities.setter
    def velocities(self, vel):
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
        return Atoms(numbers = self.atomtypes, positions = self.positions, cell = self.cell)

    def set_from_ase_atoms_object(self, atoms_obj):
        """ Imports an ase.Atoms object """
        from ase import Atoms
        self.atomtypes = atoms_obj.get_atomic_numbers()
        self.positions = atoms_obj.get_positions()
        self.velocities = atoms_obj.get_velocities()
        self.cell = atoms_obj.get_cell()


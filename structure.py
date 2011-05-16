import numpy as np

class Structure(object):
    
    def __init__(self, cell = [], atom_types = [], atom_pos = [], atom_vel = []):
        self._cell = []
        self._positions = []
        self._velocities = []
        self._types = []
        self.set_cell(cell)
        self.set_atom_pos(atom_pos)
        self.set_atom_vel(atom_vel)
        self.set_atom_types(atom_types)

    def set_cell(self, cell):
        if type(cell).__name__ == 'list':
            self._cell = np.array(cell)
        elif type(cell).__name__ == 'ndarray':
            self._cell = cell 
        else:
            print "Error: cell must be a list or numpy array"

    def set_atom_pos(self, positions):
        if type(positions).__name__ == 'list':
            self._positions = np.array(positions)
        elif type(positions).__name__ == 'ndarray':
            self._positions = positions
        else:
            print "Error: positions must be a list or numpy array"

    def set_atom_vel(self, velocities):
        if type(velocities).__name__ == 'list':
            self._velocities = np.array(velocities)
        elif type(velocities).__name__ == 'ndarray':
            self._velocities = velocities
        else:
            print "Error: velocities must be a list or numpy array"

    def set_atom_types(self, types):
        if type(types).__name__ == 'list':
            self._types = types
        else:
            print "Error: types must be a list"

    def get_ase_atoms_object(self):
        from ase import Atoms
        types = ''.join(self._types)
        return Atoms(types, positions = self._positions, cell = self._cell)


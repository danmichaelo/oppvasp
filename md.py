
import numpy as np

class Trajectory:

    def __init__(self, num_steps = 0, num_atoms = 0, filename = '', timestep = 1., atoms = []):
        if filename != '':
            self.load(filename)
        else:
            self.atoms = atoms
            self.length = num_steps
            self.num_atoms = num_atoms
            self.time = np.zeros(self.length)
            self.basis = [ np.zeros((3,3)) for i in range(self.length) ]
            self.total_energy = np.zeros(self.length)
            self.kinetic_energy = np.zeros(self.length)
            if self.num_atoms > 0:
                self.positions = np.zeros((num_steps, num_atoms, 3))
            self.set_timestep(timestep)

    def update_length(self, new_length):
        if new_length < self.length:
            print "Stripping %d empty steps" % (self.length - new_length)
            self.length = new_length 
            self.total_energy = self.total_energy[0:new_length]
            self.kinetic_energy = self.kinetic_energy[0:new_length]
            self.basis = self.basis[0:new_length]
            
            #self.traj['atoms'][i]['positions'] = self.traj['atoms'][i]['positions'][0:newlen]

    def set_timestep(self, timestep):
        self.time = np.arange(self.length) * timestep
    
    def set_basis(self, step_no, bas):
        self.basis[step_no] = bas

    def set_positions(self, step_no, pos):
        if self.num_atoms == 0:
            self.num_atoms = len(pos)
            self.positions = np.zeros((self.length, self.num_atoms, 3))
        a=np.array(pos)
        self.positions[step_no] = pos

    def get_single_particle_trajectory(self, atom_no):
        return self.positions[:,atom_no,:]

    def set_e_total(self, step_no, etot):
        self.total_energy[step_no] = etot 

    def set_e_kinetic(self, step_no, ekin):
        self.kinetic_energy[step_no] = ekin

    def save(self, filename):
        """ Saves the trajectory to a *.npz binary file using numpy.savez """
        np.savez(filename, basis = self.basis, tote = self.total_energy, kine = self.kinetic_energy, pos = self.positions, 
            time = self.time, atoms = self.atoms)
        print "Wrote trajectory (%d atoms, %d steps) to %s" % (self.num_atoms, self.length, filename)

    def load(self, filename):
        """ Loading the trajectory from a *.npz numpy binary file using numpy.load """
        npz = np.load(filename)
        self.atoms = npz['atoms']
        self.basis = npz['basis']
        self.total_energy = npz['tote']
        self.kinetic_energy = npz['kine']
        self.positions = npz['pos']
        self.time = npz['time']
        self.length = self.positions.shape[0]
        self.num_atoms = self.positions.shape[1] 
        print "Read trajectory (%d atoms, %d steps) from %s" % (self.num_atoms, self.length, filename)



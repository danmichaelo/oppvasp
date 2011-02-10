# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import sys
import numpy as np
import scipy as sp

import time
import oppvasp

# Optional:
imported = { 'progressbar' : False, 'psutil' : False, 'lxml' : False }
try:
    from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
        RotatingMarker, ReverseBar, SimpleProgress
    imported['progressbar'] = True
except ImportError:
    print "Info: Module 'progressbar' is not available"


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
            self.time = self.time[0:new_length]
            self.positions = self.positions[0:new_length]

    def set_timestep(self, timestep):
        self.time = np.arange(self.length) * timestep
    
    def set_basis(self, step_no, bas):
        self.basis[step_no] = bas

    def set_positions(self, step_no, pos):
        if self.num_atoms == 0:
            self.num_atoms = len(pos)
            self.positions = np.zeros((self.length, self.num_atoms, 3))
        self.positions[step_no] = pos

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


    def unwrap_pbc(self, remove_drift = True, init_pos = None):
        """ 
        Unwraps periodic boundary conditions (PBCs)
        Note that this procedure produces coordinates outside the unit cell. 
        
        For each step, the distances between the previous position of an atom 
        and the current position of the atom and its periodic images is 
        calculated. The shortest distance is chosen, and added to the position 
        vector of that atom.

        This function was written for the purpose of making it easier to track 
        the path of a diffusing impurity atom, but may also be useful for 
        visualisation purposes.

        Optionally, an initial set of coordinates may be specified, typically 
        from a POSCAR file. This can be useful if comparing coordinates with
        the ones found in the POSCAR file (which may be negative).
        They should be in direct coordinates, with the same basis as the first 
        step of the trajectory.

        Example:
        
        >>> pp = PoscarParser('POSCAR')
        >>> init_pos = pp.get_positions()
        >>> vp = IterativeVasprunParser('vasprun.xml')
        >>> traj = vp.get_all_trajectories()
        >>> traj.save('trajectory.npz')
        >>> traj.unwrap_pbc( init_pos = init_pos )
        >>> traj.save('trajectory_unwrapped.npz')
        
        """
        if imported['progressbar']:
            pbar = ProgressBar(widgets=['Unwrapping PBCs...',Percentage()], maxval = self.length).start()

        self.dr = np.zeros((self.length, self.num_atoms, 3), dtype='float64')
        self.r = np.zeros((self.length, self.num_atoms, 3), dtype='float64')
        if init_pos == None:
            self.r[0] = self.positions[0]
        else:
            self.r[0] = init_pos
        for stepno in range(1,self.length):
            dr = self.positions[stepno] - self.positions[stepno-1]
            c = np.abs(np.array((dr-1, dr, dr+1)))
            diff = np.array([[np.argmin(c[:,i,j]) for j in range(3)] for i in range(self.num_atoms)]) - 1  # (-1,0,1)
            dr += diff
            self.dr[stepno] = dr
            self.r[stepno] = self.r[stepno-1] + dr

            # This may not be necessary...
            drp = self.r[stepno] - self.positions[stepno]
            drift = drp - np.round(drp)
            if remove_drift:
                self.r[stepno] -= drift
            if abs(np.max(drift)) > 0.5:
                print "Warning: Max drift is",np.max(drift),". This is very high, and may lead to erroneous behaviour of the PBC unwrapping."

            if imported['progressbar']:
                pbar.update(stepno)
        self.positions = self.r
        if imported['progressbar']:
            pbar.finish()

    
    def get_single_particle_trajectory(self, atom_no, coordinates = 'direct'):
        """
        Returns the trajectory of a single atom in direct or cartesian coordinates
        as a numpy array of dimension (steps, 3).
        Parameters:
            - coordinates : either 'direct' or 'cartesian'
        """
        pos = self.positions[:,atom_no,:]
        if coordinates == 'cartesian':
            pos = oppvasp.direct_to_cartesian(pos, self.basis)
        return pos


    def get_all_trajectories(self, coordinates = 'direct'):
        """
        Returns the trajectories for all the atoms in direct or cartesian coordinates
        as a numpy array of dimension (steps, atoms, 3).
        Parameters:
            - coordinates : either 'direct' or 'cartesian'
        """
        if coordinates == 'direct':
            return self.positions
        elif coordinates == 'cartesian':
            if 'positions_cart' in dir(self):
                return self.positions_cart
            else:
                pos = oppvasp.direct_to_cartesian(self.positions, self.basis)
                # This conversion is currently stupidly slow, so we cache the cartesian coordinates
                # as a temporary solution until a faster algorithm is implemented
                self.positions_cart = pos
                return pos


    def get_displacements(self, coordinates = 'direct'):
        """
        Returns the displacements $r(t)-r(0)$ as a numpy array:
            axis[0] = time, 
            axis[1] = atom_no,
            axis[2] = coordinate no (x,y,z)
        """
        pos = self.get_all_trajectories(coordinates)
        return pos - pos[0]


    def get_displacements_squared(self, coordinates = 'direct'):
        """
        Returns the displacements squared $[r(t)-r(0)]^2), as an numpy array:
            axis[0] = time, 
            axis[1] = atom_no,
        $ \Delta r_i^2(t) = (r_i(t)-r_i(0))^2 $
        Parameters:
         - coordinates : either 'direct' or 'cartesian'
        """
        disp = self.get_displacements(coordinates)
        dispsq  = disp[:,:,0]**2 + disp[:,:,1]**2 + disp[:,:,2]**2
        return dispsq
        # transpose to get axis[0] = at no., axis[1] = time

    def get_velocities(self, coordinates = 'direct', algorithm = 'naive'):
        """
        Returns velocities
        Parameters:
         - coordinates: 'direct' returns velocities in [direct coordinates]/[timestep]
                        'cartesian' returns velocities in Å/s
        """
        timestep = self.time[1]-self.time[0]
        pos = self.get_all_trajectories(coordinates)
        vel = np.zeros(pos.shape)
        if algorithm == 'naive':
            for i in range(1,vel.shape[0]-1):
                vel[i] = (pos[i+1] - pos[i-1]) / 2.*timestep

        return vel

def pair_correlation_function(x,y,z,S,rMax,dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple 
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is 
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;) 

    From: http://www.shocksolution.com/microfluidics-and-biotechnology/calculating-the-pair-correlation-function-in-python/
    
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_x, interior_y, interior_z)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        interior_x      x coordinates of reference particles
        interior_y      y coordinates of reference particles
        interior_z      z coordinates of reference particles
    """
    from numpy import zeros, sqrt, where, pi, average, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x>rMax
    bools2 = x<(S-rMax)
    bools3 = y>rMax
    bools4 = y<(S-rMax)
    bools5 = z>rMax
    bools6 = z<(S-rMax)

    interior_indices, = where(bools1*bools2*bools3*bools4*bools5*bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax+1.1*dr, dr)
    num_increments = len(edges)-1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x)/S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index]-x)**2 + (y[index]-y)**2 + (z[index]-z)**2)
        d[index] = 2*rMax

        result, bins = histogram(d, bins=edges, normed=False)

        g[p,:] = np.array(result, dtype='float')/numberDensity
    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1])/2.        
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = average(g[:,i])/(4./3.*pi*(rOuter**3 - rInner**3))

    return (g_average, radii, x[interior_indices], y[interior_indices], z[interior_indices])
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)


# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import sys
import numpy as np
#import scipy as sp

import time
import oppvasp

# Status of optional Python modules:
imported = { 'progressbar' : False, 'psutil' : False, 'lxml' : False }
try:
    from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
        RotatingMarker, ReverseBar, SimpleProgress
    imported['progressbar'] = True
except ImportError:
    print "Info: Module 'progressbar' is not available"


class Trajectory(object):

    def __init__(self, num_steps = 0, filename = '', timestep = 1., atoms = []):
        if filename != '':
            self.load(filename)
        else:
            self.atoms = atoms
            self.length = num_steps
            self.num_atoms = atoms.shape[0]
            self._time = np.zeros(self.length)
            self.basis = [ np.zeros((3,3)) for i in range(self.length) ]
            self.total_energy = np.zeros(self.length)
            self.kinetic_energy = np.zeros(self.length)
            if self.num_atoms > 0:
                self._positions = np.zeros((num_steps, self.num_atoms, 3))
                self._forces = np.zeros((num_steps, self.num_atoms, 3))
            self.set_timestep(timestep)
            print "Allocated %.2f MB" % ((self._time.nbytes + self.total_energy.nbytes + self.kinetic_energy.nbytes + self._positions.nbytes + self._forces.nbytes)/(1024.**2))
        # Set default selection to the whole trajectory:
        self.in_point = 0
        self.out_point = self.length

    def update_length(self, new_length):
        if new_length < self.length:
            print "Stripping %d empty steps" % (self.length - new_length)
            self.length = new_length 
            self.total_energy = self.total_energy[0:new_length]
            self.kinetic_energy = self.kinetic_energy[0:new_length]
            self.basis = self.basis[0:new_length]
            self._time = self.time[0:new_length]
            self._positions = self._positions[0:new_length]
            self._forces = self._forces[0:new_length]
    
    def set_selection(self, in_point = 0, out_point = -1):
        """
        Selects a part of the trajectory. Once a part of the trajectory has been selected, 
        only that part of the trajectory is considered by the various functions in this class.
        To clear the selection, just call set_selection with no arguments.

        Parameters
        ----------
        in_point : int
            The frame number to mark the beginning of the selection
            Default is 0
        out_point : int
            The frame number to mark the end of the selection
            Default is -1

        Examples
        ----------
        >>> traj = Trajectory( filename = 'trajectory.npz' )
        >>> traj.set_selection(10000, 20000)
        >>> d = traj.get_displacements()
        >>> t = traj.time   # getter function returns the correct interval
        >>> plt.plot(t,d)

        """
        self.in_point = in_point
        if out_point == -1:
            out_point = self.length # to include the very last element
        self.out_point = out_point
    
    #------------------- Begin Properties -------------------------------

    @property
    def time(self):
        """
        nsteps array"
        """
        return self._time[self.in_point:self.out_point]

    def set_timestep(self, timestep):
        self._time = np.arange(self.length) * timestep

    @property
    def positions(self):
        """
        nsteps x natoms x 3 array containing the direct coordinates of all atoms for all steps
        """
        return self._positions[self.in_point:self.out_point]

    @positions.setter
    def positions(self, val):
        self._positions = val

    def set_positions(self, step_no, pos):
        if self.num_atoms == 0:
            self.num_atoms = len(pos)
            self._positions = np.zeros((self.length, self.num_atoms, 3))
            self._forces = np.zeros((self.length, self.num_atoms, 3))
        self._positions[step_no] = pos

    @property
    def forces(self):
        """
        nsteps x natoms x 3 array containing the forces on all atoms for all steps
        """
        return self._forces[self.in_point:self.out_point]

    @forces.setter
    def forces(self, val):
        self._forces = val

    def set_forces(self, step_no, forces):
        if self.num_atoms == 0:
            self.num_atoms = len(pos)
            self._positions = np.zeros((self.length, self.num_atoms, 3))
            self._forces = np.zeros((self.length, self.num_atoms, 3))
        self._forces[step_no] = forces

    def set_basis(self, step_no, bas):
        self.basis[step_no] = bas


    def set_e_total(self, step_no, etot):
        self.total_energy[step_no] = etot 

    def set_e_kinetic(self, step_no, ekin):
        self.kinetic_energy[step_no] = ekin

    #------------------- End Properties -------------------------------

    def save(self, filename):
        """ 
        Saves the trajectory to a *.npz binary file using numpy.savez 
        """
        np.savez(filename, basis = self.basis, tote = self.total_energy, kine = self.kinetic_energy, pos = self._positions, 
            forces = self._forces, time = self._time, atoms = self.atoms)
        print "Wrote trajectory (%d atoms, %d steps) to %s" % (self.num_atoms, self.length, filename)

    def load(self, filename):
        """ 
        Loading the trajectory from a *.npz numpy binary file using numpy.load 
        """
        npz = np.load(filename)
        self.atoms = npz['atoms']
        self.basis = npz['basis']
        self.total_energy = npz['tote']
        self.kinetic_energy = npz['kine']
        self._positions = npz['pos']
        if 'forces' in npz:
            self._forces = npz['forces']
        else:
            print "This file does not contain forces"
            self._forces = np.zeros(self._positions.shape) # waste of memory
        self._time = npz['time']
        self.length = self._positions.shape[0]
        self.num_atoms = self._positions.shape[1] 
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
        ----------
        
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
            if np.max(dr) > 1:
                print "Warning: anomalously large displacement of %.1f between step %d and %d" % (dr,stepno,stepno-1)

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

    def add_atom(self, positions):
        """
        Adds an atom to the trajectory. The new atom must have the same data length as the existing atoms.
        Note that we do not check if the new atom has periodic boundary conditions or not. This may lead to 
        inconsistency with the excisting atoms if the periodic boundary conditions have been unwrapped.

        Returns the atom id of the new atom
        """
        if positions.shape[0] != self.length:
            print "Error: The length of the position array for the new atom is %d, while the length of the original trajectory is %d" % (positions.shape[0],self.length)
            return -1
        if positions.shape[1] != 3:
            print "Error: The positions array must be an (n,3) array."
            return -1
        p = np.reshape(positions,(positions.shape[0],1,3))
        self.positions = np.append(self.positions, p, axis = 1)
        return self.positions.shape[1]-1 # atom id of the new atom


    def remove_atom(self, atom_id):
        """
        Removes an atom from the trajectory.
        Not that this will lower the id's of all the remaining atoms.
        """
        self.positions = np.delete(self.positions, atom_id, axis = 1)
        return True


    def get_geometric_center(self, atoms):
        """
        Replaces the atoms given as a tuple by a new pseudo-atom.
        Returns the id of the pseudo-atom as an int. 

        Note: This function does not take periodic boundary conditions into account!
        """
        geo = np.sum(self.positions[:,atoms],axis=1) / len(atoms)
        return geo

    
    def get_single_particle_trajectory(self, atom_no, coords = 'direct'):
        """
        Returns the trajectory of a single atom in direct or cartesian coordinates
        as a numpy array of dimension (steps, 3).

        Parameters
        ----------
            coords : either 'direct' or 'cartesian'
        """
        pos = self.positions[:,atom_no,:]
        if coords == 'cartesian':
            pos = oppvasp.direct_to_cartesian(pos, self.basis)
        return pos


    def get_all_trajectories(self, coords = 'direct'):
        """
        Returns the trajectories for all the atoms in direct or cartesian coordinates
        as a numpy array of dimension (steps, atoms, 3).
        Parameters:
            - coords : either 'direct' or 'cartesian'
        """
        if coords[0].lower() == 'd':
            return self.positions
        elif coords[0].lower() == 'c':
            if 'positions_cart' in dir(self):
                return self.positions_cart
            else:
                pos = oppvasp.direct_to_cartesian(self.positions, self.basis)
                # This conversion is currently stupidly slow, so we cache the cartesian coordinates
                # as a temporary solution until a faster algorithm is implemented
                self.positions_cart = pos
                return pos


    def get_coordinates(self, coords = 'direct'):
        """
        Returns the coordinates$r(t)$ as a numpy array of shape (#steps, #atoms, 3)
        """
        return self.get_all_trajectories(coords)
    
    def get_displacements(self, coords = 'direct'):
        """
        Returns the displacements $r(t)-r(0)$ as a numpy array of shape (#steps, #atoms, 3)
        """
        pos = self.get_all_trajectories(coords)
        return pos - pos[0]


    def get_displacements_squared(self, coords = 'direct'):
        """
        Returns the displacements squared $[r(t)-r(0)]^2), as an numpy array:
            axis[0] = time, 
            axis[1] = atom_no,
        $ \Delta r_i^2(t) = (r_i(t)-r_i(0))^2 $
        Parameters:
         - coords : either 'direct' or 'cartesian'
        """
        disp = self.get_displacements(coords)
        dispsq  = disp[:,:,0]**2 + disp[:,:,1]**2 + disp[:,:,2]**2
        return dispsq
        # transpose to get axis[0] = at no., axis[1] = time

    def get_velocities(self, coords = 'direct', algorithm = 'naive'):
        """
        Returns velocities as a (nsteps,natom,3) array
        Parameters:
         - coords : 'direct' returns velocities in [direct coordinates]/[timestep]
                        'cartesian' returns velocities in Å/s
        """
        timestep = self.time[1]-self.time[0]
        pos = self.get_all_trajectories(coords)
        vel = np.zeros(pos.shape)

        pbar = ProgressBar(widgets=['Calculating velocities...',Percentage()], maxval = vel.shape[0]).start()
        if algorithm == 'naive':
            for i in range(1,vel.shape[0]-1):
                pbar.update(i)
                vel[i] = (pos[i+1] - pos[i-1]) / 2.*timestep 
        pbar.finish()

        return vel

    def get_occupancies(self, lattice, step_size = 100):
        """
        This method generates
         (1) A (time, nsites) array with the occupancy of each lattice site at each time. 
             The occupancy of a given site is defined as the number of atoms whose positions are 
             closer to this site than any other (periodic boundary conditions are taken into account).

        Parameters:
            lattice: (nsites,3) numpy array containing the coordinates of the perfect lattice sites in direct coordinates
            step_size: (int) only include every 'step_size' step (used to speed up things when reading long trajectories)

        Returns:
            occupancies numpy array
        """

        natoms = self.num_atoms
        nsites = lattice.shape[0]
        n = self.time[::step_size].shape[0]

        pbar = ProgressBar(widgets=['Calculating occupancies...',Percentage()], maxval = n).start()
        occupancies = np.zeros((n,nsites), dtype=int)
        inhabitants = np.zeros((n,nsites,4)) # max four atoms per site

        # for each MD step
        for i in range(n):
            pbar.update(i)
            stepno = i * step_size
            
            # for each atom
            for atno in range(natoms):
                atpos = self.positions[stepno,atno] # note: direct coordinates!
                
                # Make a (nsites,3) matrix with vectors pointing from the atom to each lattice site
                dist = lattice - atpos
                
                # Make a (3,nsites,3) tensor with one periodic image in each of the directions +x,-x,+y,-y,+z,-z:
                dist_i = np.abs(np.array((dist-1,dist,dist+1))) 
                # and find the shortest distance:
                dist = np.min(dist_i, axis=0)

                # Make a (nsites) vector with the distances squared (r^2)
                dist_r2 = np.sum(dist**2, axis=1)
                
                # Find closest lattice site(s):
                closest_site = np.argmin(dist_r2)
                #closest_sites = np.argsort(dist_r2)

                # Update occupancies array:
                inhabitants[i,closest_site,occupancies[i,closest_site]] = atno
                occupancies[i,closest_site] += 1
           
        pbar.finish()
        return occupancies, inhabitants

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


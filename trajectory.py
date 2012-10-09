# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os,sys
import numpy as np
#import scipy as sp
from copy import copy

import tables

import time
import oppvasp

from oppvasp import direct_to_cartesian
from oppvasp import elements
from oppvasp.structure import Structure
from oppvasp.util import unitcell_components

# Status of optional Python modules:
imported = { 'progressbar' : False, 'psutil' : False, 'lxml' : False }
try:
    from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
        RotatingMarker, ReverseBar, SimpleProgress
    imported['progressbar'] = True
except ImportError:
    print "Info: Module 'progressbar' is not available"


class Trajectory(object):
    """
    The Trajectory class holds positions, forces, kinetic and total energy 
    for a number of atoms over a number of time frames.

    Positions and forces are stored internally as `direct' coordinates with respect 
    to a unit cell, that may change over time (the unit cell is specified per frame).

    The information is saved and read natively using the Numpy npz binary format,
    or exported to the vtf format for use with the VMD viewer.

    The class provides some simple manipulation methods.

        pass either:
        - filename:  a pytables h5table file
        or:
        - num_steps:
        - timestep:
        - atoms:
    """

    def __init__(self, filename = '', num_steps = 0, timestep = 1., atoms = []):
        self.filename = filename
        self.path = ''
        self.pbc_unwrapped = False
        if filename != '':
            self.path = os.path.abspath(filename)
        if filename != '':
            self.load(filename)
        else:
            if len(atoms) == 0:
                raise StandardError("No filename or atoms specified")
            self.atoms = atoms
            #self.num_atoms = atoms.shape[0]
            self._time = np.zeros(num_steps)
            self.basis = np.zeros((self.length,3,3))
            self.total_energy = np.zeros(self.length)
            self.kinetic_energy = np.zeros(self.length)
            self._positions = np.zeros((num_steps, self.num_atoms, 3))
            self._forces = np.zeros((num_steps, self.num_atoms, 3))
            self.set_timestep(timestep)
            #print "Allocated %.2f MB for Trajectory object" % ((self._time.nbytes + self.total_energy.nbytes + self.kinetic_energy.nbytes + self._positions.nbytes + self._forces.nbytes)/(1024.**2))
        # Set default selection to the whole trajectory:

    def update_length(self, new_length):
        if new_length < self.length:
            print "Stripping %d empty steps" % (self.length - new_length)
            #self.length = new_length 
            self._time = self.time[0:new_length] # updates length
            self.total_energy = self.total_energy[0:new_length]
            self.kinetic_energy = self.kinetic_energy[0:new_length]
            self.basis = self.basis[0:new_length]
            self._positions = self._positions[0:new_length]
            self._forces = self._forces[0:new_length]

    def get_copy(self):
        return self.get_selection()

    def get_selection(self, start = 0, stop = 0, step = 1):
        """
        Returns a part of the trajectory as a new Trajectory object. 

        Parameters
        ----------
        start : int
            Start of interval. The interval includes this value.
            The default start value is 0.
        stop : int
            End of interval.  The interval does not include this value.
            The default stop value is 0, which will evaluate as including the last frame.
        step : int
            Spacing between values.  For any output `out`, this is the distance
            between two adjacent values, ``out[i+1] - out[i]``. 
            The default step size is 1.

        Examples
        ----------
        >>> traj = Trajectory( 'trajectory.npz' )
        >>> sel = traj.get_selection(10000, 20000)
        >>> d = sel.get_displacements()
        >>> t = sel.time
        >>> plt.plot(t,d)

        """

        if stop == 0:
            stop = self.nsteps
        p = self._positions[start:stop:step]
        nsteps = p.shape[0]
        t = Trajectory( num_steps = nsteps, atoms = self.atoms ) 
        t.set_timestep((self.time[1]-self.time[0])*step)
        t.positions = copy(p)
        t.forces = self._forces[start:stop:step]
        t.basis = self.basis[start:stop:step]
        t.total_energy = self.total_energy[start:stop:step]
        t.kinetic_energy = self.kinetic_energy[start:stop:step]
        return t
    
    def set_selection(self, in_point = 0, out_point = -1):
        raise DeprecationWarning("set_selection is deprecated. Use get_selection instead")
        #self.in_point = in_point
        #if out_point == -1:
        #    out_point = self.length # to include the very last element
        #self.out_point = out_point
        




    #------------------- Begin Properties -------------------------------

    @property
    def natoms(self):
        """number of atoms"""
        return self.atoms.shape[0]

    @property
    def num_atoms(self):
        """number of atoms"""
        return self.atoms.shape[0]
        #return self._positions.shape[1]

    @property
    def nsteps(self):
        """number of frames"""
        return self._time.shape[0]

    @property
    def length(self):
        """number of frames"""
        return self._time.shape[0]

    @property
    def time(self):
        """nsteps array"""
        return self._time
    
    @time.setter
    def time(self,val):
        self._time = val

    @property
    def timestep(self):
        """assumes linear time"""
        return self._time[1] - self._time[0]

    def set_timestep(self, timestep):
        self._time = np.arange(self.length) * timestep

    @property
    def positions(self):
        """
        nsteps x natoms x 3 array containing the direct coordinates of all atoms for all steps
        """
        return self._positions

    @positions.setter
    def positions(self, val):
        self._positions = val

    def set_positions(self, step_no, pos):
        if step_no >= self.length:
            #print "Warning: This trajectory was initiated to contain %d steps. Skipping step %d." % (self.length, step_no+1)
            return
        self._positions[step_no] = pos

    @property
    def forces(self):
        """
        nsteps x natoms x 3 array containing the forces on all atoms for all steps
        in eV/Angstrom
        """
        return self._forces

    @forces.setter
    def forces(self, val):
        self._forces = val

    def set_forces(self, step_no, forces):
        if step_no >= self.length:
            #print "Warning: This trajectory was initiated to contain %d steps. Skipping step %d." % (self.length, step_no+1)
            return
        self._forces[step_no] = forces

    def set_basis(self, step_no, bas):
        if step_no >= self.length:
            print "Warning: This trajectory was initiated to contain %d steps. Skipping step %d." % (self.length, step_no+1)
            return
        self.basis[step_no] = bas


    def set_e_total(self, step_no, etot):
        if step_no >= self.length:
            #print "Warning: This trajectory was initiated to contain %d steps. Skipping step %d." % (self.length, step_no+1)
            return
        self.total_energy[step_no] = etot 

    def set_e_kinetic(self, step_no, ekin):
        if step_no >= self.length:
            #print "Warning: This trajectory was initiated to contain %d steps. Skipping step %d." % (self.length, step_no+1)
            return
        self.kinetic_energy[step_no] = ekin

    #------------------- End Properties -------------------------------

    def save(self, filename, save_forces = False):
        """ Saves the trajectory to a hdf5 binary file using PyTables """

        # not sure if 8 bytes is correct on all systems
        mbperfloat = 8./1024**2

        b1 = self.nsteps * mbperfloat
        b2 = b1 * self.natoms * 3

        print "Approximate size needed:"
        print "  positions: %.2f MB" % (b2)
        if save_forces:
            print "  forces: %.2f MB" % (b2)
        print "  ekin: %.2f MB" % b1
        print "  etot: %.2f MB" % b1
        print "  time: %.2f MB" % b1
        print "  unit cell: %.2f MB" % (b1*3)
        if save_forces:
            print " Total: %.2f MB" % (6*b1 + 2*b2)
        else:
            print " Total: %.2f MB" % (6*b1 + b2)
        
        self.filename = filename
        self.path = os.path.abspath(filename)
        ext = os.path.splitext(filename)[1]
        if ext == '.h5':
            h5file = tables.openFile(filename, mode='w', title='Oppvasp Trajectory File')
            h5file.createArray(h5file.root, 'metadata', [['path', self.path], ['pbc_unwrapped', '1' if self.pbc_unwrapped else '0']])
            h5file.createArray(h5file.root, 'basis', self.basis)
            h5file.createArray(h5file.root, 'positions', self._positions)
            if save_forces:
                h5file.createArray(h5file.root, 'forces', self._forces)
            h5file.createArray(h5file.root, 'time', self._time)
            h5file.createArray(h5file.root, 'atoms', self.atoms)
            h5file.createArray(h5file.root, 'ekin', self.kinetic_energy)
            h5file.createArray(h5file.root, 'etot', self.total_energy)
            h5file.close()
        elif ext == '.npz':
            if save_forces:
                np.savez(filename, 
                    basis = self.basis, tote = self.total_energy, kine = self.kinetic_energy, 
                    pos = self._positions, forces = forces,
                    time = self._time, atoms = self.atoms)
            else:
                np.savez(filename, 
                    basis = self.basis, tote = self.total_energy, kine = self.kinetic_energy, 
                    pos = self._positions, 
                    time = self._time, atoms = self.atoms)
        else:
            raise StandardError("Unknown file extension!")

        print "Wrote trajectory (%d atoms, %d steps) to %s" % (self.num_atoms, self.length, filename)

    def load(self, filename):
        """ Reads trajectory from a hdf5 binary file using PyTables or a Numpy npz file"""

        ext = os.path.splitext(filename)[1]
        has_forces = 'no'
        has_metadata = 'no'
        if ext == '.h5':
            tbl = tables.openFile(filename, mode='a')
            self._time = tbl.getNode(tbl.root, 'time').read()
            self.atoms = tbl.getNode(tbl.root, 'atoms').read()
            self.basis = tbl.getNode(tbl.root, 'basis').read()
            self.total_energy = tbl.getNode(tbl.root, 'etot').read()
            self.kinetic_energy = tbl.getNode(tbl.root, 'ekin').read()
            self._positions = tbl.getNode(tbl.root, 'positions').read()
            try:
                options = tbl.getNode(tbl.root, 'metadata').read()
                for key, val in options:
                    if key == 'pbc_unwrapped':
                        self.pbc_unwrapped = (val == '1')
                    elif key == 'path':
                        self.path = val
                        if not os.path.isfile(self.path):
                            print "Warning: Original vasprun.xml file not found at %s" % self.path
                            # check if file is in current directory
                            b = os.path.basename(self.path)
                            a = os.path.abspath(b)
                            if os.path.isfile(a):
                                print 'File found at %s. Would you like to update the trajectory file? [y/n]' % a
                                if raw_input()[0].lower() == 'y':
                                    for row in tbl:
                                        if row[0] == 'path':
                                            row[1] = a
                                            self.path = a
                                            row.update()
                                            print "Ok, saved new path: %s" % a
                                            break

                            # check if file is in current directory
                has_metadata = 'yes'
            except tables.NoSuchNodeError:
                pass
                #print "File does not contain metadata"
            try:
                self._forces = tbl.getNode(tbl.root, 'forces').read()
                has_forces = 'yes'
            except tables.NoSuchNodeError:
                #print "This file does not contain forces"
                self._forces = np.zeros(self._positions.shape) # waste of memory
            tbl.close()
        elif ext == '.npz':
            npz = np.load(filename)
            self._time = npz['time']
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
            #self.length = self._positions.shape[0]
        else:
            raise StandardError("Unknown file extension!")

        print "Read trajectory (%d atoms, %d steps, forces: %s, metadata: %s) from %s" % (self.num_atoms, self.length, has_forces, has_metadata, filename)

    def unwrap_pbc(self, remove_numerical_noise = False, init_pos = None):
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
        >>> traj = vp.get_positions()
        >>> traj.save('trajectory.npz')
        >>> traj.unwrap_pbc( init_pos = init_pos )
        >>> traj.save('trajectory_unwrapped.npz')
        
        """
        if self.length == 0:
            print "Warning: Trajectory is empty";
            return
        #if imported['progressbar']:
        #    pbar = ProgressBar(widgets=['Unwrapping PBCs...',Percentage()], maxval = self.length).start()

        dr = np.zeros((self.length, self.num_atoms, 3), dtype='float64')
        r = np.zeros((self.length, self.num_atoms, 3), dtype='float64')
        if init_pos == None:
            r[0] = self.positions[0]
        else:
            r[0] = init_pos
        for stepno in range(1, self.length):
            dr = self.positions[stepno] - self.positions[stepno-1]
            c = np.abs(np.array((dr-1, dr, dr+1)))
            diff = np.array([[np.argmin(c[:,i,j]) for j in range(3)] for i in range(self.num_atoms)]) - 1  # (-1,0,1)
            dr += diff

            r[stepno] = r[stepno-1] + dr
            if np.max(dr) > 1:
                print "Warning: anomalously large displacement of %.1f between step %d and %d" % (dr,stepno,stepno-1)           

            # Remove numerical noise
            if remove_numerical_noise:
                # For very long runs, I feared that summation of all the small positinal differences 
                # might not sum up exactly to the final positions due to inaccuracy in floating 
                # point arithmetic, and wrote this piece of code. 
                # However, it seems like my fear was quite unjustified. I've never noticed the noise
                # to increase above 1e14.
                drp = r[stepno] - self.positions[stepno]
                noise = drp - np.round(drp)
                r[stepno] -= noise
                if abs(np.max(drift)) > 0.5:
                    print "Warning: Max drift is",np.max(drift),". This is very high, and may lead to erroneous behaviour of the PBC unwrapping."

            #if imported['progressbar']:
            #    pbar.update(stepno)
        self.positions = r
        #if imported['progressbar']:
        #    pbar.finish()
        self.pbc_unwrapped = True


    def remove_drift(self):
        r = np.zeros((self.nsteps, self.num_atoms, 3), dtype='float64')
        r[0] = self.positions[0]
        drift = np.zeros((self.nsteps, 3), dtype='float64')
        for stepno in range(1, self.length):
            dr = self.positions[stepno] - self.positions[stepno-1]
            drift[stepno] = np.sum(dr, axis = 0) / self.natoms
            dr = dr - drift[stepno]
            r[stepno] = r[stepno-1] + dr
        self.positions = r
        return drift


    def wrap_pbc(self):
        """
        Wraps coordinates outside the unit cell into it
        """
        # Minimum image convention (may not be appropriate for extreme cell geometries!!)
        #self.positions = self.positions - (2*self.positions-1).astype(np.int) # works only in the range <-0.5, 1.5>
        self.positions = self.positions - np.floor(self.positions) # works in any range
        self.pbc_unwrapped = False


    def get_snapshot(self, frame):
        """
        Returns a Structure object for a single frame. This can in turn be saved.

        Parameters:
          - frame      : the frame id (int)

        """
        return Structure( cell = self.basis[frame], positions = self.positions[frame], coords = 'direct', atomtypes = self.atoms )

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

    
    #def get_single_particle_trajectory(self, atom_no, coords = 'direct'):
    #    """
    #    Returns the trajectory of a single atom in direct or cartesian coordinates
    #    as a numpy array of dimension (steps, 3).

    #    Parameters
    #    ----------
    #        coords : either 'direct' or 'cartesian'
    #    """
    #    pos = self.positions[:,atom_no,:]
    #    if coords == 'cartesian':
    #        pos = oppvasp.direct_to_cartesian(pos, self.basis)
    #    return pos


    def get_positions(self, coords = 'direct', atom = -1):
        """
        Returns positions for one (atom_no>=0) or all (atom_no=-1) atoms 
        in direct or cartesian coordinates as a numpy array 
        of dimension (steps, atoms, 3).
        
        Parameters:
            - coords : either 'direct' or 'cartesian'
            - atom   : (int) Set this to only return positions for a single atom.
        """
        
        if atom == -1: 
            p = self.positions.copy()
        else:
            p = self.positions[:,atom,:].copy()

        if coords[0].lower() == 'd':
            return p
        elif coords[0].lower() == 'c':
            return direct_to_cartesian(p, self.basis)


    def get_coordinates(self, coords = 'direct'):
        """
        Returns the coordinates$r(t)$ as a numpy array of shape (#steps, #atoms, 3)
        """
        raise DeprecationWarning("get_coordinates is deprecated. Use get_positions instead")
        return self.get_positions(coords)
   
    def get_displacements(self, coords = 'cart'):
        """
        Returns the displacements $r(t)-r(0)$ as a numpy array of shape (#steps, #atoms, 3)
        """
        pos = self.get_positions(coords)
        return pos - pos[0]

    def get_masses(self):
        """
        Returns array of atomic masses in atomic mass unit (amu)
        """
        return np.array([elements[i]['mass'] for i in self.atoms])  # amu

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
        Calculates velocities from the change in positions between adjacent frames.
        Returns a (nsteps,natom,3) array

        Parameters:
         - coords : 'direct' returns velocities in [direct coordinates]/fs
                        'cartesian' returns velocities in Å/fs
        """
        timestep = self.time[1]-self.time[0]
        #print timestep
        pos = self.get_positions(coords)
        vel = np.zeros(pos.shape)

        pbar = ProgressBar(widgets=['Calculating velocities...',Percentage()], maxval = vel.shape[0]).start()
        if algorithm == 'naive':
            for i in range(1,vel.shape[0]-1):
                pbar.update(i)
                vel[i] = (pos[i+1] - pos[i-1]) / 2. / timestep 
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

    #def get_nearest_sites(reference_poscar, atomid, resolution = 50, treshold = 0.2):



    def add_periodic_images(self, rx = [-1,0], ry = [-1,0], rz = [-1,0] ):

        a = self.atoms
        p = self.positions
        f = self.forces
        c = self.basis

        nx = len(rx)
        ny = len(ry)
        nz = len(rz)

        nsteps = p.shape[0] 
        natoms = p.shape[1]

        natoms2 = natoms*nx*ny*nz

        # we may run out of memory here, but better now than later
        # forces may be commented out if not used 
        # (or better, there should be an option to generally don't include forces to save mem)
        p2 = np.empty((nsteps,natoms2,3))
        #f2 = np.empty((nsteps,natoms2,3))
        f2 = np.repeat(f,nx*ny*nz,axis=1)
        #print p2.shape
        #print f2.shape
        print "%.2f MB allocated." % (p2.nbytes/(1024.**2))
        c2 = np.empty(c.shape)
        a2 = np.tile(a,nx*ny*nz)

        # This may take some time
        pbar = ProgressBar(widgets=[u'Calculating images... frame ',SimpleProgress()], maxval = nsteps).start()
        for s in range(nsteps):
            pbar.update(s)
            for i,x in enumerate(rx):
                for j,y in enumerate(ry):
                    for k,z in enumerate(rz):
                        idx = ny*nz*i + nz*j + k 
                        #print "Adding cell %d of %d" % (f+1,nx*ny*nz)
                        ids = idx*natoms
                        ide = (idx+1)*natoms
                        p2[s,ids:ide,0] = (x + p[s,:,0]) / nx
                        p2[s,ids:ide,1] = (y + p[s,:,1]) / ny
                        p2[s,ids:ide,2] = (z + p[s,:,2]) / nz
                        #f2[s,ids:ide,0] = (x + f[s,:,0]) / nx
                        #f2[s,ids:ide,1] = (y + f[s,:,1]) / ny
                        #f2[s,ids:ide,2] = (z + f[s,:,2]) / nz
            c2[s,0] = c[s,0] * nx
            c2[s,1] = c[s,1] * ny
            c2[s,2] = c[s,2] * nz

        self.atoms = a2
        self.positions = p2
        self.forces = f2
        self.basis = c2
        pbar.finish()
        print "Added %d images of the original unit cell. The number of atoms is now %d." % (nx*ny*nz-1, self.num_atoms)

    def _vtfAtomDef(self,b,e,c):
        if b == e-1:
            str = 'atom %d' % (b)
        else:
            str = 'atom %d:%d' % (b,e)
        str += ' radius %.1f atomicnumber %d\n' % (elements[c]['vdw'],c)
        # the atomic radius is not used by VMD, but is part of the vtf file spec, so I add it
        return str

    def save_as_vtf(self, vtffile):
        """
        Saves the trajectory as a vtf file, that can be read by e.g. VMD. The vtf format is easy 
        to write (and read), and flexible in that it allows for addition of custom "userdata", that 
        can be easily read by newer versions of VMD's vtfplugin (https://github.com/olenz/vtfplugin)
        and used in VMD scripts. This method will save forces as userdata, but it could be modified
        to save other data instead (or as well).

        The vtf format is documented here:
        https://github.com/olenz/vtfplugin/wiki/VTF-format
        """
        
        if (isinstance(vtffile, file)):
            f = vtffile
        else:
            f = open(vtffile,'w')
        
        basis = self.basis
        atoms = self.atoms
        nsteps = self.length # NSW
        nions = self.num_atoms
        positions = self.positions
        forces = self.forces
        total_energy = self.total_energy

        # Write atom definitions
        f.write("# Structure block\n")
        b=0; e=0; c=-1
        for at in self.atoms:
            if c != at and e > 0:
                f.write(self._vtfAtomDef(b,e,c))
                b = e
            c = at
            e += 1
        f.write(self._vtfAtomDef(b,e-1,c))

        # This currently takes longer time than necessary I think. Should look into optimization
        pbar = ProgressBar(widgets=[u'Writing vtf... frame ',SimpleProgress()], maxval = nsteps).start()
        for i in range(nsteps):
            pbar.update(i)
            f.write("\n# Energy: %.3f\n" % total_energy[i])
            f.write("timestep\n")
            f.write("pbc %.2f %.2f %.2f %.1f %.1f %.1f\n" % unitcell_components(basis[i]))
            
            pos = np.dot(positions[i], basis[i])
            force = forces[i]
            for po,fo in zip(pos,force):
                f.write("%.8f    %.8f    %.8f     %.8f    %.8f    %.8f\n" % tuple(list(po) + list(fo)))
        pbar.finish()
        f.close()
        print "Wrote %s" % (f.name)
                
                    
                    



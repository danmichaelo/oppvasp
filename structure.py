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

    def __repr__(self):
        f = ' with forces' if self.has_forces() else ' no forces'
        return "%d atoms %s" % (self._positions.shape[0],f)
    
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
    
    def has_forces(self):
        """ Does the Structure contain forces? """
        return (self._forces != None)

    def get_forces(self, coords):
        """ (n,3) array """
        if not self.has_forces():
            return None
        if coords[0].lower() == 'd':
            return self._forces
        elif coords[0].lower() == 'c':
            return direct_to_cartesian(self._forces, self._cell)

    def set_forces(self, forces, coords):
        """
        Parameters:
            forces : (n,3) numpy array of forces on all atoms
            coords : 'direct' or 'cartesian'
        """
        if type(forces).__name__ == 'ndarray':
            self._forces = forces.copy()
        elif type(forces).__name__ == 'list':
            self._forces = np.array(forces)
        elif type(forces).__name__ == 'NoneType':
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
                    #print m
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

    
    def get_neighbours(self, atom_no, rmax = 3.0):
        """
        Returns a list of all atoms within <rmax> of atom <atom_no>,
        including atoms across periodic boundary conditions.
        """

        natoms = self.get_num_atoms()
        a = self._cell[0]
        b = self._cell[1]
        c = self._cell[2]

        # Find distance between opposing faces of the cell:
        wa = abs(np.dot(a,np.cross(b,c))) / np.linalg.norm(np.cross(b,c))
        wb = abs(np.dot(b,np.cross(c,a))) / np.linalg.norm(np.cross(c,a))
        wc = abs(np.dot(c,np.cross(a,b))) / np.linalg.norm(np.cross(a,b))

        # If half the max radius exceeds one of the face-to-face distances (wa, wb or wc), 
        # we add more atoms in the given direction(s) as necessary.
        la = 1 + int(2*rmax/wa)
        lb = 1 + int(2*rmax/wb)
        lc = 1 + int(2*rmax/wc)
        print "dim:",la,lb,lc
        #p2 = self.get_supercell_positions(la,lb,lc, coords='d')

        pos = np.zeros((natoms*la*lb*lc,3))
        offsets = np.zeros((natoms*la*lb*lc,3), dtype=np.int)
        #c = np.diag((1,1,1)) # since we use direct coordinates
        for i in range(la):
            for j in range(lb):
                for k in range(lc):
                    m = i*lb*lc + j*lc + k
                    pos[m*natoms:(m+1)*natoms] = self._positions + np.array((i, j, k))
                    offsets[m*natoms:(m+1)*natoms] = [i, j, k]

        indices = np.tile(np.arange(0,natoms), la*lb*lc)

        # Center at atom <index>:
        pos -= pos[atom_no]
        # Minimum image convention with coordinates in range lb*[-0.5,0.5>
        pos[:,0] -= la*(pos[:,0]*2/la).astype('int')
        pos[:,1] -= lb*(pos[:,1]*2/lb).astype('int')
        pos[:,2] -= lc*(pos[:,2]*2/lc).astype('int')
        
        cart = direct_to_cartesian(pos,self._cell)
        dr = np.sqrt(np.sum(cart**2,1))
        cond = np.logical_and(dr < rmax, dr > 0.0)

        dr = dr[cond]
        pos = pos[cond]
        offsets = (2*pos).astype('int')
        indices = indices[cond]

        # sort:
        so = np.argsort(dr)
        dr = dr[so]
        pos = pos[so]
        offsets = offsets[so]
        indices = indices[so]

        return indices, dr, offsets, pos

    
    def get_shortest_bond(self, include_self = False, safe_mode = True, rmax = 3.0):
        """
        Returns a tupple containing the shortest bond length in Angstrom,
        followed by the indices of the two atoms making up the bond
        
        Parameters:
            include_self : bool (default False)
                Whether to check for bonds between an atom and its _own_ periodic image.
                This is only needed for very small unit cells, typically with a single atom.
            safe_mode : bool (default True)
                Much slower, but safer for non-cubic cells
            rmax : float (default 3.0)
                Maximum bond length (in Angstrom) to check for 
                (only used in safe_mode)

        """
        a = self._cell[0]
        b = self._cell[1]
        c = self._cell[2]

        pos = self._positions
        natoms = pos.shape[0]
        if natoms == 1:
            print "Note from get_shortest_bond(): include_self is enabled since only a single atom was found!"
            include_self = True

        # Limit to sphere of radius Rc
        if safe_mode:
            #if include_self:
            #    print "include_self and safe_mode are not currently compatible. sorry"
            #    return
            
            natoms = self.get_num_atoms()
            a = self._cell[0]
            b = self._cell[1]
            c = self._cell[2]

            # Find distance between opposing faces of the cell:
            wa = abs(np.dot(a,np.cross(b,c))) / np.linalg.norm(np.cross(b,c))
            wb = abs(np.dot(b,np.cross(c,a))) / np.linalg.norm(np.cross(c,a))
            wc = abs(np.dot(c,np.cross(a,b))) / np.linalg.norm(np.cross(a,b))

            # If half the max radius exceeds one of the face-to-face distances (wa, wb or wc), 
            # we add more atoms in the given direction(s) as necessary.
            la = 1 + int(2*rmax/wa)
            lb = 1 + int(2*rmax/wb)
            lc = 1 + int(2*rmax/wc)
            print "dim:",la,lb,lc
            #p2 = self.get_supercell_positions(la,lb,lc, coords='d')
            pos = np.zeros((natoms*la*lb*lc,3))
            #offsets = np.zeros((natoms*la*lb*lc,3), dtype=np.int)
            #c = np.diag((1,1,1)) # since we use direct coordinates
            for i in range(la):
                for j in range(lb):
                    for k in range(lc):
                        m = i*lb*lc + j*lc + k
                        pos[m*natoms:(m+1)*natoms] = self._positions + np.array((i, j, k))
                        #offsets[m*natoms:(m+1)*natoms] = [i, j, k]
            indices = np.tile(np.arange(natoms), la*lb*lc)

            paired = np.zeros(natoms, dtype=np.int)
            minr = np.zeros(natoms)
            for j in np.arange(natoms):
                # Center at atom <index>:
                pos -= pos[j]
                # Minimum image convention with coordinates in range lb*[-0.5,0.5>
                pos[:,0] -= la*(pos[:,0]*2/la).astype('int')
                pos[:,1] -= lb*(pos[:,1]*2/lb).astype('int')
                pos[:,2] -= lc*(pos[:,2]*2/lc).astype('int')
                
                cart = direct_to_cartesian(pos,self._cell)
                dr = np.sqrt(np.sum(cart**2,1))
                cond = np.logical_and(dr < rmax, dr > 0.0)

                dr = dr[cond]
                if dr.shape[0] == 0:
                    print "[Warning]: No neighbours found within %.2f Angstrom of atom %d" % (rmax,j)
                    minr[j] = 1001.0
                    continue

                so = np.argsort(dr)

                #pos = pos[cond]
                #offsets = (2*pos).astype('int')
                ##indices = indices[cond]

                ## sort:
                #dr = dr[so]
                #pos = pos[so]
                #offsets = offsets[so]
                #indices = indices[so]

                #indices, dr, offsets, pos = self.get_neighbours(j)
                paired[j] = indices[cond][so][0]
                minr[j] = dr[so][0]
                #pairs[j,nearest] = 
                #minpair = pairs[j,nearest] 
                #minr = dr[0]
            
            #print nearest,shortest
            minidx = np.argmin(minr)
            minr = minr[minidx]
            if minr > 1000.0:
                print "[Error]: No bonds shorter than %.2f Angstrom found. Please increase rmax" % rmax

            return (minr,minidx,paired[minidx])

        else:

            # Build array of pair indices:
            pairs = get_pairs(natoms, include_self=include_self)
            npairs = pairs.shape[0]

            # Find displacement vectors for all pairs:
            x = pos[pairs[:,0]] - pos[pairs[:,1]]
        
            # Use minimum image convention to threat bonds over PBCs
            # Note: use safe_mode with tilted unit cells in which the 
            # shortest bond is close to half of one of the cell dimensions
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


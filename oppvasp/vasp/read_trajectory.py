# encoding=utf-8

import os

from oppvasp.trajectory import Trajectory
from oppvasp.vasp.parsers import PoscarParser, IterativeVasprunParser

def read_trajectory(dir='', unwrap_pbcs=True, xml_file='vasprun.xml',
        orig_file='traj_orig.h5', unwrapped_file='traj_unwrapped.h5',
        poscar_file='POSCAR'):
    """
    Convenience function for reading trajectory data from vasprun.xml.
    - The first time a trajectory is read, the trajectory data is read from
      vasprun.xml using lxml, and saved in a HDF5 binary file using PyTables. 
    - The next time the trajectory is read, the HDF5 file is read directly.
      This is much faster.

    If `unwrap_pbcs' is set to True, periodic boundary conditions (PBCs) will
    be unwrapped, using initial positions from a specified POSCAR file.
    Initial positions are taken from a POSCAR file instead of the initial
    positions in the vasprun.xml file, because the latter are wrapped into the
    cell, making visual comparison with the original structure difficult if the
    POSCAR file contained coordinates out of the cell, such as may result from
    a cell where the atoms were relaxed from initial positions at the edge of
    the cell.

    Parameters
    ----------
    dir : str
        Directory to read files from
        Default is current directory
    unwrap_pbcs : bool
        Set to True to unwrap the periodic boundary conditions (PBCs) or False to leave them.
        Default is True
    xml_file : str
        Name of the vasprun.xml file, default is "vasprun.xml"
    orig_file : str
        Name of the HDF5 file
    unwrapped_file : str
        Name of the HDF5 file containing unwrapped PBCs. 
        Only used if unwrap_pbcs is set to True

    """

    if not os.path.isfile(dir + orig_file):
        p = IterativeVasprunParser(dir + xml_file)
        traj = p.get_trajectory()
        traj.save(dir + orig_file)
        if not unwrap_pbcs:
            return traj
    elif not unwrap_pbcs:
        return Trajectory(filename = dir + orig_file)

    if os.path.isfile(dir + unwrapped_file):
        return Trajectory(filename = dir + unwrapped_file)
    else:
        poscar = PoscarParser(dir + poscar_file)
        pos = poscar.get_positions(coords = 'direct')
        try:
            t = traj
        except:
            t = Trajectory(filename = dir + orig_file)
        print "Unwrapping using initial pos from POSCAR"
        t.unwrap_pbc(init_pos = pos)
        t.save(dir + unwrapped_file)
        return t


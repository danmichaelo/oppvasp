#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import math
import os, sys
import numpy as np
from progressbar import ProgressBar, SimpleProgress

from oppvasp.plotutils import symmetric_running_mean, symmetric_running_median, prepare_canvas
from oppvasp import direct_to_cartesian
from plot_lattice_drift import AvgDistPlotter 
from plot_translational_order import TranslationalOrderPlot

from matplotlib.patches import Polygon
from matplotlib import cm
import pickle
import hashlib

import time,datetime

class Event:
    def __init__(self):
        self.handlers = set()

    def subscribe(self, handler):
        self.handlers.add(handler)
        return self

    def unsubscribe(self, handler):
        try:
            self.handlers.remove(handler)
        except:
            raise ValueError("Handler is not handling this event, so cannot unhandle it.")
        return self

    def fire(self, *args, **kargs):
        for handler in self.handlers:
            handler(*args, **kargs)

    def getHandlerCount(self):
        return len(self.handlers)

    __iadd__ = subscribe
    __isub__ = unsubscribe
    __call__ = fire
    __len__  = getHandlerCount


class DistanceFromLatticeSites(object):
    """
    This class calculates the distances between any atom in the trajectory and the n nearest 
    lattice site(s).
    The lattice sites are taken from the ''lattice'' variable, and the number of sites defined
    need not be equal to the number of atoms found in the trajectory.
    """

    def __init__(self, trajectory, lattice, atom = 0, nn = 8, step = 50, remove_lattice_drift = True, cache = True):
        """
        DistanceFromLatticeSites class
        
        Parameters:
            - trajectory : oppvasp.Trajectory object
            - lattice    : oppvasp.Structure object describing the lattice sites
            - atom       : atom to follow
            - nn         : number of nearest sites
            - step       : step size
        """

        if trajectory.pbc_unwrapped:
            print "----"
            print "DistanceFromLatticeSites: Probably best to use a trajectory with PBCs (no unwrapping)"
            print "----"

        self.entered_site = Event()
        self.left_site = Event()
        self.bottombar_labels = {}
        self.bottombar_height = .5
        self.bottombar_colors = {}
        self.valid = True

        self.lattice = lattice
        self.traj = trajectory.get_selection(step = step)
        if self.traj.length <= 1:
            print "[DistanceFromLatticeSites ERROR] Trajectory selection is too short. Aborting"
            self.valid = False
            return

        ts = self.traj.time[1] - self.traj.time[0]
        print "[DistanceFromLatticeSites] Timestamp for trajectory selection is %.f fs" % ts

        p = trajectory.path.split('/')
        basename = p.pop()
        cache_file = '/'.join(p) + '/' + basename.rsplit('.',1)[0] + '_dist_from_lattice.pkl'
        print "cache file:",cache_file
        
        if cache and os.path.isfile(cache_file):
            pkl_file = open(cache_file, 'rb')
            data = pickle.load(pkl_file)
            pkl_file.close()
            self.x = data['x']
            self.y = data['y']
            self.sites = data['sites']
            print "[DistanceFromLatticeSites] Loaded from cache %s" % cache_file
            return
            #checks?
            
        
        if remove_lattice_drift:
            print "[DistanceFromLatticeSites] Calculating drift using AvgDistPlotter"
            avg = AvgDistPlotter(trajectory = self.traj, lattice = self.lattice, step = 1)
            drift = avg.y_smoothed

        
        pos = self.traj.get_positions(coords = 'direct')
        print "Min/max: %.2f / %.2f" % (np.min(pos), np.max(pos))

        print "[DistanceFromLatticeSites] Unwrapping PBCs before smoothing"
        self.traj.unwrap_pbc()
        pos = self.traj.get_positions(coords = 'direct')
        self.traj.wrap_pbc()
        print "Min/max: %.2f / %.2f" % (np.min(pos), np.max(pos))

        print "[DistanceFromLatticeSites] Generating symmetric running mean..."
        atpos = pos[:,atom]
        atpos = symmetric_running_mean(atpos, 250/ts)
        atpos = symmetric_running_mean(atpos, 500/ts)

        # Wrap back into cell
        atpos = atpos - np.floor(atpos)
        print "Min/max: %.2f / %.2f" % (np.min(atpos), np.max(atpos))

        np.save('temp.npz', atpos)

        cell = self.traj.basis
        latpos = lattice.get_positions(coords = 'direct')

        natoms = pos.shape[1]
        nsteps = pos.shape[0]
        nsites = latpos.shape[0]
        #print "[DistanceFromLatticeSites] Reference lattice contains %d lattice sites" % (nsites)

        dists = np.zeros((nsteps, nn), dtype=np.float)
        sites = np.zeros((nsteps, nn), dtype=np.int)

        pbar = ProgressBar(widgets=['[DistanceFromLatticeSites] Analysing step ',SimpleProgress(),'...'], maxval = nsteps).start()
        for i in np.arange(0, nsteps):
            pbar.update(i)

            # Find displacement vectors between atom and all sites:
            r = atpos[i] - latpos

            # Use minimum image convention to threat bonds over PBCs
            # Note: This will not work with *very* tilted unit cells
            r = r - (2*r).astype('int')
            #r = r - np.floor(r)

            r = direct_to_cartesian(r, cell[i])
            if remove_lattice_drift:
                r -= drift[i] # cartesian? yes

            r2 = (r**2).sum(axis=1)
            r = np.sqrt(r2)

            idx = np.argsort(r)
            dists[i] = r[idx[0:nn]]
            sites[i] = idx[0:nn]

            #print "Shortest bond is: %.3f Angstrom (between atoms %d and %d)" % (minr,minpair[0],minpair[1])
            #break
        pbar.finish()

        #dists = symmetric_running_median(dists,6)
        #dists = symmetric_running_mean(dists,2)

        self.x = self.traj.time/1000.
        self.y = dists
        self.sites = sites

        pkl_file = open(cache_file, 'wb')
        pickle.dump({'x': self.x, 'y': self.y, 'sites': self.sites, 'trajectory': self.traj.path, 'lattice': lattice }, pkl_file)
        pkl_file.close()
        print "[DistanceFromLatticeSites] Wrote cache file: %s" % cache_file
    
    def find_melting(self, gvec = (1,1,1), crop = False):
        # Find melting, and crop trajectory if necessary. Also reove 1 ps
        to = TranslationalOrderPlot(self.traj, gvec = gvec)
        mp = to.find_melting()
        if mp != None:
            if crop:
                self.traj = traj.get_selection(stop = mp)
            print " -> Melted after %d ps" % (self.x[mp])

        # Check that trajectory is still non-empty after cropping
        if crop and traj.length <= 1:
            print "!! Trajectory empty. Skipping"
            print "[DistanceFromLatticeSites ERROR] Trajectory selection (before melting) is too short. Aborting"
            self.valid = False
            return

        return mp
    

    def get_index(self,lst,val):
        """Returns first index of val in lst"""
        return [i for i,x in enumerate(lst) if x == val][0]

    def get_occupancies(self, follow_atom = 0, closest_sites = 8, make_periodic_supercell = False):
        """
        This method generates
         (1) A (time, nsites) array with the occupancy of each lattice site at each time. 
             The occupancy of a given site is defined as the number of atoms whose positions are 
             closer to this site than any other (periodic boundary conditions are taken into account).

         (2) A (time,closest_sites) array of a given number of closest sites for atom 'follow_atom'
             at any time, and the distances to those sites.

        Parameters:
            follow_atom : (int) or list of ints - the atom id(s) to follow
            closest_sites : (int) the number of closest sites to include. 
            Increasing this value does not affect the processing time very much.

        Returns:
            A (occupancies, nearest_sites, nearest_sites_r2) tuple.
            'nearest_sites' is an array containing the nearest_sites (in decreasing order) for the followed atom,
            and 'nearest_sites_r2' the corresponding distances to these sites.
        """

        n = self.traj.time[::self.step_size].shape[0]
        occupancies = np.zeros((n, self.nsites), dtype=int)
        nclosest = closest_sites # plot <> closest sites
        p_site = np.zeros((len(follow_atom),n,nclosest), dtype=int)
        p_dist_r2 = np.zeros((len(follow_atom),n,nclosest))

        pbar = ProgressBar(widgets=['Frame...',SimpleProgress(),' (step size: %d)' % (self.step_size)], maxval = n*self.step_size).start()

        # for each step
        for i in range(n):
            stepno = i * self.step_size
            pbar.update(stepno)
            
            # for each atom
            for atno in range(self.natoms):
                dist = self.lattice - self.pos[stepno,atno]
                atpos = self.pos[stepno,atno]
                
                # Make a (nsites,3) array with vectors pointing from the atom to each lattice site
                dist = self.lattice - atpos

                # a) Minimum image convention work for most unit cells:
                #    (use direct coordinates)
                #
                #dist = dist - (2*dist-1).astype(np.int)
               
                # b) Alternative (slower) routine that can be used with very skewed unit cells: 
                #    (use cartesian instead of direct coordinates)
                #
                if make_periodic_supercell:
                    # Make a (3,nsites,3) tensor with one periodic image in each of the directions +x,-x,+y,-y,+z,-z:
                    dist_i = np.abs(np.array((dist-1,dist,dist+1))) 
                    # and find the shortest distance:
                    dist = np.min(dist_i, axis=0)

                # Make a (nsites) vector with the distances squared (r^2)
                dist_r2 = np.sum(dist**2, axis=1)
                
                # Find closest lattice site(s):
                #closest_site = np.argmin(dist_r2)
                closest_sites = np.argsort(dist_r2)

                if atno in follow_atom:
                    idx = self.get_index(follow_atom,atno)
                    p_site[idx,i] = closest_sites[0:nclosest]
                    p_dist_r2[idx,i] = dist_r2[p_site[idx,i]]

                # Update occupancies array:
                occupancies[i,closest_sites[0]] += 1
           
        pbar.finish()
        # (atno, stepno, nclosest)

        #p_sites = { }
        ## loop over atno
        #for atidx, at in enumerate(p_site[:,:,0]):
        #    for site in at[:,0]: # closest site
        #        if site in p_sites:
        #            p_sites[site] += 1
        #        else:
        #            p_sites[site] = 1
        #    invtot = 100./p_site.shape[0]
        #    print "-----------"
        #    keys = p_sites.keys()
        #    for k in np.argsort(p_sites.values())[::-1]:  # reverse order 
        #        site = keys[k]
        #        print "At site %d for %.1f%% of the time (%.2f,%.2f,%.2f)" % (site,p_sites[site]*invtot,self.lattice[site][0],self.lattice[site][1],self.lattice[site][2])
        #    print "-----------"


        return occupancies, p_site, p_dist_r2



#for i in range(n):
#    s = np.argsort(occupancies[i])
#    if occupancies[i,s[0]] == 0:
#        print "Vacancy found at site",s[0],"at step",i
#        #sys.stdout.write(str(s[0])+" ")
#    if occupancies[i,s[-2]] == 2:
#        print "Double occupation at site",s[0],"at step",i



#############################################################################
# (5) Plot

    def bottombar_fill(self, plt, index, t0, t, labs, verbose = False):
        if not index in self.bottombar_colors.keys():
            self.bottombar_colors[index] = cm.jet(np.random.rand(1.))[0]

        color = self.bottombar_colors[index]
        #print index,t0,t
        y0 = - self.bottombar_height
        plt.fill( [t0,t0,t,t], [0,y0,y0,0], facecolor = color, linewidth = 0.3, alpha = 0.5 )
        
        if index not in labs and t-t0 > 0.8:
            try:
                label = str(self.bottombar_labels[index])
            except KeyError:
                label = str(index)

            if verbose:
                print " -> Adding label",label,"(atom %d)"%(index)
            plt.text(t0 + (t-t0)/2., y0/2., label, horizontalalignment = 'center', verticalalignment = 'center', fontsize = 8)
            labs.append(index)
        return labs
    

    def plot(self, axes, verbose = False):
        """
        Automatic plot method. This method is not very robust, and should in general be modified 
        to highlight the points of interest.
        
        Parameters:
            axes: the matplotlib.Axes object to plot to

        """
        if not self.valid:
            return

        #fig = plt.figure()

        #p = [0.15, 0.57, 0.05, 0.03]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        #upper = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        #upper.set_ylabel(u'$r(t)-r(0)$ [Å]')
        #upper.set_xticklabels([])

        #p = [0.15, 0.15, 0.05, 0.45]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        #lower = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])

        #print self.followed_atom_r2.shape

        axes.plot(self.x, self.y, color = 'black')
        axes.set_xlabel('Time $t$ [ps]')
        axes.set_ylabel("$r'(t)-R$ [A]")
        axes.set_ylim(-self.bottombar_height, 2.6)

        #ticks = [0.]
        #ticklabels = [0.]
        #for i in range(self.y.shape[1]):
        #    axes.axhline(y[0,i], color = 'black', linestyle='dashed', alpha = 0.5)
        #    print "Y:",self.y[0,i]
        #    ticks.append(self.y[0,i])
        #    ticklabels.append("%.2f" % (y[0,i]))
        #axes.set_yticks(ticks)

        q = axes.get_yticks()
        q = q[q>0.]
        axes.set_yticks(q)

        #axes.set_yticklabels(ticklabels)

        #r = np.sqrt(np.sum((self.pos[::self.step_size,0]-self.pos[0,0])**2, axis=1))

        #print time.shape
        #print r.shape
        #upper.plot(time, r)

        # ================== Paint bottombar ==================================
        

        prev_site = -1
        cur_site = -1
        cur_count = 0
        labelled_sites = []
        nearest = 10.

        in_point = -1
        out_point = -1
        in_buffer = 10   # minimum number of frames to stay within <maxdist> of 
                         # a site to be considered having ''entered'' the site
        out_buffer = 10  # minimum number of frames to stay at least <maxdist> away from 
                         # the site to be considered having ''left'' the site
        buf = 0
        maxdist = 0.5
        
        for t in range(self.sites.shape[0]):
            
            # --- 1
            cur_site = self.sites[t,0]
            if cur_site != prev_site and prev_site != -1:
                #if nearest < 0.5:
                labelled_sites = self.bottombar_fill(axes, prev_site, self.x[t-cur_count], self.x[t], labelled_sites, verbose = verbose)
                cur_count = 0
                nearest = 10.
            if self.y[t,0] < nearest:
                nearest = self.y[t,0]
            prev_site = cur_site
            cur_count += 1

            # ---- 2
            if in_point == -1 and self.y[t,0] < maxdist:
                if buf > in_buffer:
                    # we entered a lattice site
                    in_point = t - buf
                else:
                    buf += 1

            elif in_point != -1 and self.y[t,0] >= maxdist:
                if buf > out_buffer:
                    # we left a lattice site
                    out_point = t - buf
                    x0 = self.x[in_point]
                    x1 = self.x[out_point]
                    y0 = 0
                    y1 = maxdist
                    self.entered_site.fire({ 'frame': in_point, 'x': self.x[in_point] })
                    #self.left_site.fire({ 'frame': in_point, 'x': self.x[in_point] })
                    axes.add_patch( Polygon([[x0,y0],[x1,y0],[x1,y1],[x0,y1]], 
                        facecolor='yellow', alpha=.6, linewidth=0.) )

                    #dr = 
                    #axes.annotate('H', xy= ((x1-x0)/2., y1*2.), xytext=((x1-x0)/2., y1*2.), xycoords='data', ha='center', va='center' )
                    #axes.text(x0 + (x1-x0)/2., 2*y1, 'H', ha = 'center', va = 'center') #, fontsize = 8)
                    if verbose:
                        print "At site between %.2d and %.2d" % (x0, x1)
                    in_point = -1

                else:
                    buf += 1
            else:
                buf = 0

        if in_point != -1:
            out_point = self.x.shape[0] - 1
            if out_point - in_point > in_buffer:
                x0 = self.x[in_point]
                x1 = self.x[out_point]
                y0 = 0
                y1 = maxdist
                self.entered_site.fire({ 'frame': in_point, 'x': self.x[in_point] })
                #self.left_site.fire({ 'frame': in_point, 'x': self.x[in_point] })
                axes.add_patch( Polygon([[x0,y0],[x1,y0],[x1,y1],[x0,y1]], 
                    facecolor='yellow', alpha=.6, linewidth=0.) )
                if verbose:
                    print "At site between %.2d and %.2d" % (x0, x1)
        
        
        axes.axhline(0.0, color = 'black')

        mp = self.find_melting(gvec=(1,1,1))
        if mp != None:
            axes.axvline(self.x[mp], color='red')
            x0 = self.x[mp]
            x1 = self.x[-1]
            y0 = 0.
            y1 = 3.0
            axes.add_patch( Polygon([[x0,y0],[x1,y0],[x1,y1],[x0,y1]], 
                        facecolor='red', alpha=.3, linewidth=0.) )
        
        axes.set_xlim(self.x[0], self.x[-1])
        
        #cur_count-=1
        labelled_sites = self.bottombar_fill(axes, prev_site, self.x[-cur_count], self.x[-1], labelled_sites)

        #self.plt = plt
        #self.= lower
        #self.upper = upper

    def plot_to_file(self, outfile = '', verbose = False):

        if not self.valid:
            return
        print "[DistanceFromLatticeSites] Plotting..."

        import matplotlib
        matplotlib.use('pdf')
        prepare_canvas('10 cm')
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.plot(ax, verbose = verbose)

        if outfile == '':
            outfile = os.path.splitext(sys.argv[0])[0] + '.pdf'
        sys.stdout.write("Writing %s... " % os.path.basename(outfile))
        sys.stdout.flush()
        plt.savefig(outfile)
        sys.stdout.write("done!\n")

    def time_at_lattice_sites(self, treshold = 0.5):

        mp = self.find_melting(gvec = (1,1,1))
        if mp == None:
            mp = self.y.shape[0]

        #print self.followed_atom_r2.shape  # (65, 600, 8)
        #y = np.sqrt(self.y[:mp,0])
        #print "y",y.shape
        print self.y.shape
        y = self.y[:mp]
        return float(y[y < treshold].shape[0]) / y.shape[0]
        #print y[not y>treshold].shape


if __name__ == "__main__":

    from oppvasp import read_trajectory
    from oppvasp.vasp.parsers import PoscarParser

    #atom = 42
    atom = 0

    traj = read_trajectory( './', unwrap_pbcs = False )
    poscar = '/Users/danmichael/Documents/Studier/Master/notur/hipersol/templates/si64/POSCAR'
    pl = DistanceFromLatticeSites(trajectory = traj, 
            lattice = PoscarParser(poscar).get_structure(),
            atom = atom, nn = 1, step = 50
            )
    first = True
    def on_enter_site(obj):
        global first
        if first:
            sys.stdout.write('\033[31;1mEntered lattice site at %.2f\033[0m\n' % obj['x'])
        first = False
    pl.entered_site.subscribe(on_enter_site)

    pl.plot_to_file('%s_%s.pdf' % (os.getcwd().split('/').pop(), 'dist_from_lattice' ), verbose = False)
    #dp.analyse()


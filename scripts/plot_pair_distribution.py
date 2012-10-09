#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0
from __future__ import unicode_literals
import sys
import os
from math import sqrt, sin

import numpy as np
import hashlib
import pickle

from oppvasp import read_trajectory
from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
from oppvasp import direct_to_cartesian

class PairDistributionFunctionPlot(object):

    def __init__(self, trajectory, n_bins = 10, atoms = [10], smin = 1000, smax = 20000, cache = None, 
            step = 5):
        
        struc = trajectory.get_snapshot(0)
        vol = struc.get_volume()
        s = sqrt(np.sum(struc.get_cell()[0]**2)) # sidelength a
        #max_dist = sqrt(3*(s/2.)**2) # max distance between atoms in 3D

        self.max_dist = s / 2. # the largest sphere that can be inscribed in a cube of side L

        print "Sidelength: %.2f Å, volume: %.2f Å^3, max r: %.2f Å" % (s, vol, self.max_dist)

        if n_bins == 'auto':
            rr = 0. # mean atom-atom distance
            for i in range(trajectory.natoms):
                rr += struc.get_neighbours(i)[1].mean()
            rr /= trajectory.natoms
            n_bins = int(self.max_dist / (0.025 * rr))
            # the value 0.025 rr from J. M. Haile. Molecular dynamics simulation : elementary methods. Wiley, 1992, p. 263
            print "Average atom-atom distance: %.2f Angstrom => setting n_bins = %d" % (rr, n_bins)

        
        self.dr = self.max_dist / n_bins
        print "Number of bins: %d, makes dr = %.2f Å" % (n_bins, self.dr)
        
        nsteps = trajectory.length
        natoms = trajectory.num_atoms
        
        self.density = natoms / vol
        
        if cache != None:
            hasher = hashlib.sha1()
            f = open(trajectory.filename, 'rb')
            try:
                hasher.update(f.read())
            finally:
                f.close()
            trajhash = hasher.hexdigest()
            print "Trajectory hash:",trajhash

            if os.path.isfile(cache):
                pkl_file = open(cache, 'rb')
                data1 = pickle.load(pkl_file)
                params = data1['params']
                if params['hash'] != trajhash:
                    print 'Note: Trajectory file hash differs from value in cache file "%s". Bypassing cache.' % cache
                elif params['n_bins'] != n_bins:
                    print 'Note: n_bins value differs from value in cache file "%s". Bypassing cache.' % cache
                elif params['atoms'] != atoms:
                    print 'Note: atoms value differs from value in cache file "%s". Bypassing cache.' % cache
                elif params['smin'] != smin:
                    print 'Note: smin value differs from value in cache file "%s". Bypassing cache.' % cache
                elif params['smax'] != smax:
                    print 'Note: smax value differs from value in cache file "%s". Bypassing cache.' % cache
                else:
                    # Hooray, we shall use the cache
                    self.x = data1['r']
                    self.y = data1['g']
                    pkl_file.close()
                    print 'Read data from cache "%s", trajectory has is "%s"' % (cache, trajhash)
                    return
                pkl_file.close()



        #d = 2 # super cell size 
        #fac = 2 * np.pi * d             # 2pi/0.125 = 2pi*8  (0.125 is the smallest distance between atoms in direct coordinates for the 2x2x2 Si cell)


        pos = trajectory.get_positions( coords = 'direct' )
        
        frames = np.arange(smin, np.min([nsteps,smax]), step)
        hist = np.zeros((frames.shape[0], n_bins))
        print hist.shape[0],"steps"
        
        for i, probe_at in enumerate(atoms):
            print "Probing atom %d of %d" % (i+1, len(atoms))

            for i, n in enumerate(frames):

                cpos = pos[n]
                ccell = trajectory.basis[n]

                r = cpos[[x for x in range(cpos.shape[0]) if x != probe_at]] - cpos[probe_at] # use probe_at as center

                # Use minimum image convention to threat bonds over PBCs
                # Note: This will not work with *very* tilted unit cells
                r = r - (2*r).astype('int')

                r = direct_to_cartesian(r, ccell)

                r2 = (r**2).sum(axis=1)
                r = np.sqrt(r2)

                cbin = np.floor(r / self.dr).astype('int')
                # could we vectorise this somehow?
                for x in cbin:
                    if x < n_bins:
                        hist[i,x] += 1


        # Normalize
        hist /= len(atoms)
        #norm = 2.0 * np.pi * dr * density

        self.x = np.zeros(n_bins)
        self.y = np.zeros(n_bins)

        for i in np.arange(0, n_bins):
            shell_r = i * self.dr                                                # distance to beginning of shell
            shell_vol = 4./3. * np.pi * ((shell_r + self.dr)**3 - shell_r**3)    # volume of the shell
            time_avg_n = np.sum(hist[:,i]) / hist.shape[0]                  # time average
            val = time_avg_n / shell_vol / self.density                          # local density to system density ratio

            #natoms / norm / ((rr**2) + (dr**2) / 12.0)
            self.x[i] = shell_r + 0.5*self.dr  # use middle of the shell
            self.y[i] = val

        #self.x = np.arange(n_bins) * dr

        #self.y = np.mean(hist[0:10000], axis = 0)
        ##self.y = np.mean(hist, axis = 0)

        #vmd_r = [0.05, 0.15000000000000002, 0.25, 0.35000000000000003, 0.45, 0.55, 0.65, 0.75, 0.8500000000000001, 0.9500000000000001, 1.05, 1.1500000000000001, 1.25, 1.35, 1.4500000000000002, 1.55, 1.6500000000000001, 1.75, 1.85, 1.9500000000000002, 2.0500000000000003, 2.15, 2.25, 2.35, 2.45, 2.5500000000000003, 2.6500000000000004, 2.75, 2.85, 2.95, 3.0500000000000003, 3.1500000000000004, 3.25, 3.35, 3.45, 3.5500000000000003, 3.6500000000000004, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.3500000000000005, 4.45, 4.55, 4.65, 4.75, 4.8500000000000005, 4.95, 5.050000000000001, 5.15, 5.25, 5.3500000000000005, 5.45, 5.550000000000001, 5.65, 5.75, 5.8500000000000005, 5.95]
        #vmd_g = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007519646352273967, 0.16803063700833396, 0.8705159618606413, 1.758370840475422, 2.283940404874828, 1.983953361507034, 1.4766080953857323, 1.1769551771463873, 0.9292292714634547, 0.8134830099213072, 0.7275062820612057, 0.7369909416061019, 0.7738143781341961, 0.8264792577807392, 0.9459870991621928, 1.008875524678367, 1.0764976921218232, 1.1637934824221774, 1.1568744636400938, 1.1417407274585984, 1.1471446817553927, 1.12710440892153, 1.081071932165155, 1.0372008450448278, 0.9681381745490185, 0.9209937945818178, 0.9260629623139609, 0.9387623185735751, 0.9248015478764121, 0.9193931865534244, 0.9593854795073877, 0.969338021724002, 1.0044774942407746, 1.0028481628862973, 1.0456729324756333, 1.0467550666063061, 1.0812861053364635, 1.0873968724468983, 1.102965726490678, 1.0960581045516773, 1.0686904640968937]


        #plt.plot(vmd_r, vmd_g, 'x-')
        #plt.show()
        
        if cache != None:
            # should we ask to overwrite?
            pkl_file = open(cache, 'wb')
            pickle.dump({
                'params': { 'hash': trajhash, 'n_bins': n_bins, 'atoms': atoms, 'smin': smin, 'smax': smax },
                'r': self.x,
                'g': self.y
            }, pkl_file)
            pkl_file.close()
            print 'Wrote cache "%s", trajectory has is "%s"' % (cache, trajhash)
            

    def getStructureFactor(self):

        self.q = np.arange(1.0, 8, 0.1)
        sq = np.zeros(self.q.shape[0])

        for i, q in enumerate(self.q):
            sqi = 0
            for r, g in zip(self.x, self.y):
                #print r, g
                sqr = (g - 1) * sin(q*r) / (q*r) * r**2 * self.dr
                sqi += sqr
            sq[i] = 4*np.pi*self.density * sqi

        self.sq = sq


    
    def plot_rdf(self, axes, data_style = {}, interpol_style = {}, **kwargs):
        """
        Automatic plot method. This method is not very robust, and should in general be modified 
        to highlight the points of interest.
        
        Parameters:
            axes: the matplotlib.Axes object to plot to

        """
        from scipy.interpolate import InterpolatedUnivariateSpline

        dst = { 'linestyle': '', 'marker': '.' } # default values
        ist = { 'linestyle': '-' } # default values

        for k,v in kwargs.items():
            dst[k] = v
            ist[k] = v
        for k,v in data_style.items():
            dst[k] = v
        for k,v in interpol_style.items():
            ist[k] = v


        print ist
        print dst

        iy = InterpolatedUnivariateSpline(self.x, self.y, k=3)
        xp = np.linspace(self.x[0], self.x[-1], 500)
        
        axes.axhline(1.0, linestyle = '--', color='gray')
        ip = axes.plot(xp, iy(xp), **ist)
        line = axes.plot(self.x, self.y, **dst)

        axes.set_xlabel('$r$ [Å]')
        axes.set_ylabel(u'$g(r)$')
        axes.set_xlim(0, self.max_dist)

        return line[0], ip[0]
    
    def plot_sq(self, axes):
        """
        Automatic plot method. This method is not very robust, and should in general be modified 
        to highlight the points of interest.
        
        Parameters:
            axes: the matplotlib.Axes object to plot to

        """
        axes.plot(self.q, self.sq, color = 'black', linestyle='-', marker='.')
        axes.axhline(0.0, linestyle = '--')
        axes.set_xlabel('$q$ [Å$^{-1}$]')
        axes.set_ylabel(u'$S(q)$')
        axes.set_xlim(0, self.q[-1])

    def plot_to_file(self, outfile = ''):

        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        self.plot_rdf(ax1)
        self.getStructureFactor()
        self.plot_sq(ax2)

        if outfile == '':
            outfile = os.path.splitext(sys.argv[0])[0] + '.pdf'
        sys.stdout.write("Writing %s... " % os.path.basename(outfile))
        sys.stdout.flush()
        plt.savefig(outfile)
        sys.stdout.write("done!\n")

if __name__ == '__main__':

    print " --> Reading trajectory"
    traj = read_trajectory(unwrap_pbcs = False)

    #for at in range(0,65):
    at = [10]
    pl = PairDistributionFunctionPlot(traj, n_bins = 'auto', atoms = at)
    pl.plot_to_file('%s_%s_%d.pdf' % (os.getcwd().split('/').pop(), 'rdf', at[0] ))

 

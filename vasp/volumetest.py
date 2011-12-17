#############################################################################
#
# @file volumetest.py @version 3
# This file should be called by <jobfile.sh>
# Last modified: Dec 17, 2011 21:53:00
#
# Example usage:
#
#   import os
#   import numpy as np
#   from oppvasp.vasp.volumetest import VolumeTestCubicUnitCell
#
#   job = VolumeTestCubicUnitCell(
#       lattice_parameters = np.arange(5.2, 5.8, 0.05),
#       basedir = './',
#       workdir = os.environ['SCRATCH'],
#       vaspcmd = 'vasp.x'
#   )
#   job.start(dry_run = False, first_step = 1)
#
#############################################################################
import os,sys
from job import BatchJob, BatchStep
from oppvasp import util

class VolumeTestCubicUnitCell(BatchJob):

    def __init__(self, lattice_parameters, **kwargs):
        """
        Create a VolumeTestCubicUnitCell object. 

        Parameters
        ----------
        lattice_parameters : array of lattice parameters to test.
           A set of POSCAR files named POSCAR.1, POSCAR.2, ... 
           will be created based on the file POSCAR with the 
           second line replaced with each of these values.

        The remaining parameters are passed over to BatchJob.
        Please see BatchJob for help on these.
        """
        BatchJob.__init__(self, **kwargs)

        self.paramName = 'VOL' # summary file header
         
        POSCAR = open('POSCAR', 'r')
        plines = POSCAR.readlines()
        POSCAR.close()

        for idx,param in enumerate(lattice_parameters):
            plines[1] = '%.4f\n' % float(param)
            ifile = open('POSCAR.%d' % (idx+1), 'w')
            ifile.writelines(plines)
            ifile.close()
            step = VolumeTestStep(idx+1, float(param))
            for template_name in BatchStep.input_files.keys():
                if template_name != 'POSCAR':
                    # use index-less 'template' file
                    step[template_name] = template_name
            self.add_step(step)

        self.print_info()

class VolumeTestStep(BatchStep):
    
    def __init__(self, index, paramValue):
        BatchStep.__init__(self, index)
        self.paramValue = paramValue

    def __str__(self):
        return "a = %.3f" % (self.paramValue)
    
    def get_name(self):
        return "%.3f" % (self.paramValue)


#############################################################################
#
# @file ktest.py @version 4
# This file should be called by <jobfile.sh>
# Last modified: Nov 26, 2010 19:20:17
#
# Example usage:
#
#   import os
#   from oppvasp.vasp.ktest import KTest 
#
#   analyzeOnly = ('vaspcommand' not in os.environ)
#   if analyzeOnly:
#       print "Environment variable 'vaspcommand' not set. Entering analyze-only mode."
#       basedir = os.path.curdir
#       vaspcmd = "ls" #dummy
#       workdir = '/dev/null' #dummy
#   else:
#       basedir = os.environ['SUBMITDIR']
#       vaspcmd = os.environ['vaspcommand']
#       workdir = os.environ['SCRATCH']
#
#   job = KTest(basedir,workdir,vaspcmd)
#   job.start(analyzeOnly)
#
#############################################################################
import os,sys,re
from batchjob import BatchJob, ManualBatchStep
from oppvasp import utils

class KTest(BatchJob):

    def __init__(self,basedir,workdir,vaspcmd):
        BatchJob.__init__(self,basedir,workdir,vaspcmd)
        self.paramName = 'K' # summary file header

        # Read parameter input
        self.parameterfile = 'ktest.in'
        os.chdir(self.basedir)
        if not os.path.isfile(self.parameterfile):
            print "Parameter-file '%s' not found!" % (self.parameterfile)
            sys.exit(1)
        f = open(self.parameterfile,'r')
        lines = f.readlines()
        f.close()
        ktype = lines[0]
        paramValues = []
        for l in lines[1:]:
            m = re.match('^[ \t]*([0-9])?[ \t]*([0-9])?[ \t]*([0-9])?', l)
            if m:
                paramValues.append(m.group(0))
        
        KPOINTS= open('KPOINTS', 'r')
        plines = KPOINTS.readlines()
        KPOINTS.close()

        for i in range(len(paramValues)):
            plines[3] = paramValues[i] + '\n'
            ifile = open('KPOINTS.%d' % (i), 'w')
            ifile.writelines(plines)
            ifile.close()
            self.addStep(KTestStep(i,paramValues[i]))

        self.info()

class KTestStep(ManualBatchStep):
    
    def __init__(self,index,paramValue):
        ManualBatchStep.__init__(self,index)
        self.paramValue = paramValue

    def __str__(self):
        return "K = %s" % (self.paramValue)
    
    def get_name(self):
        return "%s" % (self.paramValue)


#############################################################################
#
# @file convergencetest.py @version 3
# This file should be called by <jobfile.sh>
# Last modified: Dec 03, 2010 22:39:47
#
# Example usage:
#
#   import os
#   from oppvasp.vasp.convergencetest import ConvergenceTest
#
#   analyzeOnly = ('vaspcommand' not in os.environ)
#   if analyzeOnly:
#      print "Environment variable 'vaspcommand' not set. Entering analyze-only mode."
#      basedir = os.path.curdir
#      vaspcmd = "ls" #dummy
#      workdir = '/dev/null' #dummy
#   else:
#      basedir = os.environ['SUBMITDIR']
#      vaspcmd = os.environ['vaspcommand']
#      workdir = os.environ['SCRATCH']
#
#   job = ConvergenceTest(basedir,workdir,vaspcmd)
#   job.start(analyzeOnly)
#
#############################################################################
import os,shutil,sys,re,math
from batchjob import BatchJob, BatchStep

class ConvergenceTest(BatchJob):

    def __init__(self,basedir,workdir,vaspcmd, parameterfile = 'convergencetest.in', distributecmd = 'cp -Rupf'):
        BatchJob.__init__(self,basedir,workdir,vaspcmd,distributecmd)
        
        # Read parameter input
        self.parameterfile = parameterfile 
        os.chdir(self.basedir)
        if not os.path.isfile(self.parameterfile):
            print "Parameter-file '%s' not found!" % (self.parameterfile)
            sys.exit(1)
        f = open(self.parameterfile,'r')
        lines = f.readlines()
        f.close()
        param = lines[0].strip() #trim
        paramValues = [l.strip() for l in lines[1:]]
        self.paramName = param # summary file header
        
        # Add batch steps
        for i in range(len(paramValues)):
            self.addStep(ConvergenceTestStep(i,param,paramValues[i]))
       
        self.info()

class ConvergenceTestStep(BatchStep):

    def __init__(self,index,param,paramValue):
        BatchStep.__init__(self,index)
        self.param = param
        self.paramValue = paramValue

    def preProcess(self):
        # Update INCAR for the current run:
        f = open(self['INCAR'],'r+'); fc = f.read()
        fcr = re.sub(
            r'%s([ \t]*)=([ \t]*)([.\w]*)' % (self.param),
            '%s\\1=\\2 %s' % (self.param, self.paramValue),
            fc)
        if fc == fcr:
            fcr = fc+"\n %s = %s\n" % (self.param,self.paramValue)
        f.seek(0); f.write(fcr); f.truncate(); f.close();

    def postProcess(self):
        pass

    def __str__(self):
        return '%s=%s' % (self.param,self.paramValue)

    def getName(self):
        return self.paramValue

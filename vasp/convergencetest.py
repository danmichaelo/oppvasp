#############################################################################
#
# @file convergencetest.py @version 3
# This file should be called by <jobfile.sh>
# Last modified: Nov 17, 2010 18:46:42
#
# Usage:
#
#    analyzeOnly = ('vaspcommand' not in os.environ)
#    if analyzeOnly:
#        print "Environment variable 'vaspcommand' not set. Entering analyze-only mode."
#        basedir = os.path.curdir
#        vaspcmd = "ls" #dummy
#        workdir = '/dev/null' #dummy
#    else:
#        basedir = os.environ['SUBMITDIR']
#        vaspcmd = os.environ['vaspcommand']
#        workdir = os.environ['SCRATCH']
#
#   job = ConvergenceTest(basedir,workdir,vaspcmd)
#   job.summaryfile = 'summary.out'
#   job.parameterfile = 'convergencetest.in' 
#   job.start(analyzeOnly)
#
#
#############################################################################
import os,shutil,sys,re,math
from batchjob import BatchJob

class ConvergenceTest(BatchJob):

    def __init__(self,basedir,workdir,vaspcmd):
        BatchJob.__init__(self,basedir,workdir,vaspcmd)
        
        # Read parameter input
        self.parameterfile = 'convergencetest.in'
        os.chdir(self.basedir)
        if not os.path.isfile(self.parameterfile):
            print "Parameter-file '%s' not found!" % (self.parameterfile)
            sys.exit(1)
        f = open(self.parameterfile,'r')
        lines = f.readlines()
        f.close()
        param = lines[0].strip() #trim
        paramValues = [l.strip() for l in lines[1:]]
        
        # Add batch steps
        for i in range(len(paramValues)):
            self.addStep(ConvergenceTestStep(i,param,paramValues[i]))
       
        self.info()

class ConvergenceTestStep:

    def __init__(self,index,param,paramValue):
        self.index = index
        self.param = param
        self.paramValue = paramValue

    def preProcess(self):
        # Update INCAR for the current run:
        f = open('INCAR','r+'); fc = f.read()
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

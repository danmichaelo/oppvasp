#############################################################################
#
# @file convergencetest.py @version 3
# This file should be called by <jobfile.sh>
# Last modified: Dec 06, 2010 12:44:24
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
import numpy as np
import glob # for finding files using wildcards
from oppvasp.utils import query_yes_no
from oppvasp.vasp.parsers import vasprunParser
__docformat__ = "restructuredtext en"

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

class ConvergenceTestData():
    """
    This class parses all vasprun.xml files in the current directory 
    with the filename pattern 'vasprun*.xml'. A numpy array is created
    that holds the convergence parameter and the total energy from all the xml files.
    The array is sorted on volume and can be plotted using the GenericPlot class  
    or exported as CSV.

    Example usage:

    >>> data = ConvergenceTestData(parameter = 'ENCUT')
    >>> plot = GenericPlot('')
    >>> plot.addData(data)
    >>> plot.setXLabel('$E$')
    >>> plot.plot('out.pdf')
    >>> data.exportCSV('out.csv')
    
    """

    def __init__(self, directory = '', parameter = 'ENCUT'):

        #
        # Find and parse vasprun xml-files:
        #
        xmlFiles = glob.glob(directory+'vasprun*.xml')
        print "Found %d files matching 'vasprun*.xml'" % (len(xmlFiles))
        xy = np.zeros((2,len(xmlFiles)))
        for i in range(len(xmlFiles)):
            p = vasprunParser(xmlFiles[i])
            xy[0,i] = p.getIncarProperty(parameter)
            xy[1,i] = p.getTotalEnergy()


        #
        # Sort self.xy by increasing parameter value 
        #
        idx = np.argsort(xy[0])
        xy = np.array([[xy[j,i] for i in idx] for j in [0,1]])
        
        #
        # Print to screen
        #
        print "------------------------------------------------"
        print "File\t\t"+parameter+"\t\tEnergy"
        print "------------------------------------------------"
        for i in range(xy.shape[1]):
            print "%s\t%.2f\t\t%.4f" % (xmlFiles[idx[i]],xy[0,i],xy[1,i])
        print 

        self.xy = xy

    def exportCSV(self, filename = 'convergencetestplot.csv'):
        """
        Exports the internal x,y numpy array as a .csv file, 
        that can be imported into Excel or what have you.
        """
        sys.stdout.write("\nSaving table to %s... " % filename)
        sys.stdout.flush()
        f = open(filename,'w')
        for i in range(self.xy.shape[1]):
            f.write('%.6f\t%.6f\n' % (self.xy[0,i],self.xy[1,i]))
        f.close()
        sys.stdout.write("done!\n\n")

    def getData(self):
        return self.xy


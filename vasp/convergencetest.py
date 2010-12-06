import os,shutil,sys,re,math
from batchjob import BatchJob, BatchStep
import numpy as np
import glob # for finding files using wildcards
from oppvasp.utils import query_yes_no
from oppvasp.vasp.parsers import vasprunParser, outcarParser
__docformat__ = "restructuredtext en"

class ConvergenceTest(BatchJob):
    """
    Example usage:

    >>> import os
    >>> from oppvasp.vasp.convergencetest import ConvergenceTest
    >>>
    >>> analyzeOnly = ('vaspcommand' not in os.environ)
    >>> if analyzeOnly:
    >>>     print "Environment variable 'vaspcommand' not set. Entering analyze-only mode."
    >>>     basedir = os.path.curdir
    >>>     vaspcmd = "ls" #dummy
    >>>     workdir = '/dev/null' #dummy
    >>> else:
    >>>     basedir = os.environ['SUBMITDIR']
    >>>     vaspcmd = os.environ['vaspcommand']
    >>>     workdir = os.environ['SCRATCH']
    >>>
    >>> job.start(analyzeOnly)
    """

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

    def preprocess(self):
        """Updates INCAR for the current run"""
        f = open(self['INCAR'],'r+'); fc = f.read()
        fcr = re.sub(
            r'%s([ \t]*)=([ \t]*)([.\w]*)' % (self.param),
            '%s\\1=\\2 %s' % (self.param, self.paramValue),
            fc)
        if fc == fcr:
            fcr = fc+"\n %s = %s\n" % (self.param,self.paramValue)
        f.seek(0); f.write(fcr); f.truncate(); f.close();

    def postprocess(self):
        pass

    def __str__(self):
        return '%s=%s' % (self.param,self.paramValue)

    def get_name(self):
        return self.paramValue

class ConvergenceTestData():
    """
    This class will scan a directory for vasprun.xml or OUTCAR files and generate
    an numpy array holding the convergence parameter and the total energy from all the files.
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

    outcarPattern = 'OUTCAR*'
    vasprunPattern = 'vasprun*.xml'

    def __init__(self, directory = '', parameter = 'ENCUT', verbose = False, useVasprunXML = False):
        """
        :Parameters:
            directory : string
                directory to look for files in
            parameter : string
                INCAR parameter to look for (default is 'ENCUT')
            verbose : bool
                Print info about each step to screen 
            useVasprunXML : bool
                True to scan for vasprun.xml files, False to scan for OUTCAR files. 
        """
        self.directory = directory
        self.parameter = parameter
        self.verbose = verbose

        if useVasprunXML:
            self.readXMLfiles()
        else:
            self.readOUTCARfiles()
        
        # Sort self.xy and self.cpu by increasing parameter value 
        idx = np.argsort(self.xy[0])
        self.xy = np.array([[self.xy[j,i] for i in idx] for j in [0,1]])
        if hasattr(self,'cpu'):
            self.cpu = np.array([self.cpu[i] for i in idx])

        if verbose:
            print "------------------------------------------------"
            print "File\t\t"+parameter+"\t\tEnergy\t\tCPU"
            print "------------------------------------------------"
            for i in range(self.xy.shape[1]):
                if cpu in self:
                    print "%s\t%.2f\t\t%.4f\t\t%.f" % (xmlFiles[idx[i]],self.xy[0,i],self.xy[1,i],self.cpu[i])
                else:
                    print "%s\t%.2f\t\t%.4f\t\t-" % (xmlFiles[idx[i]],self.xy[0,i],self.xy[1,i])
            print 


    def readOUTCARfiles(self):
        """
        Finds and parses OUTCAR files 
        using the pattern defined in `ConvergenceTestData.outcarPattern`
        """
        files = glob.glob(self.directory + ConvergenceTestData.outcarPattern)
        print "Found %d files matching '%s%s'" % (len(files), self.directory, ConvergenceTestData.outcarPattern )
        self.xy = np.zeros((2,len(files)))
        self.cpu = np.zeros((len(files)))
        for i in range(len(files)):
            p = outcarParser(files[i])
            self.xy[0,i] = p.getIncarProperty(self.parameter)
            self.xy[1,i] = p.getTotalEnergy()
            self.cpu[i] = p.getCPUTime()

    def readXMLfiles(self):
        """
        Finds and parses OUTCAR files 
        using the pattern defined in `ConvergenceTestData.vasprunPattern`
        """
        files = glob.glob(self.directory + ConvergenceTestData.vasprunPattern)
        print "Found %d files matching '%s%s'" % (len(files), self.directory, ConvergenceTestData.vasprunPattern)
        self.xy = np.zeros((2,len(files)))
        for i in range(len(files)):
            p = vasprunParser(files[i])
            self.xy[0,i] = p.getIncarProperty(self.parameter)
            self.xy[1,i] = p.getTotalEnergy()

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

    def getTable(self):
        """Return xy numpy array"""
        return self.xy

    def getCPUTimes(self):
        """
        Returns numpy array of CPU times. Note that CPU times are not available from vasprun.xml files.
        """
        return self.cpu


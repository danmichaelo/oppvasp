import os,shutil,sys,re,math
from batchjob import BatchJob, BatchStep
import numpy as np
from oppvasp.utils import query_yes_no
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


class ConvergenceTestData(BatchJobDataExtractor):
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

    def __init__(self, directory = '', parameter = 'ENCUT', verbose = False, use_vasprun_xml = False, use_cache = True):
        """
        :Parameters:
            directory : string
                directory to look for files in
            parameter : string
                INCAR parameter to look for (default is 'ENCUT')
            verbose : bool
                Print info about each step to screen 
            use_vasprun_xml: bool
                True to scan for vasprun.xml files, False to scan for OUTCAR files. 
            use_cache: bool
                Cache as csv file
        """

        data_spec = [
            ['get_incar_property',parameter],
            ['get_total_energy'],
            ['get_cpu_time']
        ]
        BatchJobDataExtractor.__init__(self, directory = directory, verbose = verbose, use_vasprun_xml = use_vasprun_xml, 
                use_cache = use_cache, data_spec = data_spec )

        BatchJobDataExtractor.sort(self,0)

    
    #def export_csv(self, filename = 'convergencetestdata.csv'):
    #    """
    #    Exports the internal x,y numpy array as a .csv file, 
    #    that can be imported into Excel or what have you.
    #    """
    #    sys.stdout.write("Saving table to %s... " % filename)
    #    sys.stdout.flush()
    #    f = open(filename,'w')
    #    f.write(self.parameter+"\tEnergy\tCPU\n")
    #    for i in range(self.xy.shape[1]):
    #        if hasattr(self,'cpu'):
    #            f.write('%.6f\t%.6f\t%.f\n' % (self.xy[0,i],self.xy[1,i],self.cpu[i]))
    #        else:
    #            f.write('%.6f\t%.6f\t%.f\n' % (self.xy[0,i],self.xy[1,i],0.))
    #    f.close()
    #    sys.stdout.write("done!\n")

    #def import_csv(self, filename = 'convergencetestdata.csv'):
    #    f = open(filename,'r')
    #    all_lines = f.readlines()
    #    lines = all_lines[1:]   # skip first line
    #    self.xy = np.zeros((2,len(lines)))
    #    self.cpu = np.zeros(len(lines))
    #    for i in range(len(lines)):
    #        s = lines[i].split()
    #        self.xy[0,i] = float(s[0])
    #        self.xy[1,i] = float(s[1])
    #        self.cpu[i] = float(s[2])
    #    print "%d datapoints read from cache" % (len(lines))

    #def getTable(self):
    #    """Return xy numpy array"""
    #    return self.xy

    #def getCPUTimes(self):
    #    """
    #    Returns numpy array of CPU times. Note that CPU times are not available from vasprun.xml files.
    #    """
    #    return self.cpu


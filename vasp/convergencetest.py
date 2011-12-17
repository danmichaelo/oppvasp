import os,shutil,sys,re,math
from job import BatchJob, BatchStep, BatchJobDataExtractor
import numpy as np
from oppvasp.util import query_yes_no
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

    def __init__(self, parameters, **kwargs):
            
        BatchJob.__init__(self, **kwargs)

        self.params = parameters
        self.paramNames = parameters.keys()

        # Check lengths:
        npar = len(self.paramNames)
        nval = np.array([len(self.params[k]) for k in self.params])
        if not np.all(nval == nval[0]):
            print "ConvergenceTest initialization error: \n" + \
                " The number of values must be the same for all the parameters specified!"
            sys.exit(1)
        
        # Add batch steps
        for i in range(1,nval[0]+1):
            print "Adding step",i
            p = {}
            for k in self.paramNames:
                p[k] = str(self.params[k][i-1])
            step = ConvergenceTestStep(i,p)
            for template_name in BatchStep.input_files.keys():
                # use index-less 'template' file
                step[template_name] = template_name
            self.add_step(step)
       
        self.print_info()

class ConvergenceTestStep(BatchStep):
    """
    A single convergence job step.
    """

    def __init__(self, index, params):
        """
        Initializes a convergence test step.
        params is a dict of parameters with values, like { 'ENCUT': 200, 'SIGMA': 0.1 }
        """
        BatchStep.__init__(self, index)
        self.params = params
    
    def preprocess_info(self):
        """
        This method is called from the print_info method
        """
        print "  -> Update INCAR:",str(self)


    def preprocess(self):
        """
        Updates the INCAR file to prepare for the execution of this step 
        """
        f = open(self['INCAR'],'r+'); 
        incar = f.read() # read the whole file
        for key in self.params:
            val = self.params[key]
            search_str = r'%s(?P<ws1>[ \t]*)=(?P<ws2>[ \t]*)([.\w]*)' % (key)
            matches = re.search(search_str, incar)
            if matches:
                incar_mod = re.sub(
                    search_str,
                    r'%s\g<ws1>=\g<ws2>%s' % (key, val),
                    incar)
            else:
                # the parameter was not found. let's add it
                incar_mod = incar + '\n %s = %s\n' % (key, val)
            incar = incar_mod
        f.seek(0); 
        f.write(incar); 
        f.truncate(); 
        f.close();
        #print incar

    def postprocess(self):
        pass

    def __str__(self):
        return ','.join(('%s=%s' % (k,self.params[k]) for k in self.params.keys()))

    def get_name(self):
        return ','.join((str(self.params[k]) for k in self.params.keys()))



class ConvergenceTestData(BatchJobDataExtractor):
    """
    This class will scan a directory for vasprun.xml or OUTCAR files and generate
    an numpy array holding the convergence parameter and the total energy from all the files.
    The array is sorted on volume and can be plotted using the GenericPlot class  
    or exported as CSV.
    
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


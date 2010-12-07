import os,sys,re
from oppvasp import utils
from batchjob import BatchJob, ManualBatchStep

class KTest(BatchJob):
    """
    Sets up a k-point convergence test batch job using k-points specified in ktest.in
    Example usage:

    >>> import os
    >>> from oppvasp.vasp.ktest import KTest 

    >>> basedir = os.environ['SUBMITDIR']
    >>> vaspcmd = os.environ['vaspcommand']
    >>> workdir = os.environ['SCRATCH']

    >>> job = KTest(basedir,workdir,vaspcmd)
    >>> job.start()
    """

    def __init__(self, basedir, workdir, vaspcmd, parameterfile = 'convergencetest.in', distributecmd = 'cp -Rupf'):
        BatchJob.__init__(self,basedir,workdir,vaspcmd,distributecmd)
        self.paramName = 'K' # summary file header

        # Read parameter input
        self.parameterfile = parameterfile
        os.chdir(self.basedir)
        if not os.path.isfile(self.parameterfile):
            print "Parameter-file '%s' not found!" % (self.parameterfile)
            sys.exit(1)
        f = open(self.parameterfile,'r')
        lines = f.readlines()
        f.close()
        ktype = lines[0]
        param_values = []
        for l in lines[1:]:
            m = re.match('^[ \t]*([0-9])?[ \t]*([0-9])?[ \t]*([0-9])?', l)
            if m:
                param_values.append(m.group(0))
        
        f = open('KPOINTS', 'r')
        plines = f.readlines()
        f.close()

        for i in range(len(param_values)):
            plines[3] = param_values[i] + '\n'
            ifile = open('KPOINTS.%d' % (i), 'w')
            ifile.writelines(plines)
            ifile.close()
            self.addStep(KTestStep(i,param_values[i]))

        self.info()

class KTestStep(ManualBatchStep):
    
    def __init__(self,index,paramValue):
        ManualBatchStep.__init__(self,index)
        self.paramValue = paramValue

    def __str__(self):
        return "K = %s" % (self.paramValue)
    
    def get_name(self):
        return "%s" % (self.paramValue)


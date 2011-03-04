#############################################################################
#
# @file batchjob.py @version 2
# This file should be called by <jobfile.sh>
# Last modified: Mar 04, 2011 15:37:50
#
# Example usage:
#
#   import os
#   from oppvasp.vasp.batchjob import ManualBatchJob 
#
#   analyze_only = ('vaspcommand' not in os.environ)
#   if analyze_only:
#       print "Environment variable 'vaspcommand' not set. Entering analyze-only mode."
#       basedir = os.path.curdir
#       vaspcmd = "ls" #dummy
#       workdir = '/dev/null' #dummy
#   else:
#       basedir = os.environ['SUBMITDIR']
#       vaspcmd = os.environ['vaspcommand']
#       workdir = os.environ['SCRATCH']
#
#   job = ManualBatchJob(basedir,workdir,vaspcmd)
#   job.start(analyze_only)
#
#############################################################################
import os,shutil,sys,re,datetime
import numpy as np
import glob # for finding files using wildcards
from operator import itemgetter
from parsers import VasprunParser, OutcarParser, PoscarParser, IterativeVasprunParser
__docformat__ = "restructuredtext en"

class BatchJob:

    def __init__(self,basedir,workdir,vaspcmd,distributecmd = "cp -Ruf"):
        self.basedir = basedir
        self.workdir = workdir
        self.vaspcmd = vaspcmd
        self.distributecmd = distributecmd
        self.summaryfile = 'summary.out'
        self.steps = []
        self.paramName = 'Run' # summary file header

    def addStep(self,step):
        self.steps.append(step)
        step.basedir = self.basedir
        step.workdir = self.workdir
        step.vaspcmd = self.vaspcmd
        step.distributecmd = self.distributecmd
    
    def info(self):
        print "Total number of runs: %d " % (len(self.steps))
        print " Index\t Step info"
        for step in self.steps: 
            print " %4d\t %s" % (step.index,step)
        print


    def start(self, analyze_only = False, first_step = 0):

        os.chdir(self.basedir)

        # Initialize summary file
        if not os.path.isfile(self.summaryfile):
            sf = open(self.summaryfile,'w')
            sf.write('# File: '+os.path.abspath(self.summaryfile)+'\n')
            sf.write('# Created: '+datetime.datetime.today().strftime("%F %T")+'\n')
            sf.write(str(self.paramName)+"\t# pts\tDist.\tToten\t\tCPU\tForce\tPress.\tDrift\n")
            sf.close()

        # (4) Go!
        for step in self.steps:

            if (not analyze_only) and (step.index >= first_step):
                step.execute()

            # Analyze output and print results to summary file

            #vasprun = vasprunParser('vasprun.%d.xml' % (step.index))
            if not os.path.isfile(step['OUTCAR']):
                print "No output file '%s' to analyze. Exciting..." % (step['OUTCAR'])
                sys.exit(1)

            poscar = PoscarParser(step['POSCAR'])
            outcar = OutcarParser(step['OUTCAR'], selective_dynamics = poscar.selective_dynamics)
            outcar.readItAll()
            
            try:
                pressure = "%.4f" % outcar.get_max_pressure()
            except AttributeError:
                pressure = "  -  "


            summaryline = "%s\t%d\t%.3f\t%.4f\t%.0f\t%.4f\t%s\t%.4f" % (
                step.get_name(),
                outcar.get_num_kpoints(),
                outcar.dist,
                outcar.get_total_energy(),
                outcar.get_cpu_time(),
                outcar.get_max_force(),
                pressure,
                outcar.get_max_drift()
            )
            #print summaryline
            sf = open(self.summaryfile,'a')
            sf.write(summaryline+"\n")
            sf.close()


class BatchStep:

    def __init__(self,index):

        self.index = index
        self.inlist = ['INCAR','POSCAR','POTCAR','KPOINTS']
        self.outlist = ['OUTCAR','vasprun.xml','XDATCAR']
        # Set default values
        self.files = {}
        for i in self.inlist:
            self.files[i] = i
        self.files['OUTCAR'] = 'OUTCAR.%d' % (self.index)
        self.files['vasprun.xml'] = 'vasprun%d.xml' % (self.index)
        self.files['XDATCAR'] = 'XDATCAR.%d' % (self.index)
    
    def __setitem__(self,key,val):
        if not key in (self.inlist + self.outlist):
            print "%s is not a valid input/output filename" % (key)
            sys.exit(1)
        self.files[key] = val

    def __getitem__(self,key):
        if not key in (self.inlist + self.outlist):
            print "%s is not a valid input/output filename" % (key)
            sys.exit(1)
        return self.files[key]

    def preprocess(self):
        pass

    def postprocess(self):
        pass

    def execute(self): 

        # Prepare for the current run
        self.preprocess()

        # (Re-)Distribute files to all nodes:
        os.system("%s %s/* %s/" % (self.distributecmd,self.basedir, self.workdir))

        # Run VASP:
        os.chdir(self.workdir)
        os.system(self.vaspcmd)
        os.chdir(self.basedir)

        # Move files back
        os.system('cp -Rupf %s/* %s/ ' % (self.workdir, self.basedir))

        # Save some output files
        for f in self.outlist:
            if os.path.isfile(f):
                print "Renaming %s to %s" % (f, self[f])
                os.rename(f,self[f])

        # May be used to save more files
        self.postprocess()


class ManualBatchJob(BatchJob):

    def __init__(self,basedir,workdir,vaspcmd):
        BatchJob.__init__(self,basedir,workdir,vaspcmd)
        

        # Set up the different steps. If e.g. INCAR.0, INCAR.1 and INCAR.2 are found, 3 steps will be set up.
        while True: # bad habit...
            step = ManualBatchStep(len(self.steps))
            if step.necessary:
                self.addStep(step)
            else:
                break

        self.info()

class ManualBatchStep(BatchStep):
    
    def __init__(self,index):
        BatchStep.__init__(self,index)
        self.necessary = False
        for f in self.inlist:
            self[f] = f
            if os.path.exists(f+'.0'):
                d = 0
                while d <= index:
                    if os.path.exists('%s.%d' % (f,d)):
                        fi = '%s.%d' % (f,d)
                        if d == index:
                            self.necessary = True
                    d += 1
                self[f] = fi
            elif not os.path.exists(f):
                print "Neither a file '%s' nor a file '%s' was found. I'm terribly confused." % (f,f+'.0')
                sys.exit(1)


    def __str__(self):
        return " ".join(self.files.values())

    def preprocess(self):
        # Make available the files for the current run
        for f in self.inlist:
            os.system("cp %s %s 2>/dev/null" % (self[f],f)) # if the files are identical, an error is sent to /dev/null

    def get_name(self):
        return self.index

class SingleJob(BatchStep):
    """
    Single job (in contrast to BatchJob). Basically the same as a single batch step.
    """

    def __init__(self,basedir,workdir,vaspcmd,distributecmd = "cp -Ruf"):
        BatchStep.__init__(self,0)
        self.basedir = basedir
        self.workdir = workdir
        self.vaspcmd = vaspcmd
        self.distributecmd = distributecmd
        self['OUTCAR'] = 'OUTCAR'
        self['vasprun.xml'] = 'vasprun.xml'
    
    def start(self, analyze_only = False):
        if (not analyze_only): 
            self.execute()


class BatchJobDataExtractor:
    """
    Generic class to extract data from a batch job using the data extraction functions in the VasprunParser or OutcarParser classes.
    The data is ordered in a dataset and can be sorted by any of the columns.
    
    Example 1:
    >>> extractor = BatchJobDataExtractor('./', 'data_spec' = [
    >>>     ['get_final_atom_position',1],
    >>>     ['get_force_on_atom',1]
    >>> ])

    Example 2:
    >>> data_spec = [
    >>>        ['get_incar_property',parameter],
    >>>        ['get_total_energy'],
    >>>        ['get_cpu_time']
    >>>     ]
    >>> extractor = BatchJobDataExtractor(directory = './', use_vasprun_xml = True, data_spec = data_spec )
    >>> # Sort by column 1 (ENCUT):
    >>> extractor.sort(self,0)

    """

    outcar_pattern = 'OUTCAR.([0-9]*)'
    vasprun_pattern = 'vasprun([0-9]*).xml'
    cache_filename = 'oppvasp_cache.csv'

    def __init__(self, directory = '', verbose = False, use_vasprun_xml = False, use_cache = True, data_spec = []):
        """
        Generates a dataset with each entry in <data_spec> making up a column.

        :Parameters:
            directory : string
                directory to look for files in
            verbose : bool
                Print info about each step to screen 
            use_vasprun_xml: bool
                True to scan for vasprun.xml files, False to scan for OUTCAR files. 
            use_cache: bool
                Cache as csv file
            data_spec: array
                Specification of the data columns to extract
        """
        self.directory = directory
        self.verbose = verbose
        self.data_spec = data_spec 
        self.use_xml = use_vasprun_xml

        # Import data
        cache_file = directory + BatchJobDataExtractor.cache_filename
        if use_cache and os.path.isfile(cache_file):
            self.import_csv(cache_file)
        else:
            self.extract_data()
        
        # Write cache
        #if use_cache and not os.path.isfile(cache_file):
        #    self.export_csv(cache_file)
    
    def sort(self, column):
        """
        Sorts the dataset according to the values of the column <column> of the dataset. 
        The columns are just the values specified for <data_spect> in the constructor.
        The data in the sort column must be of a type that can be sorted with numpy.
        """
        tmp = np.array(self.get_data_column(column))
        idx = np.argsort(tmp)
        self.data = [self.data[i] for i in idx]


    def extract_data(self):
        """
        Finds and parses either
        (1) all files matching the filename pattern defined in 
            `BatchJobDataExtractor.vasprun_pattern` as vasprun.xml files or
        (2) all files matching the filename pattern defined in 
            `BatchJobDataExtractor.outcar_pattern` as OUTCAR files
        """
        
        filetypes = (
            ('vasprun',BatchJobDataExtractor.vasprun_pattern),
            ('outcar',BatchJobDataExtractor.outcar_pattern)
        )
        jobs = {}
        for filetype,pattern in filetypes:
            filelist = glob.glob(self.directory + pattern.replace('([0-9]*)','*'))
            for filename in filelist:
                m = re.search(pattern,filename)
                if m:
                    index = m.group(1)
                    if not index in jobs:
                        jobs[index] = { }
                    jobs[index][filetype] = filename
        print "Found output files from %d jobs" % (len(jobs))
        if self.verbose:
            for key,val in jobs.items():
                print key+':'
                for key2,val2 in val.items():
                    print "    "+val2

        num_cols = len(self.data_spec)
        num_rows = len(jobs)
        self.data = [[0 for i in range(num_cols)] for j in range(num_rows)]
        i = 0
        for index, files in sorted(jobs.items()):
            if self.verbose:
                print 
            parser = None
            for j in range(num_cols):
                fn = self.data_spec[j]['fn']
                if not fn in dir(parser):
                    oldparser = parser
                    parser = self._get_parser(files, fn)
                    if parser is None:
                        print "Uh oh, the extractor method '%s' was not found " % (fn) \
                            + "in any of the available parsers for the files:"
                        print files
                        sys.exit(1)
                    if self.verbose:
                        print "switched from %s to %s for %s " % (oldparser,parser,fn)
                extract_function = getattr(parser, self.data_spec[j]['fn'])
                if 'params' in self.data_spec[j]:
                    self.data[i][j] = extract_function(*self.data_spec[j]['params']) 
                else:
                    self.data[i][j] = extract_function() 
                if 'type' in self.data_spec[j]:
                    if self.data_spec[j]['type'] == 'float':
                        self.data[i][j] = float(self.data[i][j])
            i+=1
        if self.verbose:
            print 

    def _get_parser(self, files, method):
        """
        Returns a suitable parser to extract 'method' given a list of 'files'.
        If no suitable parser is found, the method returns None.
        """
        if 'vasprun' in files and method in dir(IterativeVasprunParser):
            return IterativeVasprunParser(files['vasprun'], verbose = self.verbose)
        elif 'vasprun' in files and method in dir(VasprunParser):
            return VasprunParser(files['vasprun'], verbose = self.verbose)
        elif 'outcar' in files and method in dir(OutcarParser):
            return OutcarParser(files['outcar'], verbose = self.verbose)
        return None

    def get_data(self):
        """Returns data"""
        return self.data

    def get_data_column(self, column):
        """Returns one of the columns in the dataset"""
        try:
            column = int(column)
            return [row[column] for row in self.data]
        except ValueError:
            for i in range(len(self.data_spec)):
                if 'name' in self.data_spec[i] and self.data_spec[i]['name'] == column:
                    return [row[i] for row in self.data]
            return ValueError

    def get_data_length(self):
        """Returns length of data array"""
        return len(self.data)

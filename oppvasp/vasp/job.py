"""
The VASP job module contains the SingleJob and BatchJob classes for running
single and batch VASP jobs respectively, and the BatchJobDataExtractor class
for extracting results from a batch job. 

The main intention of the BatchJob class is to act as a general superclass,
to be extended by other classes that adds more specific functionality for 
a given job case. Such subclasses include the ManualBatchJob class,
and the ConvergenceTest class.
"""
import os,shutil,sys,re,datetime
import numpy as np
import glob # for finding files using wildcards
from operator import itemgetter
from time import strftime
from parsers import VasprunParser, IterativeVasprunParser, OutcarParser

__docformat__ = "restructuredtext en"


class BatchJob(object):
    """
    Basic rapper class for running VASP, containing a set of base methods inherited by classes such as
    ManualBatchJob, SingleJob, etc. The BatchJob class is not intended for standalone use.
    """

    def __init__(self, basedir = '', workdir = '', vaspcmd = 'vasp', distributecmd = 'cp -Ruf', verbose = False, summaryfile = 'summary.out'):
        """
        Create a BatchJob object.

        Parameters
        ----------
        basedir: str
            Directory to get files from. 
            Default is current directory.
        workdir: str
            Directory to run VASP in (and to copy files to and forth)
            Default is to leave the files were they are (basedir).
        vaspcmd: str
            Command used to run VASP.
            Defaults to 'vasp'
        distributecmd: str
            Command used to copy files to workdir.
            Default is "Cp -Ruf"
        verbose: Boolean
            Be more verbose
            Default is False
        summaryfile: str
            Filename of summary file. Leave blank to skip writing a summary file. 
            Default is 'summary.out'.

        Examples
        ----------

        >>> from oppvasp.vasp.job import ManualBatchJob
        >>> job = ManualBatchJob(workdir = '/work/0001/', vaspcmd = 'vasp.5.21', verbose = True)
        >>> job.start(dry_run = false)

        """
        print ">>> BatchJob initiated at %s" % strftime("%Y-%m-%d %H:%M:%S")
        self.verbose = verbose
        self.basedir = basedir
        if self.basedir == '':
            self.basedir = os.path.abspath(os.path.curdir)
        self.workdir = workdir
        self.vaspcmd = vaspcmd
        self.distributecmd = distributecmd
        self.summaryfile = summaryfile
        self.steps = []
        self.paramName = 'Run' # summary file header

    def add_step(self, step):
        """
        Adds a step (a single VASP run) to the job.

        Parameters
        ----------
        steps: BatchStep
            Objects of type BatchStep
        """
        self.steps.append(step)
        step.basedir = self.basedir
        step.workdir = self.workdir
        step.vaspcmd = self.vaspcmd
        step.distributecmd = self.distributecmd
        step.verbose = self.verbose
    
    def print_info(self):
        """ Deprecated! Just print the object instead """
        print self

    def __repr__(self):
        """
        Print out some info about the job.
        """
        s = "<--------------------------[ BatchJob ]------------------------------>\n"
        s += "Basedir: %s\nWorkdir: %s" % (self.basedir,self.workdir)
        for step in self.steps: 
            s += "[Step %d of %d]\n" % (step.index, len(self.steps))
            s += step.preprocess_info()
            prepfiles = []
            for f in step.input_files.keys():
                if f != step[f]:
                    prepfiles.append("%s -> %s" % (step[f],f))
            if len(prepfiles) != 0:
                s += "  -> Copy %s\n" % ', '.join(prepfiles)
            s += "  -> %s\n" % (step.vaspcmd)
            prepfiles = []
            for f in step.outlist:
                if f != step[f]:
                    prepfiles.append("%s -> %s" % (f, step[f]))
            if len(prepfiles) != 0:
                s += "  -> Rename %s\n" % ', '.join(prepfiles)
            for task in step.post_processing_tasks:
                s += "  -> %s\n" % (task['desc'])
        s += "<-------------------------------------------------------------------->"
        return s


    def start(self, dry_run = False, first_step = 1):
        """
        Starts the job.

        Parameters
        ----------
        dry_run: Boolean 
            Go through the job without actually running VASP. 
            Useful for testing and for running post-analysis steps, if any.
            Default is False
        first_step: int
            Step number to start at. Set this if you want to continue a job.
            Default is 1 (the first step).
        """

        print "[DEBUG] Chdir to %s" % (self.basedir)
        os.chdir(self.basedir)

        # Initialize summary file
        if self.summaryfile != '' and not os.path.isfile(self.summaryfile):
            sf = open(self.summaryfile,'w')
            sf.write('# File: '+os.path.abspath(self.summaryfile)+'\n')
            sf.write('# Created: '+datetime.datetime.today().strftime("%F %T")+'\n')
            sf.write(str(self.paramName)+"\t# pts\tDist.\tToten\t\tCPU\tForce\tPress.\tDrift\n")
            sf.close()
        
        # Loop through the batch steps 
        for step in self.steps:
            print ">>> BatchJob step %i started at %s" % (step.index, strftime("%Y-%m-%d %H:%M:%S"))

            # Run BatchStep.execute
            if not dry_run and step.index >= first_step:
                step.execute()

            # Then update summaryfile
            if self.summaryfile != '' and os.path.isfile(self.summaryfile):
                self.update_summaryfile(step)

        print ">>> BatchJob complete at %s" % strftime("%Y-%m-%d %H:%M:%S")

    def update_summaryfile(self, step):
        if not os.path.isfile(step['vasprun.xml']):
            print "No output file '%s' to analyze" % (step['vasprun.xml'])
            return
        try:
            vasprun = VasprunParser(step['vasprun.xml'])
        except:
            print "BatchJob: Failed to parse vasprun.xml. Did VASP crash?"
            sys.exit(1)
        
        #try:
        #    pressure = "%.4f" % outcar.get_max_pressure()
        #except AttributeError:
        #    pressure = "  -  "

        final_struct = vasprun.get_final_structure()
        shortest_bond = final_struct.get_shortest_bond()
        forces = final_struct.get_forces('c')
        if forces != None:
            fx = np.sum(forces[:,0])
            fy = np.sum(forces[:,1])
            fz = np.sum(forces[:,2])
            maxforce = np.max( np.sum(forces**2,axis=1) )
            drift = (fz,fy,fz)
            print drift
        else:
            print "forces not found"
            maxforce = -1

        try:
            pressure = vasprun.get_final_pressure()
        except:
            print "ERROR: Could not read pressure from vasprun"
            pressure = -1

        try:
            maxdrift = vasprun.get_maxdrift()
        except:
            print "ERROR: Could not read max drift"
            maxdrift = -1 

        cputime,realtime = vasprun.get_time_spent()
        print cputime

        summaryline = "%s\t%d\t%.3f\t%.4f\t%.0f\t%.4f\t%s\t%.4f" % (
            step.get_name(),
            vasprun.get_num_kpoints(),
            shortest_bond[0],
            vasprun.get_total_energy(),
            cputime,
            maxforce,
            pressure,
            maxdrift
        )
        sf = open(self.summaryfile,'a')
        sf.write(summaryline+"\n")
        sf.close()

    def keep(self, filename):
        """
        This will add <filename> to the list of files to be renamed 
        after a step has been carried out
        """
        BatchStep.output_files[filename] = filename+'#'
        for step in self.steps:
            step._re_check_output_files()


class SingleJob(BatchJob):
    """
    Single job (in contrast to BatchJob). Basically the same as a single BatchStep.
    """

    def __init__(self, basedir = '', workdir = '', vaspcmd = 'vasp', distributecmd = 'cp -Ruf', verbose = False, summaryfile = ''):
        """
        Create a SingleJob object.

        Parameters
        ----------
        basedir: str
            Directory to get files from. 
            Default is current directory.
        workdir: str
            Directory to run VASP in (and to copy files to and forth)
            Default is to leave the files were they are (basedir).
        vaspcmd: str
            Command used to run VASP.
            Defaults to 'vasp'
        distributecmd: str
            Command used to copy files to workdir.
            Default is "Cp -Ruf"
        verbose: Boolean
            Be more verbose
            Default is False
        """
        BatchJob.__init__(self, basedir = basedir, workdir = workdir, vaspcmd = vaspcmd, distributecmd = distributecmd, verbose = verbose)

        step = BatchStep(1)
        
        # loop over input files (INCAR, KPOINTS, ...)
        for template_name in BatchStep.input_files.keys():
            indexedname = step[template_name]
            if os.path.exists(indexedname):
                # a unique file was found. That makes this step necessary.
                pass
            elif os.path.exists(template_name):
                # use index-less 'template' file
                step[template_name] = template_name 
            else:
                raise Exception("Neither a file '%s' nor a file '%s' were found. I'm not sure how to deal with this situation." % (template_name, indexedname))
           
        # and over output files (OUTCAR, vasprun.xml, ..)
        for template_name in BatchStep.output_files.keys():
            step[template_name] = template_name 

        # Don't rename any output files
        step.outlist = {}

        self.add_step(step)


class ManualBatchJob(BatchJob):
    """
    Class for carrying out a `manual' batch job, that is a batch job in which a set of 
    input files is prepared manually (without the aid of this package) and named 
    '[FILENAME]_[INDEX]', where [FILENAME] is any input file name,
    such as 'INCAR', 'POSCAR', 'POTCAR' or 'KPOINTS', and [INDEX] is the integer index
    of the job step (starting at 1). [INDEX] can be prefixed by any number of zeros, 
    so the files 'INCAR_1' and 'INCAR_001' are considered equally.
    
    Examples
    ----------
    Preparing the files, 
        INCAR_1, INCAR_2, INCAR_3, KPOINTS_1, KPOINTS_2, POSCAR, POTCAR,
    will generate the job:
        Step    Files
        1       INCAR_1 KPOINTS_1 POSCAR POTCAR
        2       INCAR_2 KPOINTS_2 POSCAR POTCAR
        3       INCAR_3 KPOINTS_2 POSCAR POTCAR
    """

    def __init__(self, **kwargs):
        """
        Create a ManualBatchJob object. 

        Parameters
        ----------
        All parameters are passed over to BatchJob.
        Please see BatchJob for help on the parameters.
        """
        BatchJob.__init__(self, **kwargs)

        os.chdir(self.basedir)
        
        # Add steps to the job with increasing index, 
        # until no files with the given index are found 
        index = 1
        while True: # bad habit...
            step = ManualBatchStep(index)
            if index == 1:
                step_necessary = True     # include at least one step 
            else:
                step_necessary = False    # be pessimistic 
            # loop over input files (INCAR, KPOINTS, ...)
            for template_name in BatchStep.input_files.keys():
                indexedname = step[template_name]
                if os.path.exists(indexedname):
                    # a unique file was found. That makes this step necessary.
                    step_necessary = True
                elif len(self.steps) > 0:
                    # use the file from the previous step
                    previous_step_file = self.steps[len(self.steps)-1][template_name]
                    step[template_name] = previous_step_file
                elif os.path.exists(template_name):
                    # use index-less 'template' file
                    step[template_name] = template_name
                else:
                    raise Exception("Neither a file '%s' nor a file '%s' were found. I'm not sure how to deal with this situation." % (template_name, indexedname))
           
            if not step_necessary:
                break

            self.add_step(step)
            index += 1
        if self.verbose:
            self.print_info()


class BatchStep(object):
    """
    This class defines a single VASP run, typically being one step in a BatchJob.
    
    Examples:
    ----------

    >>> from oppvasp.vasp.job import Job, ManualBatchJob
    >>> jobs = ManualBatchJob()
    >>> jobs.add_step(BatchStep())
    
    To include CHGCAR in the list of files to be saved:
    ('#' is replaced by the '_i' where i is the index of the job)
    >>> step = BatchStep(...)
    >>> step.keep('CHGCAR')

    """

    input_files = {
            'INCAR': 'INCAR#', 
            'POSCAR': 'POSCAR#', 
            'POTCAR': 'POTCAR#',
            'KPOINTS': 'KPOINTS#'
        }
    output_files = {
            'OUTCAR': 'OUTCAR#',
            'vasprun.xml': 'vasprun#.xml',
            'XDATCAR': 'XDATCAR#'
        } 

    def __init__(self, index):
        """
        Create a BatchStep object. 

        Parameters
        ----------
        index: int
            The index of the batch step. This should be an increasing number, 
            starting at 1 for the first step.
            Index 0 is given is this is a single job, not a batch job.
        """

        self.index = index
        self.files = {}
        self.post_processing_tasks = []

        filling = None
        for template_name in BatchStep.input_files.keys():
            if index != 0:
                f = self.re_matchfile('.',BatchStep.input_files[template_name],index)
                if f != False:
                    filling = f
        if filling == None:
            filling = '_'

        for template_name in BatchStep.input_files.keys():
            if index == 0:
                indexed_name = template_name
            else:
                indexed_name = BatchStep.input_files[template_name].replace('#','%s%d' % (filling,index))
            #print "    [%s] => %s" % (template_name,indexed_name)
            self.files[template_name] = indexed_name

        self._re_check_output_files()

    def _re_check_output_files(self):

        # necessary? no, not really
        self.outlist = BatchStep.output_files.keys()
        
        filling = None
        for template_name in BatchStep.input_files.keys():
            if self.index != 0:
                f = self.re_matchfile('.',BatchStep.input_files[template_name],self.index)
                if f != False:
                    filling = f
        if filling == None:
            filling = '_'

        for template_name in BatchStep.output_files.keys():
            if self.index == 0:
                indexed_name = template_name
            else:
                indexed_name = BatchStep.output_files[template_name].replace('#','%s%d' % (filling,self.index))
            self.files[template_name] = indexed_name

    def re_matchfile(self,dir,pattern,index):
        names = os.listdir(dir)  # we are in basedir
        re_pattern = '('+pattern.replace('#',')([_.0]*)%d(' % (index))+')'
        rs = re.compile(re_pattern)
        for name in names:
            m = rs.match(name)
            if m:
                basename = m.group(1)
                filling = m.group(2)
                ending = m.group(3)
                return filling
        return False

    
    def __setitem__(self, key, val):
        if not key in (BatchStep.input_files.keys() + BatchStep.output_files.keys()):
            print "%s is not a valid input/output filename" % (key)
            sys.exit(1) # perhaps a bit strict
        self.files[key] = val

    def __getitem__(self, key):
        if not key in (BatchStep.input_files.keys() + BatchStep.output_files.keys()):
            print "%s is not a valid input/output filename" % (key)
            print 
            sys.exit(1) # perhaps a bit strict 
        return self.files[key]

    def __str__(self):
        return " ".join(self.files.values())
    
    def preprocess_info(self):
        """
        This method is called from the print_info method
        """
        return ''

    def preprocess(self):
        """
        This method is called before VASP is executed.
        """
        pass

    def postprocess(self):
        """
        This method is called after VASP is executed.
        """
        pass

    def get_name(self):
        return str(self.index)

    def execute(self): 
        """
        This method performs the job defined by this BatchStep object:
        1) carry out preprocessing tasks (if any)
        2) distribute files to workdir (if necessary)
        3) run VASP
        4) copy files back to workdir (if necessary)
        5) rename output files we want to save
        6) carry out postprocessing tasks (if any)
        """

        # Prepare for the current run
        print "[DEBUG] Preprocess" 
        self.preprocess()

        # (Re-)Distribute files to all nodes:
        if self.workdir != '':
            cmd = '%s %s/* %s/' % (self.distributecmd,self.basedir, self.workdir)
            print "[DEBUG] Distribute files: "+cmd
            os.system(cmd)

        # Move to workdir:
        if self.workdir != '':
            print "[DEBUG] Chdir to %s" % (self.workdir)
            os.chdir(self.workdir)

        # Run VASP:
        print "[DEBUG] Run vasp as: %s" % (self.vaspcmd)
        exit_code = os.system(self.vaspcmd)
        if exit_code != 0:
            print "VASP excited with exit code",exit_code
            sys.exit(1)
        
        # Move back to basedir:
        print "[DEBUG] Chdir to %s" % (self.basedir)
        os.chdir(self.basedir)

        # Move files back
        if self.workdir != '':
            cmd = 'cp -Rupf %s/* %s/ ' % (self.workdir, self.basedir)
            print "[DEBUG] Move files back: "+cmd
            os.system(cmd)

        # Rename output files we want to save
        for f in self.outlist:
            if os.path.isfile(f):
                if self.verbose:
                    print "[DEBUG] Renaming %s to %s" % (f, self[f])
                os.rename(f,self[f])

        # More postprocess tasks, if any
        self.postprocess()
        for task in self.post_processing_tasks:
            task['func'](self)

    def add_post_processing_task(self, function, description):
        """
        Adds a function to be called after the step has been carried out
        A reference to the step object will be sent to the function as an argument.
        """
        self.post_processing_tasks.append({'func':function,'desc':description})


class ManualBatchStep(BatchStep):
    """
    A manual batch step that can be added to a BatchJob. 
    This class inherits methods from the BatchStep class.
    """
    
    def __init__(self, index):
        """
        Create a ManualBatchStep object with a given index. 
        Search for input files with the given index.
        
        Parameters
        ----------
        index: int
            The index of the batch step. This should be an increasing number, 
            starting at 1 for the first step.
            Index 0 is given is this is a single job, not a batch job.
        """
        BatchStep.__init__(self, index)

    def preprocess_info(self):
        if not os.path.isfile(self['POSCAR']):
            if self.index == 1:
                return '  -> Copy POSCAR -> %s\n' % self['POSCAR']
            else:
                return '  -> Copy CONTCAR -> %s\n' % self['POSCAR']
        return ''

    def preprocess(self):
        """
        Copies the input files for the current BatchStep object (like 'INCAR.1') 
        into files with correct filenames for use with VASP (like 'INCAR').
        """
        if not os.path.isfile(self['POSCAR']):
            if self.index == 1:
                shutil.copy2('POSCAR',self['POSCAR'])
                print "Copying POSCAR to %s" % self['POSCAR']
            else:
                shutil.copy2('CONTCAR',self['POSCAR'])
                print "Copying CONTCAR to %s" % self['POSCAR']
        for f in BatchStep.input_files.keys():
            if self[f] != f:
                print "[DEBUG] Copy %s -> %s" % (self[f], f)
                shutil.copy2(self[f],f)
            #os.system("cp %s %s 2>/dev/null" % (self[f],f)) # if the files are identical, an error is sent to /dev/null
    def postprocess(self):
        """
        Cleans up
        """
        for f in BatchStep.input_files.keys():
            if f != self[f]:
                print "[DEBUG] Unlink %s" % (f)
                os.unlink(f) # if we copy, say INCAR.1 to INCAR before a run, remove INCAR afterwards

    def get_name(self):
        return str(self.index)


class BatchJobDataExtractor(object):
    """
    Generic class to extract data from a batch job using the data extraction 
    functions in the VasprunParser or OutcarParser classes.
    The data is ordered in a dataset and can be sorted by any of the columns.

    Examples
    ----------
    
    >>> extractor = BatchJobDataExtractor('./', 'data_spec' = [
    >>>     ['get_final_atom_position',1],
    >>>     ['get_force_on_atom',1]
    >>> ])

    >>> data_spec = [
    >>>        ['get_incar_property',parameter],
    >>>        ['get_total_energy'],
    >>>        ['get_cpu_time']
    >>>     ]
    >>> extractor = BatchJobDataExtractor(directory = './', use_vasprun_xml = True, data_spec = data_spec )
    >>> extractor.sort(self,0) # Sort by column 1 (ENCUT)

    """

    # If changing, remember to change the search pattern in glob.glob ~ 55 lines below too 
    outcar_pattern = 'OUTCAR.*?([0-9]+)'
    vasprun_pattern = 'vasprun.*?([0-9]+).xml'

    def __init__(self, directory = '', verbose = False, data_spec = []):
        """
        Generates a dataset with each entry in <data_spec> making up a column.

        Parameters
        ----------
            directory : string
                directory to look for files in
            verbose : bool
                Print info about each step to screen 
            data_spec : array
                Specification of the data columns to extract

        """
        self.directory = directory
        self.verbose = verbose
        self.data_spec = data_spec 

        # Import data
        self.extract_data()
        
    
    def sort(self, cname):
        """
        Sorts the dataset according to the values of the column <column> of the dataset. 
        The columns are just the values specified for <data_spect> in the constructor.
        The data in the sort column must be of a type that can be sorted with numpy.
        """
        
        args = np.argsort([np.linalg.norm(a) for a in self.data[cname]]) # requires numerical values!! what bout strings?
        for k in self.data:
            self.data[k] = self.data[k][args]
        #tmp = np.array(self.get_data_column(column))
        #idx = np.argsort(tmp)
        #self.data = [self.data[i] for i in idx]


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
            #print self.directory + pattern.replace('([0-9]*)','*')
            filelist = glob.glob(self.directory + pattern.replace('.*?([0-9]+)','*'))
            for filename in filelist:
                #print pattern,' :: ',filename
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
        self.num_rows = num_rows
        self.data = {}
        
        i = 0
        for index, files in sorted(jobs.items()):
            if self.verbose:
                print 
            parser = None
            for spec in self.data_spec:
                fn = spec['fn']
                sname = spec['name']
                if not fn in dir(parser):
                    oldparser = parser
                    parser = self._get_parser(files, fn)
                    #print "Parser is",str(parser)
                    if parser is None:
                        print "Uh oh, the extractor method '%s' was not found " % (fn) \
                            + "in any of the available parsers for the files:"
                        print files
                        sys.exit(1)
                    if self.verbose:
                        print "switched from %s to %s for %s " % (oldparser,parser,fn)
                extract_function = getattr(parser, fn)
                if 'params' in spec:
                    res = extract_function(*spec['params']) 
                else:
                    res = extract_function()
                if not sname in self.data:
                    # initialize array with correct dimension:
                    if type(res) == float or type(res) == np.float64:
                        shape = (num_rows)
                    elif type(res) == np.ndarray:
                        shape = (num_rows, len(res))
                    else:
                        raise StandardError("Unknown type %s returned" % str(type(res)))
                    self.data[spec['name']] = np.zeros(shape)
                
                #print res
                self.data[sname][i] = res
                #if 'type' in spec:
                #    if self.data_spec[j]['type'] == 'float':
                #        self.data[i][j] = float(self.data[i][j])
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
        return (self.num_rows)


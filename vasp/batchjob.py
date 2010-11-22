#############################################################################
#
# @file batchjob.py @version 2
# This file should be called by <jobfile.sh>
# Last modified: Nov 22, 2010 22:21:10
#
# Example usage:
#
#   import os
#   from oppvasp.vasp.batchjob import ManualBatchJob 
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
#   job = ManualBatchJob(basedir,workdir,vaspcmd)
#   job.start(analyzeOnly)
#
#############################################################################
import os,shutil,sys,re
from parsers import poscarParser, outcarParser

class BatchJob:

    def __init__(self,basedir,workdir,vaspcmd):
        self.basedir = basedir
        self.workdir = workdir
        self.vaspcmd = vaspcmd
        self.summaryfile = 'summary.out'
        self.steps = []
        self.paramName = 'Run' # summary file header

    def addStep(self,step):
        self.steps.append(step)
    
    def info(self):
        print "A total number of %d runs will be carried out." % (len(self.steps))
        print " Index\t Step info"
        for step in self.steps: 
            print " %4d\t %s" % (step.index,step)
        print


    def start(self,analyzeOnly = False):

        os.chdir(self.basedir)

        # Initialize summary file
        if not os.path.isfile(self.summaryfile):
            sf = open(self.summaryfile,'w')
            sf.write(self.paramName+"\t# pts\tDist.\tToten\t\tCPU\tForce\tPress.\tDrift\n")
            sf.close()

        # (4) Go!
        for step in self.steps:
            
            if not analyzeOnly:
                # Prepare for the current run
                step.preProcess()
                
                # (Re-)Distribute files to all nodes:
                os.system("srun --ntasks=$SLURM_JOB_NUM_NODES cp -Rupf %s/* %s/" % (self.basedir, self.workdir))
                
                # Run VASP:
                os.chdir(self.workdir)
                os.system(self.vaspcmd)
                os.chdir(self.basedir)
                
                # Move files back
                os.system('cp -Rupf %s/* %s/ ' % (self.workdir, self.basedir))

                # Save some output files
                os.rename('OUTCAR',step['OUTCAR'])
                os.rename('vasprun.xml',step['vasprun.xml'])
                
                # May be used to save more files
                step.postProcess()

            # Analyze output and print results to summary file

            #vasprun = vasprunParser('vasprun.%d.xml' % (step.index))
            if not os.path.isfile(step['OUTCAR']):
                print "No output file '%s' to analyze. Exciting..." % (step['OUTCAR'])
                sys.exit(1)

            poscar = poscarParser(step['POSCAR'])
            outcar = outcarParser(step['OUTCAR'], poscar.selective)

            summaryline = "%s\t%d\t%.3f\t%.4f\t%.0f\t%.4f\t%.4f\t%.4f" % (
                step.getName(),
                outcar.kpoints,
                outcar.dist,
                outcar.toten,
                outcar.cpu,
                outcar.maxforce,
                outcar.maxpressure,
                outcar.maxdrift)
            #print summaryline
            sf = open(self.summaryfile,'a')
            sf.write(summaryline+"\n")
            sf.close()



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

class BatchStep:

    def __init__(self,index):

        self.index = index
        self.inlist = ['INCAR','POSCAR','POTCAR','KPOINTS']
        self.outlist = ['OUTCAR','vasprun.xml']
        # Set default values
        self.files = {}
        for i in self.inlist:
            self.files[i] = i
        self.files['OUTCAR'] = 'OUTCAR.%d' % (self.index)
        self.files['vasprun.xml'] = 'vasprun%d.xml' % (self.index)
    
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

    def preProcess(self):
        # Make available the files for the current run
        for f in self.inlist:
            os.system("cp %s %s 2>/dev/null" % (self[f],f)) # if the files are identical, an error is sent to /dev/null

    def postProcess(self):
        pass

    def getName(self):
        return self.index




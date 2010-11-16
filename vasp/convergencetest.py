#!/usr/bin/env python
#############################################################################
#
# @file convergencetest.py @version 2
# This file should be called by <jobfile.sh>
# Last modified: Nov 16, 2010 22:19:48
#
#############################################################################
import os,shutil,sys,re,math

import vasprunparser

inp = ['INCAR','POSCAR','POTCAR','KPOINTS']

class ConvergenceTest:

    def __init__(self,basedir,workdir,vaspcmd):
        self.basedir = basedir
        self.workdir = workdir
        self.vaspcmd = vaspcmd
        
        # Read parameter input
        os.chdir(self.basedir)
        if not os.path.isfile('convergencetest.in'):
            print "convergencetest.in not found!"
            sys.exit(1)
        f = open('convergencetest.in','r')
        lines = f.readlines()
        f.close()
        self.param = lines[0].strip() #trim
        self.paramValues = [l.strip() for l in lines[1:]]
        self.numSteps = len(self.paramValues)
        print "A total number of %d runs will be carried out." % (self.numSteps)

    def start(self,analyzeOnly = False):

        os.chdir(self.basedir)
        
        # Create results dir and initialize summary file:
        summaryfilename = 'summary.out'
        if not os.path.isfile(summaryfilename):  # does file exist?
            summaryfile = open(summaryfilename,'w')
            summaryfile.write(self.param+"\t# pts\tDist.\tToten\t\tCPU\tForce\tPress.\tDrift\n")
            summaryfile.close()

        # Loop 
        for stepNo in range(self.numSteps):
            
            # Update INCAR for the current run:
            f = open('INCAR','r+'); fc = f.read()
            fcr = re.sub(
                r'%s([ \t]*)=([ \t]*)([.\w]*)' % (self.param),
                '%s\\1=\\2 %s' % (self.param,self.paramValues[stepNo]),
                fc)
            if fc == fcr:
                fcr = fc+"\n %s = %s\n" % (self.param,self.paramValues[stepNo])
            f.seek(0); f.write(fcr); f.truncate(); f.close();
            
            if not analyzeOnly:
                
                # (Re-)Distribute files to all nodes:
                os.system("srun --ntasks=$SLURM_JOB_NUM_NODES cp -Rupf %s/* %s/" % (self.basedir, self.workdir))
                
                # Run VASP:
                os.chdir(self.workdir)
                os.system(self.vaspcmd)
                os.chdir(self.basedir)
                
                # Move files back
                os.system('cp -Rupf %s/* %s/ ' % (self.workdir, self.basedir))


            # Analyze output and print results to summary file

            #vasprun = vasprunParser('vasprun.xml')
            poscar = poscarParser('POSCAR')
            outcar = outcarParser('OUTCAR', poscar.selective)

            summaryline = "%s\t%d\t%.3f\t%.4f\t%.0f\t%.4f\t%.4f\t%.4f" % (
                self.paramValues[stepNo],
                outcar.kpoints,
                outcar.dist,
                outcar.toten,
                outcar.cpu,
                outcar.maxforce,
                outcar.maxpressure,
                outcar.maxdrift)
            print summaryline
            summaryfile = open(summaryfilename,'a')
            summaryfile.write(summaryline+"\n")
            summaryfile.close()
            
            # Save the following files
            os.rename('OUTCAR','OUTCAR.%d' % (stepNo))
            os.rename('DOSCAR','DOSCAR.%d' % (stepNo))
            os.rename('vasprun.xml','vasprun.%d.xml' % (stepNo))

        # End of loop

class BatchStep:
    
    def __init__(self,index):
        self.files = {}
        self.index = index
        self.necessary = False
        for f in inp:
            self[f] = f
            if os.path.exists(f) and os.path.exists(f+'.0'):
                print "Both a file '%s' and a file '%s' exist. I'm terribly confused." % (f,f+'.0')
                sys.exit(1)
            elif os.path.exists(f+'.0'):
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

    def __setitem__(self,key,val):
        if not key in inp:
            print "%s is not a valid input/output filename" % (key)
            sys.exit(1)
        self.files[key] = val

    def __getitem__(self,key):
        if not key in inp:
            print "%s is not a valid input/output filename" % (key)
            sys.exit(1)
        return self.files[key]

    def __str__(self):
        return " ".join(self.files.values())




class poscarParser:
    
    def __init__(self, poscarname='POSCAR'):
        self.selective = 0 
        self.readPOSCAR(poscarname)

    def readPOSCAR(self,poscarname):
        poscarfile = open( poscarname, 'r')  # r for reading
        commentline = poscarfile.readline()
        scaleline = poscarfile.readline()
        vec1line = poscarfile.readline()
        vec2line = poscarfile.readline()
        vec3line = poscarfile.readline()
        sixthline = poscarfile.readline()  # Test for vasp5 syntax
        try:
            dummy = int(sixthline.split()[0])
            atomnumberline = sixthline
        except:
            atomnumberline = poscarfile.readline()
        atomnumbers = map(int,atomnumberline.split())
        natoms = 0
        for i in range(len(atomnumbers)):
            natoms = natoms + atomnumbers[i]
        seventhline = poscarfile.readline()

        if seventhline[0] == 'S' or seventhline[0] == 's':
            self.selective = 1
            seventhline = poscarfile.readline()
        else:
            self.selective = 0

        if self.selective:
           x = []; y = []; z = []
           for j in range(natoms):
            line = poscarfile.readline()  # read a line
            if not line: break
            xyz = line.split()
            x.append(xyz[3])
            y.append(xyz[4])
            z.append(xyz[5])
        
        poscarfile.close()

class outcarParser:

    def __init__(self, outcarname = 'OUTCAR', selective = 0):

        self.kpoints = 0
        self.dist = 0
        self.toten = 0
        self.planewaves = 0
        self.cpu = 0
        self.memory = 0
        self.maxforce = 0
        self.maxdrift = 0
        self.maxpressure = 0

        outfile = open(outcarname, 'r')
        while 1:
            line = outfile.readline()
            kmatch = re.search('(\d+) +irreducible', line)
            ematch = re.search('free  energy   TOTEN  = +(-*\d+.\d+)', line)
            cpumatch = re.search('Total CPU time used \(sec\): +(\d+.\d+)', line)
            distmatch = re.search('Following cartesian coordinates:', line)

            #k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=     96
            #number of dos      NEDOS =    301   number of ions     NIONS =      8
            #non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8
            #total plane-waves  NPLWV =  32768

            planewavematch = re.search('NPLWV[ \t]*=[ \t]*([0-9])', line)
            nbandsmatch = re.search('NPLWV[ \t]*=[ \t]*([0-9])', line)
            if planewavematch:
                self.planewaves = int(planewavematch.group(1))
            if kmatch:
                self.kpoints = int(kmatch.group(1))
            elif ematch:
                self.toten= float(ematch.group(1))
            elif cpumatch:
                self.cpu = float(cpumatch.group(1))
            elif distmatch:
                if self.kpoints > 1:
                    tmpline = outfile.readline()
                    firstline = outfile.readline()
                    secondline = outfile.readline()
                    k1x,k1y,k1z,dummy = map(float, firstline.split())
                    k2x,k2y,k2z,dummy = map(float, secondline.split())
                    self.dist = math.sqrt((k2x-k1x)*(k2x-k1x)+(k2y-k1y)*(k2y-k1y)+(k2z-k1z)*(k2z-k1z))
                else:
                    self.dist = 0
            elif re.search(r'external pressure', line): 
                tmp,tmp,tmp,pressure,tmp,tmp,tmp,tmp,tmp,tmp = line.split()
                self.maxpressure = float(pressure)
            elif re.search(r'TOTAL\-FORCE', line):
                i=0
                line = outfile.readline()
                maxdrift = 0.0
                maxforce = 0.0
                while 1:
                    line = outfile.readline()
                    if re.search(r'----', line):
                        line = outfile.readline()
                        a,b,driftx,drifty,driftz = line.split()
                        if abs(float(driftx)) > maxdrift:
                            maxdrift = abs(float(driftx))
                        if abs(float(drifty)) > maxdrift:
                            maxdrift = abs(float(drifty))
                        if abs(float(driftz)) > maxdrift:
                            maxdrift = abs(float(driftz))
                        break
                    posx,posy,posz,forx,fory,forz = map(float, line.split())
                    if selective:
                        if (abs(forx) > maxforce) and (x[i] == 'T' or x[i] == 't'):
                            maxforce = abs(forx)
                            maxi = i
                        if (abs(fory) > maxforce) and (y[i] == 'T' or y[i] == 't'):
                            maxforce = abs(fory)
                            maxi = i
                        if (abs(forz) > maxforce) and (z[i] == 'T' or z[i] == 't'):
                            maxforce = abs(forz)
                            maxi = i
                    else:
                        if abs(forx) > maxforce:
                            maxforce = abs(forx)
                            maxi = i
                        if abs(fory) > maxforce:
                            maxforce = abs(fory)
                            maxi = i
                        if abs(forz) > maxforce:
                            maxforce = abs(forz)
                            maxi = i
                        i = i+1
                self.maxforce = maxforce
                self.maxdrift = maxdrift
            if not line:
                break
        outfile.close()

# Example use:
#if __name__ == '__main__':
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

#    job = ConvergenceTest(basedir,workdir,vaspcmd)
#    job.start(analyzeOnly)


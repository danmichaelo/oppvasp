import os
from batchjob import BatchJob
import oppvasp.utils

class VolumeTestCubicUnitCell(BatchJob):

    def __init__(self,basedir,workdir,vaspcmd):
        BatchJob.__init__(self,basedir,workdir,vaspcmd)

        # Read parameter input
        self.parameterfile = 'volumetest.in'
        os.chdir(self.basedir)
        if not os.path.isfile(self.parameterfile):
            print "Parameter-file '%s' not found!" % (self.parameterfile)
            sys.exit(1)
        f = open(self.parameterfile,'r')
        lines = f.readlines()
        f.close()
        aMin = float(lines[0].strip())
        aMax = float(lines[1].strip())
        aStep = float(lines[2].strip())

        paramValues = list(utils.frange6(aMin,aMax,aStep))
         
        POSCAR = open('POSCAR', 'r')
        plines = POSCAR.readlines()
        POSCAR.close()

        for i in range(len(paramValues)):
            plines[1] = '%.4f\n' % (paramValues[i])
            ifile = open('POSCAR.%d' % (i), 'w')
            ifile.writelines(plines)
            ifile.close()
            self.addStep(ManualBatchStep(['POSCAR'],i))

        self.info()


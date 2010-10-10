'''
Unit test for plotbands.py
Created on 22. mars 2010

@author: danmichael
'''

import abpytho.vasp.parsebandstructure as vpb
import abpytho.espresso.parsebandstructure as epb
import abpytho.plotbands as pb
import unittest

class TestPlotBands(unittest.TestCase):
    
    def setUp(self):
        self.plot = pb.PlotBands({ 'basis': 'reciprocal', 'points': [
             ['\\Gamma', [1.0, 1.0,  1.0]],
             ['X'      , [0.5, 0.5,  0.0]],
             ['L'      , [0.5, 0.5,  0.5]],
             ['W'      , [0.5, 0.75, 0.25]]
         ]})

    def tearDown(self):
        #self.widget.dispose()
        #self.widget = None
        pass

    def testPlotbands(self):
        self.plot.addBands(epb.ParseBandStructure('si.band','si.band.out'), fermiLevel = 6.54, color='blue', label='Espresso', linestyle='dashed')
        self.plot.addBands(vpb.ParseBandStructure('PROCAR','OUTCAR'), fermiLevel = 6.54, color='red', label='VASP')
        self.plot.plotAndSaveToFile(outFile = 'bands.pdf', yMin = -13.0, yMax = 4.0, basis='cartesian')


if __name__ == "__main__":
    unittest.main()

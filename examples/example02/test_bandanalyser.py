'''
Unit test for bandanalyser.py
Created on 22. mars 2010

@author: danmichael
'''

import abpytho.vasp.parsebandstructure as vpb
import abpytho.espresso.parsebandstructure as epb
import abpytho.bandanalyser as ba
import unittest

class TestBandAnalyser(unittest.TestCase):
    
    def setUp(self):
        specialpoints = { 'basis': 'reciprocal', 'points': [
             ['\\Gamma', [1.0, 1.0,  1.0]],
             ['X'      , [0.5, 0.5,  0.0]],
             ['L'      , [0.5, 0.5,  0.5]],
             ['W'      , [0.5, 0.75, 0.25]]
        ]}
        bandstructure = epb.ParseBandStructure('si.band','si.band.out')
        self.analyser = ba.BandAnalyser(bandstructure, specialpoints, fermiLevel = 6.54)

    def tearDown(self):
        pass

    def testBandAnalyser(self):
        self.analyser.analyse()


if __name__ == "__main__":
    unittest.main()

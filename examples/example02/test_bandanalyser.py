'''
Unit test for bandanalyser.py
Created on 22. mars 2010

@author: danmichael
'''

import oppvasp.vasp.parsebandstructure as vpb
import oppvasp.espresso.parsebandstructure as epb
import oppvasp.bandanalyser as ba

specialpoints = { 'basis': 'reciprocal', 'points': [
     ['\\Gamma', [1.0, 1.0,  1.0]],
     ['X'      , [0.5, 0.5,  0.0]],
     ['L'      , [0.5, 0.5,  0.5]],
     ['W'      , [0.5, 0.75, 0.25]]
]}

bandstructure = epb.ParseBandStructure('si.band','si.band.out')
analyser = ba.BandAnalyser(bandstructure, specialpoints, fermiLevel = 6.54)
analyser.analyse()


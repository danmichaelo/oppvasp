'''
Created on 22. mars 2010
@author: danmichael
'''

import oppvasp.vasp.parsebandstructure as vpb
import oppvasp.espresso.parsebandstructure as epb
import oppvasp.plotbands as pb

plot = pb.PlotBands({ 'basis': 'reciprocal', 'points': [
     ['\\Gamma', [1.0, 1.0,  1.0]],
     ['X'      , [0.5, 0.5,  0.0]],
     ['L'      , [0.5, 0.5,  0.5]],
     ['W'      , [0.5, 0.75, 0.25]]
 ]})

plot.addBands(
    epb.ParseBandStructure('si.band','si.band.out'), 
    fermiLevel = 6.54, color='blue', label='Espresso', linestyle='dashed'
    )

plot.addBands(
    vpb.ParseBandStructure(), 
    fermiLevel = 6.54, color='red', label='VASP'
    )

plot.plotAndSaveToFile(
    outFile = 'bands.pdf', basis='cartesian',
    yRange = [-13.0, 4.0], yTicks = [-13.0, 4.0, 1.0] 
)


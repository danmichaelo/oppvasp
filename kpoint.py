# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4
#
# Created on 22. mars 2010
# @author: danmichael

import copy
import numpy as np

class KPoint():

    def __init__(self,x,y = 0,z = 0,eigenvals=[]):
        if type(x).__name__ == 'ndarray':
            if x.ndim != 1 or x.size != 3:
                raise StandardError("Vector must have dimension 1 and size 3")
            self.k = x
        else:
            self.k = np.array([x,y,z])
        self.eigenvals = eigenvals
        
    def appendEigenval(self, eigenval):
        self.eigenvals.append(eigenval)
    
    def getVector(self):
        return copy.deepcopy(self.k)
    
    def setVector(self, k):
        if type(k).__name__ != 'ndarray':
            raise StandardError("Not a vector")
        if k.ndim != 1 or k.size != 3:
            raise StandardError("Vector must have dimension 1 and size 3")
        self.k = k
    
    def __eq__(self,other):
        return np.all(self.k==other.getVector())
    
    def __str__(self):
        return "[%.2f,%.2f,%.2f]" % tuple(self.k)

""" Oslo Python Package for VASP """

from sys import version_info
ver = version_info[0]*10 + version_info[1]
if ver < 26:
    raise Exception("oppvasp requires python version 2.6 or higher")

import numpy as np
import time
from element_data import elements, get_atomic_number_from_symbol
from util import direct_to_cartesian, cartesian_to_direct
from read_trajectory import read_trajectory


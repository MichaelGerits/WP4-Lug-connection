import numpy as np
import math as m

"""
This document includes the loads and force calculations
"""

#resulting force [Px, Py, Pz]
P = np.zeros(3) #if a component is zero, take 10% of the total

#resulting Moment [Mx, My, Mz]
T = np.zeros(3)


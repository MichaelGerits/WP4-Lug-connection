import PartDefinition as PD
import Loads as LD
import math as m
import numpy as np

arr = np.zeros(4)
loads = np.array([1,2,3])
arr[0:3] = loads/2

print(arr)
#beep boop bop 
a=0
b=0

checkResult = (1,0)
print(np.abs(np.array(checkResult) - 1))
a, b += np.abs(np.array(checkResult) - 1) * 0.001
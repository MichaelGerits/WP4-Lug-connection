import math as m
import numpy as np
import PartDefinition as PD
import Loads

# 4.3
hinge = PD.Hinge(sigmaY=250)
Loads.F1 = F1

D1 = 0.001
w = D1
t = 0.001 #m
mValuesA = []
DValuesA = []
wValuesA = []

while t < 0.25:
    D1 = (P[1] + F1) * 1.5 / 4 / hinge.sigmaY / t #Bearing stress (1.5 is MS)
    K_t = -0.05 * w/D1 + 3.05
    w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t + D1 #Tension of net section
    m = t * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
    mValuesA.append(m)
    DValuesA.append(D1)
    wValuesA.append(w)
    t += 0.001


D1 = 0.001
w = D1
t = 0.001 #m
mValuesB = []
DValuesB = []
wValuesB = []

while t < 0.25:
    A_frac = 6 / (D1 * (4/(0.5*w-m.sqrt(2)*0.25 * D1) + 2/(0.5*(w-D1)))) # p18 fig D1.15
    K_bending = 1.2143 * A_frac #Fig D1.15 page 18
    D = Loads.P[0] / 4 / (hinge.K_bending * hinge.sigmaY) / t #Bending
    m = t * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
    mValuesB.append(m)
    DValuesB.append(D1)
    wValuesB.append(w)
    t += 0.001

mMinList = [min(mValuesA), min(mValuesB)]
mMax = max(mMinList)
mMaxIndex = mMinList.index(mMax)
if mMaxIndex = 0:
    mMinIndex = mValuesA.index(mMax)
    t = mMinIndex * 0.001 + 0.001
    D1 = DValuesA[mMinIndex]  # Bearing stress (1.5 is MS)
    w = wValuesA[mMinIndex]
else:
    mMinIndex = mValuesB.index(mMax)
    t = mMinIndex * 0.001 + 0.001
    D1 = DValuesB[mMinIndex]  # Bearing stress (1.5 is MS)
    w = wValuesB[mMinIndex]

print(f"t = {t}", f"hole diameter = {D1}", f"width = {w}")




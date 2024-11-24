import math
import numpy as np
import PartDefinition as PD
import Loads

"""
Below are a list of functions which sequentially calculate the dimensions and the stresses of the hinge
These are made functions such that they can be ran again through itteration.

The result of each calculation will be Saved in the hinge object, such that one reference poijt exists.
"""
#initial definition of the hinge obect
hinge = PD.Hinge(t1=0.001, D1=0.01, w=0.01 ,sigmaY=250)

# 4.3------------------------------------------------------------------------------------------------------------------------------------
Loads.F1 = F1 #TODO what does this mean?
def CalcLugDim(t1_init, D1_init, w_init):
    #TODO: explain how the algorithm works please
    """

    """

    D1 = D1_init
    w = w_init
    t1 = t1_init #m
    mValuesA = []
    DValuesA = []
    wValuesA = []

    while t < 0.25:
        D1 = (P[1] + F1) * 1.5 / 4 / hinge.sigmaY / t1 #Bearing stress (1.5 is MS)
        K_t = -0.05 * w/D1 + 3.05
        w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t1 + D1 #Tension of net section
        m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
        mValuesA.append(m)
        DValuesA.append(D1)
        wValuesA.append(w)
        t1 += 0.001

    #TODO: comment why the same this has been done twice please
    """

    """

    D1 = D1_init
    w = w_init
    t1 = t1_init #m
    mValuesB = []
    DValuesB = []
    wValuesB = []

    while t1 < 0.25:
        A_frac = 6 / (D1 * (4/(0.5*w-m.sqrt(2)*0.25 * D1) + 2/(0.5*(w-D1)))) # p18 fig D1.15
        K_bending = 1.2143 * A_frac #Fig D1.15 page 18
        D = Loads.P[0] / 4 / (hinge.K_bending * hinge.sigmaY) / t #Bending
        m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
        mValuesB.append(m)
        DValuesB.append(D1)
        wValuesB.append(w)
        t1 += 0.001

    #TODO: explain how the optimum is chosen

    mMinList = [min(mValuesA), min(mValuesB)]
    mMax = max(mMinList)
    mMaxIndex = mMinList.index(mMax)
    if mMaxIndex == 0:
        mMinIndex = mValuesA.index(mMax)
        t = mMinIndex * 0.001 + 0.001
        D1 = DValuesA[mMinIndex]  # Bearing stress (1.5 is MS)
        w = wValuesA[mMinIndex]
    else:
        mMinIndex = mValuesB.index(mMax)
        t1 = mMinIndex * 0.001 + 0.001
        D1 = DValuesB[mMinIndex]  # Bearing stress (1.5 is MS)
        w = wValuesB[mMinIndex]

    #updates the new values to the hinge object
    hinge.w = w
    hinge.D1 = D1
    hinge.t1 = t1
    print(f"t = {t1}", f"hole diameter = {D1}", f"width = {w}")
#-------------------------------------------------------------------------------------------------------------------------------------------
#runs the function for the first time
CalcLugDim(hinge.t1, hinge.D1, hinge.w)

#4.4-----------------------------------------------------------------------------------------------------------------
def CalcBasePlateDim(e1Fac=1.5, e2Fac=1.5, holeSepFac=2.5):
    """
    calculates the dimensions of the baseplate with the width and the factors of seperation
    """

    hinge.D2 = hinge.w/(2+ 2 * e1Fac + holeSepFac)
    hinge.e1 = e1Fac*hinge.D2
    hinge.e2 = e2Fac*hinge.D2

    hinge.depth = 2*(2* hinge.e2 + hinge.t1) + hinge.h
#---------------------------------------------------------------------------------------------------------------------------
#runs the function for the first time
CalcBasePlateDim()

#TODO:calculate h
import math
import numpy as np
import PartDefinition as PD
import Loads

hinge = PD.Hinge(sigmaY=250)
Loads.F1 = F1 #TODO what does this mean?
# 4.3------------------------------------------------------------------------------------------------------------------------------------
def CalcLugDim(t1_init = 0.001, D1_init = 0.01):
    #TODO: explain what you're doing here please
    """

    """

    D1 = D1_init
    w = D1
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
    w = D1
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
CalcLugDim()

#4.4-----------------------------------------------------------------------------------------------------------------
def CalcBasePlateDim(e1Fac=1.5, e2Fac=1.5, holeSepFac=2.5):
    """
    calculates the dimensions of the baseplate with the initial hole dimension 
    """

    hinge.D2 = hinge.w/(2+ 2 * e1Fac + holeSepFac)
    hinge.e1 = e1Fac*hinge.D2
    hinge.e2 = e2Fac*hinge.D2

    hinge.depth = 2*(2* hinge.e2 + hinge.t1) + hinge.h
#---------------------------------------------------------------------------------------------------------------------------
#runs the function for the first time
CalcBasePlateDim()
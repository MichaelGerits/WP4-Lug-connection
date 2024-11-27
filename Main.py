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
F1 = Loads.F1 #The force caused by the moment, defined in the Loads docuent
P = Loads.P

def CalcLugDim(t1_init, D1_init, w_init):
    #TODO: explain how the algorithm works please
    """
    The program runs through a range of values for thickness(t), and finds hole diameter(D1) and width(w)
    that will meet the bearing and bending stress (hence the two different for loops. one checks for
    bearing, the other for bending. Then, the code finds a "mass(m)" that isnt really mass, but this value
    is proportional to mass. Thus, this value takes on a minimum when mass is a minimum. Then the code finds

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
        A_frac = 6 / (D1 * (4/(0.5*w-math.sqrt(2)*0.25 * D1) + 2/(0.5*(w-D1)))) # p18 fig D1.15
        K_bending = 1.2143 * A_frac #Fig D1.15 page 18
        D1 = P[0] / 4 / (K_bending * hinge.sigmaY) / t #Bending
        K_t = -0.05 * w / D1 + 3.05
        w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t1 + D1
        m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
        mValuesB.append(m)
        DValuesB.append(D1)
        wValuesB.append(w)
        t1 += 0.001

    #TODO: explain how the optimum is chosen
    dList = []
    mList = []
    for i in range(len(DValuesA)):
        dList.append([DValuesA[i], DValuesB[i], mValuesA[i], mValuesB[i]])
        if dList[i][0] >= dList[i][1]:
            dList[i].append("A")
            mList.append(dList[i][2])
        else:
            dList[i].append("B")
            mList.append(dList[i][3])

    mMin = min(mList)
    mMinIndex = mList.index(mMin)

    if dList[mMinIndex][4] == "A":
        t = mMinIndex * 0.001 + 0.001
        D1 = DValuesA[mMinIndex]  # Bearing stress (1.5 is MS)
        w = wValuesA[mMinIndex]
    else:
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
def CalcBasePlateDim(e1Fac=1.5, e2Fac=1.5, holeSepFac=2.5, fastenerAmount=PD.fastenerAmount):
    """
    calculates the dimensions of the baseplate with the width and the factors of seperation
    """

    hinge.D2 = hinge.w/(fastenerAmount/2 + 2 * e1Fac + holeSepFac*(fastenerAmount/2-1))
    hinge.e1 = e1Fac*hinge.D2
    hinge.e2 = e2Fac*hinge.D2

    hinge.depth = 2*(2* hinge.e2 + hinge.t1) + hinge.h
#---------------------------------------------------------------------------------------------------------------------------
#runs the function for the first time
CalcBasePlateDim()

#TODO:calculate h
#4.5------------------------------------------------------------------------------------------------------------------------

def CalcCG(fastenerAmount = PD.fastenerAmount):
    """
    calculate the centre of gravity of the fasteners
    
    The fasteners are defined first.
    The locacal coordinate system is the centre of the baseplate
    
    an even amount of fasteners is needed
    """
    #You can fil in a custom amount of fasteners

    #initialise the list of fasteners
    Fasteners = list(fastenerAmount)

    #define the list of z positions
    posZs = np.linspace(hinge.w/2 - hinge.e1, -hinge.w/2 + hinge.e1, fastenerAmount/2)

    #define the fasteners and their positions
    for i in range(fastenerAmount):
        #deside on which side the fastener is
        xSign = 1 if i%2==0 else -1

        #calculate positions
        posX = xSign * (hinge.depth/4 + hinge.h/4 + hinge.t1/2)
        posZ = posZs[math.floor(i/2)]

        #add to the list of fasteners
        Fasteners[i] = PD.Fastener(D_h = hinge.D2, xPos=posX, zPos=posZ)


    numSum = 0
    denomSum = 0
    for fast in Fasteners:
        #area of the hole
        A = fast.D_h**2 * math.pi * 0.25
        numSum += A * fast.xPos
        denomSum += A
    
    #calculate cg in X
    cgX = numSum/denomSum

    numSum = 0
    denomSum = 0
    for fast in Fasteners:
        #area of the hole
        A = fast.D_h**2 * math.pi * 0.25
        numSum += A * fast.zPos
        denomSum += A

    #calculate CG in Z
    cgZ = numSum/denomSum

    return(cgX, cgZ)
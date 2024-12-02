import math
import numpy as np
import PartDefinition as PD
import Loads
import itertools

"""
Below are a list of functions which sequentially calculate the dimensions and the stresses of the hinge
These are made functions such that they can be ran again through iteration.

The result of each calculation will be Saved in the hinge object, such that one reference poijt exists.
"""

# 4.3------------------------------------------------------------------------------------------------------------------------------------


#Function that actually does the optimising
def CalcLugDim(hinge):
    #explaination of what this thing does
    """
    The program runs through a range of values for thickness(t1), and finds hole diameter(D1) and width(w)
    that will meet the bearing and bending stress (hence the two different for loops. one checks for
    bearing, the other for bending. Then, the code finds a "mass(m)" that isnt really mass, but this value
    is proportional to mass. Thus, this value takes on a minimum when mass is a minimum.
    The equations show that D1 is inversely proportional to both bearing stress and bending stress, for
    constant w and t1. So, the optimiser finds both values of D1 for a given t, and puts the resulting tuples
    in a list(dList) along with their corresponding masses. An if statement checks which one is bigger, and
    appends either an A or a B. Then, a list(mList) is created with the corresponding m values. i.e. if the
    D1 required for bearing is bigger than the D1 required for bending, then an A is appended, and mList
    contains the m from the first for loop.
    The minimum(mMin) m from mList is used as the optimal value. The index(mMinIndex) of this minimum is
    used to find the optimal t, D1, and w values. The A or B appended earlier helps find out which list to
    take these values from.
    """

    t_step = 0.000001

    # Taking values from Loads and giving them nicer names
    P = Loads.P
    F1 = Loads.F1
    #Declare initial values
    D1 = hinge.D1
    w = hinge.w
    t1 = hinge.t1 #m
    Purple = True

    #Make lists to append bearing values of m, D1, and w to later
    mValuesA, DValuesA, wValuesA = [np.array([])] * 3


    #Solves for Bearing stress for a range of t values
    while Purple == True:
        if t1 < 0.01:
            D1 = (P[1] + F1) * 1.5 / (4 * hinge.sigmaY * t1) #Bearing stress (1.5 is MS)
            K_t = -0.05 * w/D1 + 3.05 #from appendix A tab D1.3, Curve 1 (W/D up to 3)
            w = (P[1] + F1) * 1.5 / (4 * K_t * hinge.sigmaY * t1) + D1 #Tension of net section
            m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
            mValuesA = np.append(mValuesA, m)
            DValuesA = np.append(DValuesA, D1)
            wValuesA = np.append(wValuesA, w)
            t1 += t_step
        else:
            Purple = False


    #Reset initial values (crucial for the cursed K's)
    D1 = hinge.D1
    w = hinge.w
    t1 = hinge.t1 #m
    Blue = True

    #Make lists to append bending values of m, D1, and w to later
    mValuesB, DValuesB, wValuesB = [np.array([])] * 3

    #Solves for bending stress for a range of t values
    while Blue == True:
        if t1 < 0.01:
            if w > D1:
                A_frac = 3*(1/D1) / (8 / (2*w -2*D1 +math.sqrt(2)*D1) + 2 / (w - D1)) # p18 fig D1.15
                K_bending = 1.3 / 1.4 * A_frac  # Fig D1.15 page 18
                D1 = P[0] / (4 * K_bending * hinge.sigmaY * t1)  # Bending
                K_t = -0.05 * w / D1 + 3.05  # from appendix A tab D1.3, Curve 1 (W/D up to 3)
                w = (P[1] + F1) * 1.5 / (4 * K_t * hinge.sigmaY * t1) + D1
                m = t1 * (w ** 2 - D1 ** 2)  # this value is not actually mass, but it is proportional to mass
                mValuesB = np.append(mValuesB, m)
                DValuesB = np.append(DValuesB, D1)
                wValuesB = np.append(wValuesB, w)
                print(w ** 2 - D1 ** 2)
            else:
                A_frac = 1  # p18 fig D1.15
                K_bending = 1.3 / 1.4 * A_frac  # Fig D1.15 page 18
                D1 = 4  # Bending
                K_t = -0.05 * w / D1 + 3.05  # from appendix A tab D1.3, Curve 1 (W/D up to 3)
                w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t1 + D1
                mValuesB = np.append(mValuesB, 31415)
                DValuesB = np.append(DValuesB, 0)
                wValuesB = np.append(wValuesB, 31415)
                print("I hate this")
            t1 += t_step
        else:
            Blue = False

    dTest = (DValuesA >= DValuesB)
    mList = [a if c==True else b for a,b,c in zip(mValuesA, mValuesB, dTest)]
    mMinIndex = np.argmin(mList)
    print(mList[761], mList[762], mList[763], mList[764])

    if dTest[mMinIndex] == True:
        t1 = mMinIndex * t_step + hinge.t1
        D1 = DValuesA[mMinIndex]
        w = wValuesA[mMinIndex]
        m = mValuesA[mMinIndex]
    else:
        t1 = mMinIndex * t_step + hinge.t1
        D1 = DValuesB[mMinIndex]
        w = wValuesB[mMinIndex]
        m = mValuesB[mMinIndex]

    #updates the new values to the hinge object
    hinge.t1 = t1
    hinge.w = w
    hinge.D1 = D1
    print(f"t = {t1}", f"hole diameter = {D1}", f"width = {w}", m)

def CalcLugDimTwo(hinge):
    #explaination of what this thing does
    """
    The program runs through a range of values for thickness(t1), and finds hole diameter(D1) and width(w)
    that will meet the bearing and bending stress (hence the two different for loops. one checks for
    bearing, the other for bending. Then, the code finds a "mass(m)" that isnt really mass, but this value
    is proportional to mass. Thus, this value takes on a minimum when mass is a minimum.
    The equations show that D1 is inversely proportional to both bearing stress and bending stress, for
    constant w and t1. So, the optimiser finds both values of D1 for a given t, and puts the resulting tuples
    in a list(dList) along with their corresponding masses. An if statement checks which one is bigger, and
    appends either an A or a B. Then, a list(mList) is created with the corresponding m values. i.e. if the
    D1 required for bearing is bigger than the D1 required for bending, then an A is appended, and mList
    contains the m from the first for loop.
    The minimum(mMin) m from mList is used as the optimal value. The index(mMinIndex) of this minimum is
    used to find the optimal t, D1, and w values. The A or B appended earlier helps find out which list to
    take these values from.
    """

    t_step = 0.000001

    # Taking values from Loads and giving them nicer names
    P = Loads.P
    F1 = Loads.F1
    #Declare initial values
    D1 = hinge.D1
    w = hinge.w
    t1 = hinge.t1 #m
    Purple = True

    #Make lists to append bearing values of m, D1, and w to later
    mValuesA, DValuesA, wValuesA, tValues = [np.array([])] * 4

    while w < 1:
        D1 = hinge.D1
        while D1 < w:
            t1 = hinge.t1
            while t1 < 0.25:
                A_frac = 3 * (1 / D1) / (8 / (2 * w - 2 * D1 + math.sqrt(2) * D1) + 2 / (w - D1))  # p18 fig D1.15
                K_bending = 1.3 / 1.4 * A_frac  # Fig D1.15 page 18
                K_t = -0.05 * w / D1 + 3.05  # from appendix A tab D1.3, Curve 1 (W/D up to 3)
                check1 = 4 * K_bending * D1 * t1 * hinge.sigmaY
                check2 = 4 * K_t * hinge.sigmaY * t1 * (w - D1)
                if check1 > P[0]*1.5 and check2 > (P[1] + F1) * 1.5:
                    m = t1 * (w ** 2 - D1 ** 2)  # this value is not actually mass, but it is proportional to mass
                    mValuesA = np.append(mValuesA, m)
                    DValuesA = np.append(DValuesA, D1)
                    wValuesA = np.append(wValuesA, w)
                    tValues = np.append(tValues, t1)
                    # print(m, D1, w, t1)
                t1 += 0.001
            D1 += 0.001
            # print("d")
        w += 0.001
        # print("w")

    mMinIndex = np.argmin(mValuesA)

    t1 = tValues[mMinIndex]
    D1 = DValuesA[mMinIndex]
    w = wValuesA[mMinIndex]
    m = mValuesA[mMinIndex]


    #updates the new values to the hinge object
    hinge.t1 = t1
    hinge.w = w
    hinge.D1 = D1
    print(f"t = {t1}", f"hole diameter = {D1}", f"width = {w}", m)

def CalcLugDimThree(arr):
    t1, D1, w = arr
    P = Loads.P
    F1 = Loads.F1
    A_frac = 3 * (1 / D1) / (8 / (2 * w - 2 * D1 + math.sqrt(2) * D1) + 2 / (w - D1))  # p18 fig D1.15
    K_bending = 1.3 / 1.4 * A_frac  # Fig D1.15 page 18
    K_t = -0.05 * w / D1 + 3.05  # from appendix A tab D1.3, Curve 1 (W/D up to 3)
    check1 = 4 * K_bending * D1 * t1 * 4.14e7
    check2 = 4 * K_t * 4.14e7 * t1 * (w - D1)
    if check1 > P[0] * 1.5 and check2 > (P[1] + F1) * 1.5 and w > D1 + t1:
        return t1 * math.pi*(w ** 2 - D1 ** 2)/4
    else:
        return 31415
#-------------------------------------------------------------------------------------------------------------------------------------------


#4.4-----------------------------------------------------------------------------------------------------------------
def CalcBasePlateDim(hinge, e1Fac=1.5, e2Fac=1.5, holeSepFac=2.5, fastenerAmount=PD.fastenerAmount, fastenerColumns = PD.fastenerColumns):
    """
    calculates the dimensions of the baseplate with the width and the factors of seperation
    """

    hinge.D2 = hinge.w/(fastenerAmount/fastenerColumns + 2 * e1Fac + holeSepFac*(fastenerAmount/fastenerColumns-1))
    hinge.e1 = e1Fac*hinge.D2
    hinge.e2 = e2Fac*hinge.D2

    hinge.depth = 2*(2* hinge.e2 + hinge.t1 + fastenerColumns/2 * holeSepFac*hinge.D2) + hinge.h
#---------------------------------------------------------------------------------------------------------------------------


#TODO:calculate h
#4.5------------------------------------------------------------------------------------------------------------------------
def CalcFastenerPos(hinge, fastenerAmount = PD.fastenerAmount, columnAmount = PD.fastenerColumns):
    """
    calculates the position of the fasteners
    """

    #You can fil in a custom amount of fasteners

    #initialise the list of fasteners
    Fasteners = [] * fastenerAmount

    #define the list of z positions
    posZs = np.linspace((hinge.w/2 - hinge.e1), (-hinge.w/2 + hinge.e1), (fastenerAmount/columnAmount))
    #ceates the positive x positions
    posXs = np.linspace((hinge.depth/2 - hinge.e2), (hinge.t1 + hinge.h/2 + hinge.e2), columnAmount/2)
    #mirrors for the negative x positions
    posXs = np.append(posX, -posX[::-1])
    #creates a list with all the positions
    posTup = list(itertools.product(posXs, posZs))
    #define the fasteners and their positions


    for i in range(fastenerAmount):
        #calculate positions
        posX, posZ = posTup[i]

        #add to the list of fasteners
        Fasteners[i] = PD.Fastener(D_h = hinge.D2, xPos=posX, zPos=posZ)

    return Fasteners

def CalcCG(Fasteners):
    """
    calculate the centre of gravity of the fasteners
    
    The fasteners are defined first.
    The locacal coordinate system is the centre of the baseplate
    
    an even amount of fasteners is needed
    """

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

    return np.array([cgX, cgZ])
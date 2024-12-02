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

    t_step = 0.0001

    # Taking values from Loads and giving them nicer names
    P = Loads.P
    F1 = Loads.F1
    #Declare initial values
    D1 = hinge.D1
    w = hinge.w
    t1 = hinge.t1 #m

    #Make lists to append bearing values of m, D1, and w to later
    mValuesA, DValuesA, wValuesA = [np.array([])] * 3


    #Solves for Bearing stress for a range of t values
    while t1 < 0.25:
        D1 = (P[1] + F1) * 1.5 / 4 / hinge.sigmaY / t1 #Bearing stress (1.5 is MS)
        K_t = -0.05 * w/D1 + 3.05 #from appendix A tab D1.3, Curve 1 (W/D up to 3)
        w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t1 + D1 #Tension of net section
        m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
        mValuesA = np.append(mValuesA, m)
        DValuesA = np.append(DValuesA, D1)
        wValuesA = np.append(wValuesA, w)
        t1 += t_step


    #Reset initial values (crucial for the cursed K's)
    D1 = hinge.D1
    w = hinge.w
    t1 = hinge.t1 #m

    #Make lists to append bending values of m, D1, and w to later
    mValuesB, DValuesB, wValuesB = [np.array([])] * 3

    #Solves for bending stress for a range of t values
    while t1 < 0.25:
        if w > D1:
            A_frac = 6 / (D1 * (4/(0.5*w-math.sqrt(2)*0.25 * D1) + 2/(0.5*(w-D1)))) # p18 fig D1.15
        else:
            A_frac = 6 / (D1 * (4/(0.5*w-math.sqrt(2)*0.25 * D1) + 2/(0.5*(0.001))))
        K_bending = 1.3/1.4 * A_frac #Fig D1.15 page 18
        D1 = P[0] / 4 / (K_bending * hinge.sigmaY) / t1 #Bending
        K_t = -0.05 * w / D1 + 3.05 #from appendix A tab D1.3, Curve 1 (W/D up to 3)
        w = (P[1] + F1) * 1.5 / 4 / (K_t * hinge.sigmaY) / t1 + D1
        m = t1 * (w ** 2 - D1 ** 2) #this value is not actually mass, but it is proportional to mass
        mValuesB = np.append(mValuesB, m)
        DValuesB = np.append(DValuesB, D1)
        wValuesB = np.append(wValuesB, w)
        t1 += t_step

    dTest = (DValuesA >= DValuesB)
    mList = [a if c==True else b for a,b,c in zip(mValuesA, mValuesB, dTest)]
    mMinIndex = np.argmin(mList)

    if dTest[mMinIndex] == 1:
        t1 = mMinIndex * 0.001 + hinge.t1
        D1 = DValuesA[mMinIndex]
        w = wValuesA[mMinIndex]
    else:
        t1 = mMinIndex * 0.001 + hinge.t1
        D1 = DValuesB[mMinIndex]
        w = wValuesB[mMinIndex]

    #updates the new values to the hinge object
    hinge.t1 = t1
    hinge.w = w
    hinge.D1 = D1
    print(f"t = {t1}", f"hole diameter = {D1}", f"width = {w}")
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
    posZs = np.linspace((hinge.w/2 - hinge.e1), (-hinge.w/2 + hinge.e1), int(fastenerAmount/columnAmount))
    #ceates the positive x positions
    posXs = np.linspace((hinge.depth/2 - hinge.e2), (hinge.t1 + hinge.h/2 + hinge.e2), int(columnAmount/2))
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

#4.6----------------------------------------------------------
def CalcCGForces(Fasteners, CG):
    """
    Calculates the forces on each fastener
    """
    Fx, Fz = Loads.P[0], Loads.P[2]

    #Force on the cg
    Fcg = np.array([Fx, Fz])
    #Due to the definition of the axis system and the fact that the fastener pattern is symmetrical, there is no resulting moment at the CG of the fasteners, but this is calculated regardless
    FastenerLoads = np.empty(len(Fasteners))

    #moment on the cg
    Mycg = np.cross(CG, Fcg)

    #get the sum of the area and distance to calculate the moment load later on
    Sum = 0
    for Fast in Fasteners:
        posi = np.array([Fast.xPos, Fast.zPos])
        di = posi - CG
        Sum += Fast.D_h**2 * math.pi * 0.25 * (di**2)    


    #Calculate the loads on each fasteners
    i=0
    for Fast in Fasteners:
        load = np.zeros(3)
        pos = np.array([Fast.xPos, Fast.zPos])
        d = pos - CG

        load[0:2] = Fcg/(len(Fasteners))
        load[2] = (Mycg * Fast.D_h**2 * math.pi * 0.25 * d)/Sum

        #saves the load vector to the fastener object
        Fast.loadsInPlane = load
        FastenerLoads[i] = load
        i+=1
    #returns the loads as a matrix for debugging
    return FastenerLoads
    
#4.7-------------------------------------------------------------------------
def CheckBearing(hinge, Fasteners):
    result = (1,1)
    for Fast in Fasteners:
        P = np.linalg.norm(Fast.loadsInPlane)

        #bearing check for the baseplate thickness
        if P/(hinge.D2*hinge.t2) > hinge.SigmaB:
            result[0] = 0

        #bearing check for the spacecraft wall
        if P/(hinge.D2*hinge.t3) > hinge.SigmaB:
            result[1] = 0
        
    return result


#4.8-----------------------------------------------------------------------------------------------------------------

#this function calculates the load that could push or pull the fasteners through
def calcPullThroughLoad(yforce, zmoment, Fasteners):
    pullforce = yforce/len(Fasteners)
    cg =CalcCG(Fasteners)
    Sum = 0
    for fastener in Fasteners:
        posi = np.array([fastener.xPos, fastener.zPos])
        di = posi - cg
        Sum += fastener.D_h**2 * math.pi * 0.25 * (di**2)
    loads = []
    for fastener in Fasteners:
        momentload = zmoment * fastener.D_h**2 * math.pi * math.sqrt((fastener.xPos - cg[0])**2+(fastener.zPos - cg[1])**2) / Sum
        if fastener.zPos >= 0:
            loads.append(-momentload)
        else:
            loads.append(momentload)
    for load in loads:
        load = load + pullforce
    return loads
import math
import numpy as np
import scipy
import PartDefinition as PD
import Loads
import itertools


"""
Below are a list of functions which sequentially calculate the dimensions and the stresses of the hinge
These are made functions such that they can be ran again through iteration.

The result of each calculation will be Saved in the hinge object, such that one reference poijt exists.
"""

# 4.3------------------------------------------------------------------------------------------------------------------------------------
#set the minimum values o ore machinable values
def CalcLugDimOne(hinge, w_min=0.002):
    resulto = scipy.optimize.minimize(CalcLugDimTwo, [hinge.t1, hinge.D1, hinge.w], bounds=scipy.optimize.Bounds([0.0015, 0.008, w_min], [0.25, 0.5, 0.5]))
    print("\n",resulto,"\n")
    hinge.t1 = resulto.x[0]
    hinge.D1 = resulto.x[1]
    hinge.w = resulto.x[2]
    print(f"t1 is {hinge.t1}", f"w is {hinge.w}", f"D1 is {hinge.D1}")


def CalcLugDimTwo(arr):
    t1, D1, w = arr
    P = Loads.P
    F1 = Loads.F1
    sigma = 4.14e7 # Pa

    A_frac = 6 / (D1 * (4 / (0.5 * w - math.sqrt(2) * 0.25 * D1) + 2 / (0.5 * (w - D1))))  # p18 fig D1.15
    K_bending = 1.3 / 1.4 * A_frac  # Fig D1.15 page 18
    K_t = -0.05 * w / D1 + 3.05  # from appendix A tab D1.3, Curve 1 (W/D up to 3)
    check1 = 4 * K_bending * D1 * t1 * sigma
    check2 = 4 * K_t * sigma * t1 * (w - D1)

    if check1 > P[0] * 1.5 and check2 > (P[1] + F1) * 1.5 and w > D1 + t1:
        #prints the margin of safety
        print("MS LUG: ",check1/(P[0] * 1.5) - 1, check2/((P[1] + F1) * 1.5) - 1)
        return t1 * math.pi * (w ** 2 - D1 ** 2) / 4
    else:
        return 31415


#------------------------------------------------------------------------------------------------------------------------------------------


#4.4-----------------------------------------------------------------------------------------------------------------
def CalcBasePlateDim(hinge, e1Fac=1.5, e2Fac=1.5, holeSepFac=2, fastenerAmount=PD.fastenerAmount, fastenerColumns = PD.fastenerColumns, minD2 = 0):
    """
    calculates the dimensions of the baseplate with the width and the factors of seperation
    """
    
    #calc the minimum width
    w_min = minD2 * (2 * e1Fac + holeSepFac*(fastenerAmount/fastenerColumns-1))
    #updates the width and quits to rerun the first function
    if w_min > hinge.w:
        print("\nAdjusting Width size\n")
        hinge.w = w_min
        return (False, w_min)
    
    hinge.D2 = hinge.w/(2 * e1Fac + holeSepFac*(fastenerAmount/fastenerColumns-1))

    hinge.e1 = e1Fac*hinge.D2
    hinge.e2 = e2Fac*hinge.D2

    hinge.depth = 2*(2* hinge.e2 + hinge.t1 + (fastenerColumns/2 - 1) * holeSepFac*hinge.D2) + hinge.h
    return (True, 0.002)
#---------------------------------------------------------------------------------------------------------------------------

#4.5------------------------------------------------------------------------------------------------------------------------
def CalcFastenerPos(hinge, fastenerAmount = PD.fastenerAmount, columnAmount = PD.fastenerColumns):
    """
    calculates the position of the fasteners
    """

    #You can fil in a custom amount of fasteners

    #initialise the list of fasteners
    Fasteners = [None] * fastenerAmount

    #define the list of z positions
    posZs = np.linspace((hinge.w/2 - hinge.e1), (-hinge.w/2 + hinge.e1), int(fastenerAmount/columnAmount))
    #ceates the positive x positions
    posXs = np.linspace((hinge.depth/2 - hinge.e2), (hinge.t1 + hinge.h/2 + hinge.e2), int(columnAmount/2))
    #mirrors for the negative x positions
    posXs = np.append(posXs, -posXs[::-1])
    #creates a list with all the positions
    posTup = list(itertools.product(posXs, posZs))
    #define the fasteners and their positions

    for i in range(fastenerAmount):
        #calculate positions
        posX, posZ = posTup[i]

        #add to the list of fasteners TODO: define geometry clearly
        Fasteners[i] = PD.Fastener(d_uh_brg = 1.7*hinge.D2, D_h = hinge.D2, L = hinge.t2 + hinge.t2, L_h=0.64*hinge.D2, L_n=0.019, xPos=posX, zPos=posZ)

    return Fasteners

def CalcCG(Fasteners):
    """
    calculate the centre of gravity of the fasteners
    
    The fasteners are defined first.
    The local coordinate system is the centre of the baseplate
    
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
    Fx, Fz = Loads.P[0]/2 + Loads.T[1]*Loads.H, Loads.P[2]/2

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
        Sum += Fast.D_h**2 * math.pi * 0.25 * (np.linalg.norm(di)**2)    
    


    #Calculate the loads on each fasteners
    for Fast in Fasteners:
        load = np.zeros(3)
        pos = np.array([Fast.xPos, Fast.zPos])
        d = pos - CG

        load[0:2] = Fcg/(len(Fasteners))
        load[2] = (Mycg * Fast.D_h**2 * math.pi * 0.25 * np.linalg.norm(d))/Sum

        #saves the load vector to the fastener object
        Fast.loadsInPlane = load
    #returns the loads as a matrix for debugging
    return FastenerLoads
    
#4.7-------------------------------------------------------------------------
def CheckBearing(hinge, Fasteners):
    result = [1,1]
    MS = []
    for Fast in Fasteners:
        P = np.linalg.norm(Fast.loadsInPlane) * 1.5

        #bearing check for the baseplate thickness
        if P/(hinge.D2*hinge.t2) > hinge.SigmaB:
            result[0] = 0

        #bearing check for the spacecraft wall
        if P/(hinge.D2*hinge.t3) > hinge.SigmaB:
            result[1] = 0
    MS = [hinge.SigmaB/(P/(hinge.D2*hinge.t2)) - 1, hinge.SigmaB/(P/(hinge.D2*hinge.t3)) - 1]
    return (result, MS)


#4.8-----------------------------------------------------------------------------------------------------------------

"""this function calculates the load that could push or pull the fasteners through, the pull through force is composed of the force in y direction (force acting along the center
of the bolt) and the moment around the z axis. The load on each bolt due to the force in y direction is equal for every bolt and is simply the force in y direction divided
by the number of bolts. The moment depends on the position of the fastener so it needs to be determined for each fastener separately. For this the below function will iterate 
through all fasteners, and return a list of the loads on all the fasteners."""

def calcPullThroughLoad(Fasteners):
#changed to add the moment couple
    pullforce = (Loads.P[1])/(len(Fasteners)*2) + (Loads.F1)/(len(Fasteners))
    #calculate force on each bolt due the force in y direction
    cg =CalcCG(Fasteners)   #calculate the cg of all the fasteners
    Sum = 0     #initialise variable for the sum
    for fastener in Fasteners:  #iterate through all fasteners, find their distance to cg and sum up their area times the square of their distance
        posi = np.array([fastener.xPos, fastener.zPos])
        di = posi - cg
        Sum += fastener.D_h**2 * math.pi * 0.25 * (di[0]**2)

    for fastener in Fasteners:  #iterate through fasteners again and determine the force due to the moment and its sign and edit the out of plane force for each fastener
        posi = np.array([fastener.xPos, fastener.zPos])
        di = posi - cg
        momentload = Loads.T[2] * fastener.D_h**2 * math.pi* 0.25 * di[0] / (Sum*2)
        fastener.loadsOutPlane = (pullforce + momentload)
    return

#4.9---------------------------------------------------------------------------------------------------------------------

"""this function performs a pull through test. the variable load is the force that is vertically applied to the load. it causes the fastener head to apply a compressive
force on the hinge wall/spacecraft wall and in the skin material there will be a shear force. In order to be safe and conservative we will consider the failure of one of the walls
to be failure, therefore the vonmises stress will be calculated for both walls and the greater stress will be checked for failure."""

def CheckPullThrough(Fasteners, hinge):
    check = []
    MS = []
    for fastener in Fasteners:   #iterate through all fastener objects
        areabolthead = (fastener.D_fo/2)**2 * math.pi - (fastener.D_fi/2)**2 * math.pi    #calculate area on which compressive stress acts
        sigmay = fastener.loadsOutPlane/areabolthead      #calculate compressive stress
        areat2 = math.pi * fastener.D_fo * hinge.t2   #calculate areas over which the shear stress will act
        areat3 = math.pi * fastener.D_fo * hinge.t3
        tau2 = fastener.loadsOutPlane/areat2      #calculate shear stresses
        tau3 = fastener.loadsOutPlane/areat3
        if areat2 <= areat3:   #calculate von mises stress for greater shear stress
            vonmises = math.sqrt(sigmay**2 + 3 * tau2**2)
        else:
            vonmises = math.sqrt(sigmay**2 + 3 * tau3**2)
        
        if vonmises < hinge.sigmaY:  #if vonmises stress is below tensile yield stress, test is passed and add true, if vonmises stress is higher add false to list
            check.append(1)
            MS.append(hinge.sigmaY/vonmises - 1)

        else:
            check.append(0)
            MS.append(hinge.sigmaY/vonmises - 1)
        
    return check, MS

#4.10--------------------------------------------------------------------------
def CalcComplianceA(hinge, t, Fast):
    """
    calculates the compliance of the sheet that it is conneted to the fastener
    """
    delA = (4* t *hinge.E)/(math.pi*(Fast.D_fo**2 - Fast.D_fi**2)) #assume the youngs modulus is similar to the backplate one
    
    return delA 
#compliance B is calculated in the fastener object

import Main
import PartDefinition as PD
import numpy as np
from pprint import pprint

#specify the minimum bolt diameter


# initial definition of the hinge object
#change t2, t3 to more resonable begin values
hinge = PD.Hinge(t1=0.001, t2=0.0005, t3=0.0005, D1=0.01, w=0.02, sigmaY=4.14e7, SigmaB=297e6)

# runs the functions for the first time
min_check = (False, 0.002)
while min_check[0] == False:
    Main.CalcLugDimOne(hinge, min_check[1])
    min_check = Main.CalcBasePlateDim(hinge, minD2=PD.minBoltD)
Fasteners = Main.CalcFastenerPos(hinge) #fatsener dimensions are also defined here
FastCG = Main.CalcCG(Fasteners)
Main.CalcCGForces(Fasteners, FastCG)

print("\n----------------------------------------------------------\n")
#if bearingcheck returns false, we should increase the thickness
checkResult, MS = Main.CheckBearing(hinge,Fasteners)
print(checkResult, "\nBearing check MS: ", MS, "\n")
while 0 in checkResult:
    updateVal = np.abs(np.array(checkResult) - 1) * 0.0005
    hinge.t2 += updateVal[0]
    hinge.t3 += updateVal[1]
    checkResult, MS = Main.CheckBearing(hinge, Fasteners)
    print(checkResult, "\nBearing check MS: ", MS, "\n")

print("\n----------------------------------------------------------\n")
#Pullthrough check
Main.calcPullThroughLoad(Fasteners)
checkResult, MS = Main.CheckPullThrough(Fasteners, hinge)
print(checkResult, "\npullthrough check MS: ", MS, "\n")
while 0 in checkResult:
    hinge.t2 += 0.0015
    hinge.t3 += 0.0015
    Main.calcPullThroughLoad(Fasteners)
    checkResult, MS = Main.CheckPullThrough(Fasteners, hinge)
    print(checkResult, "\npullthrough check MS: ", MS, "\n")


"""
calculate some lengths of the fasteners
"""

DelA_bp = Main.CalcComplianceA(hinge, hinge.t2, Fasteners[0]) #calculate the compliance of the backplate and the spacecraft wall
DelA_w = Main.CalcComplianceA(hinge, hinge.t3, Fasteners[0]) #calculate the compliance of the backplate and the spacecraft wall
Phi = []
for Fast in Fasteners:
    Fast.CalcComplianceB() #calculates the compliance of each bolt
    Phi.append([DelA_bp/(DelA_bp + Fast.comp), DelA_w/(DelA_w + Fast.comp)])

print("\n_____________________________FINAL DIMENSIONS____________________________\n")
pprint(vars(hinge))
print("---------------")
pprint(vars(Fasteners[0]))
print(Phi)



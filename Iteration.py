import Main
import PartDefinition as Pd
import numpy as np
from pprint import pprint

#initial definition of the hinge obect
hinge = PD.Hinge(t1=0.001, t2=0.01, t3=0.01, D1=0.01, w=0.02, sigmaY=2.5e8, SigmaB=297e6)

# runs the functions for the first time
Main.CalcLugDimOne(hinge)
print("finishedCalclugdim\n")
Main.CalcBasePlateDim(hinge)
Fasteners = Main.CalcFastenerPos(hinge)
FastCG = Main.CalcCG(Fasteners)
Main.CalcCGForces(Fasteners, FastCG)
print("Bearingcheck")
checkResult = Main.CheckBearing(hinge,Fasteners)
#if bearingcheck returns false, we should increase the thickness
while 0 in checkResult:
    print(checkResult)
    
    updateVal = np.abs(np.array(checkResult) - 1) * 0.001
    hinge.t2 += updateVal[0]
    hinge.t3 += updateVal[1]
    checkResult = Main.CheckBearing(hinge, Fasteners)

"""
calculate some lengths of the fasteners
"""
for Fast in Fasteners:
    Fast.Lj = hinge.t2 + hinge.t3
    
pprint(vars(hinge))
print("---------------")
pprint(vars(Fasteners[0]))


DelA_bp = Main.CalcComplianceA(hinge, hinge.t2, Fasteners[0]) #calculate the compliance of the backplate and the spacecraft wall
DelA_w = Main.CalcComplianceA(hinge, hinge.t3, Fasteners[0]) #calculate the compliance of the backplate and the spacecraft wall


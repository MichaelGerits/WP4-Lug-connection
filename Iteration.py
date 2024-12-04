import Main
import PartDefinition as Pd
import numpy as np

# initial definition of the hinge object
hinge = Pd.Hinge(t1=0.001, t2=0.01, t3=0.01, D1=0.01, w=0.02, sigmaY=250000000)

# runs the functions for the first time
Main.CalcLugDimOne(hinge)
Main.CalcBasePlateDim(hinge)
Fasteners = Main.CalcFastenerPos(hinge)
FastCG = Main.CalcCG(Fasteners)
Main.CalcCGForces(Fasteners, FastCG)

checkResult = Main.CheckBearing(hinge, Fasteners)
# if bearing check returns false, we should increase the thickness
while 0 in checkResult:
    
    updateVal = np.abs(np.array(checkResult) - 1) * 0.001
    hinge.t2 += updateVal[0]
    hinge.t3 += updateVal[1]
    checkResult = Main.CheckBearing(hinge, Fasteners)

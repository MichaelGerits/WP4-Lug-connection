import Main
import PartDefinition as PD
import Loads
import scipy
import numpy as np

#initial definition of the hinge obect
hinge = PD.Hinge(t1=0.001, t2=0.01, t3=0.01, D1=0.01, w=0.02, sigmaY=250000000)

#runs the functions for the first time
#TODO:Save optimised functions
print(scipy.optimize.minimize(Main.CalcLugDimThree, [0.01, 0.01, 0.02], bounds=scipy.optimize.Bounds([0.001, 0.001, 0.002], [0.25, 0.5, 0.5])))
Main.CalcBasePlateDim(hinge)
Fasteners = Main.CalcFastenerPos(hinge)
FastCG = Main.CalcCG(Fasteners)
Main.CalcCGForces(Fasteners, FastCG)

checkResult = Main.CheckBearing(hinge,Fasteners)
#if bearingcheck returns false, we should increase the thickness
while 0 in checkResult:
    
    hinge.t2, hinge.t3 = np.abs(np.array(checkResult) - 1) * 0.001
    hinge.t2 += 0.001
    checkResult = Main.CheckBearing(hinge, Fasteners)

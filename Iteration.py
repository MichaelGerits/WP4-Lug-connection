import Main
import PartDefinition as PD
import Loads

#initial definition of the hinge obect
hinge = PD.Hinge(t1=0.001, t2 = 0.01, D1=0.01, w=0.02 ,sigmaY=250000000)

#runs the functions for the first time
Main.CalcLugDim(hinge)
Main.CalcBasePlateDim(hinge)
Fasteners = Main.CalcFastenerPos(hinge)
FastCG = Main.CalcCG(Fasteners)
Main.CalcCGForces(Fasteners, FastCG)

BearingCheck = Main.CheckBearing(hinge,Fasteners)
#if bearingcheck returns false, we should increase the thickness
while BearingCheck == False:
    hinge.t2 += 0.001
    BearingCheck(hinge, Fasteners)
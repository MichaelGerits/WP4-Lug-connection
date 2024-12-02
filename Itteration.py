import Main
import PartDefinition as PD
import scipy
import Loads

#initial definition of the hinge obect
hinge = PD.Hinge(t1=0.0011, D1=0.0011, w=0.0021, sigmaY=414000000)

#runs the functions for the first time
print(scipy.optimize.minimize(Main.CalcLugDimThree, [0.01, 0.01, 0.02], bounds=scipy.optimize.Bounds([0.001, 0.001, 0.002], [0.25, 0.5, 0.5])))
# Main.CalcLugDimTwo(hinge)
# Main.CalcBasePlateDim(hinge)
# Fasteners = Main.CalcFastenerPos(hinge)
# FastCG = Main.CalcCG(Fasteners)
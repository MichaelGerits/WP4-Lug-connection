import math as m
import numpy as np
import Loads

"""
Below the definition goes as follows: A hinge is the entirety and the Lugs are the flanges of said hinge
"""

fastenerAmount = 4 #the amount of fasteners per Hinge
fastenerColumns = 2

class Hinge:
    """
    define the geometry and materials properties of the Hinge
    as well as any calculation
    """
    def __init__(self, h=Loads.H/3, w=0, D1=0, D2=0, t1=0, t2=0, t3=0, e1=0, e2=0, E=70e9, G=0, rho=0, sigmaY=0, SigmaB=0, depth=0) -> None:
        self.h = h
        self.D1 = D1
        self.D2 = D2
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.e1 = e1 #limits depending on the material choice
        self.e2 = e2
        self.w = w
        self.depth = depth

        self.E = E
        self.G = G
        self.rho = rho
        self.sigmaY = sigmaY
        self.SigmaB = SigmaB
        pass

class Fastener:
    """
    define geometry material properties and forces on a bolt
    geometry taken from WP4 p20-21
    """
    def __init__(self, d_uh_brg=0, L_h=0, D_h=0, d_sha=0, L_n=0, d=0, sw=0, L1=0, L2=0, L3=0, E_b=0, E_n=70e9, G=0, sigmaY=0, rho=0, xPos = 0, zPos = 0) -> None:
        self.d_uh_brg = d_uh_brg
        self.L_h = L_h
        self.D_h = D_h
        self.d_sha = d_sha
        self.d = d
        self.D_fo = d_uh_brg
        self.D_fi = D_h
        self.sw = sw

        self.L_n = L_n
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.Lj = L1+L2+L3 #sum of the previous three lengths as per page 21
        self.L_h_sub = 0.5*d #TODO: decide on nut geometry to find factor
        self.L_eng_sub = 0.4*d #TODO: decide on nut geometry to find factor Table 7.1 page 22
        self.L_n_sub = 0.4*d #TODO: decide on nut geometry to find factor

        self.E_b = E_b
        self.E_n = E_b if E_n==0 else E_n
        self.G = G
        self.sigmaY = sigmaY
        self.rho = rho

        self.xPos = xPos
        self.zPos = zPos

        self.loadsInPlane = np.empty(3) #due to Fx, due to Fz, due to My
        self.loadsOutPlane = 0 #due to Fx, due to Fz, due to My
        pass

    def CalcComplianceB(self):
        self.comp = 0

        #TODO: sum up the elongations and divide by the youngs modulus
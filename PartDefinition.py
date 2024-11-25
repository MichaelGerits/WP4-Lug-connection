import math as m
import numpy as np

"""
Below the definition goes as follows: A hinge is the entirety and the Lugs are the flanges of said hinge
"""

fastenerAmount = 4 #the amount of fasteners per Hinge

class Hinge:
    """
    define the geometry and materials properties of the Hinge
    as well as any calculation
    """
    def __init__(self, h=0, w=0, D1=0, D2=0, t1=0, t2=0, e1=0, e2=0, E=0, G=0, rho=0, sigmaY=0) -> None:
        self.h = h
        self.w = w
        self.D1 = D1
        self.D2 = D2
        self.t1 = t1
        self.t2 = t2
        self.e1 = e1 #limits depending on the material choice
        self.e2 = e2

        self.E = E
        self.G = G
        self.K_t = -0.05 * w/D1 + 3.05 #from appendix A tab D1.3, Curve 1 (W/D up to 3)
        self.rho = rho
        self.sigmaY = sigmaY
        self.K_bry =1.85 #P17 eq9 says so

        """
        Variables for bending of the lug 
        """
        self.A_frac = 6 / (D1 * (4/(0.5*w-m.sqrt(2)*0.25 * D1) + 2/(0.5*(w-D1)))) # p18 fig D1.15
        K_bending = 
        pass

class Fastener:
    """
    define geometry material properties and forces on a bolt
    geometry taken from WP4 p20-21
    """
    def __init__(self, d_uh_brg=0, L_h=0, D_h=0, d_sha=0, L_n=0, d=0, sw=0, L1=0, L2=0, L3=0, E_b=0, E_n=0, G=0, sigmaY=0, rho=0) -> None:
        self.d_uh_brg = d_uh_brg
        self.L_h = L_h
        self.D_h = D_h
        self.d_sha = d_sha
        self.d = d
        self.sw = sw

        self.L_n = L_n
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.Lj = L1+L2+L3 #sum of the previous three lengths as per page 21
        self.L_h_sub = x*d #TODO: decide on nut geometry to find factor
        self.L_eng_sub = x*d #TODO: decide on nut geometry to find factor Table 7.1 page 22
        self.L_n_sub = x*d #TODO: decide on nut geometry to find factor

        self.E_b = E_b
        self.E_n = E_b if E_n==0 else E_n
        self.G = G
        self.sigmaY = sigmaY
        self.rho = rho
        pass

    def CalcCompliance(self):
        comp = 0
        #TODO: sum up the elongations and divide by the youngs modulus
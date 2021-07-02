import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *

class Photon:

    def __init__(self, 
                 amin,      # scale factor at emssion 
                 amax,          # scale factor at receiving 
                 ):
        Result = self.Dc(amin, amax)
        self.a = Result[0]
        self.t = Result[1]
        return

    def Dc(self,
           amin,        # The scale factor of when photon was launched
           amax,            # scale factor at receiving 
           N = int(1e3)     # Sample points
           ):
        
        a = np.linspace(amin, amax, N)
       
        # Hubble parameter evolve with scale factor
        def Hubble_Parameter(a):
            # Assuming k = 0, standard LCDM cosmology
            return H0*np.sqrt(Omega_R*a**-4 + Omega_M*a**-3 + Omega_L*a**0)
        
        # The comoving distance photon travels
        def Cov_Dis(a, d):
            H = Hubble_Parameter(a) # Hubble Parameter
            return c/(a**2*H)*ray_direction
        
        # Cosmic time
        def cosmic_time(a, t):
            H = Hubble_Parameter(a) # Hubble Parameter
            return 1/(a*H)
        
        
        # Relation between a and t
        source_t = solve_ivp(   fun=cosmic_time, 
                                t_span=(amin, amax),
                                t_eval=a,
                                y0=[0])
        
        t = source_t.y[0]        # Cosmic time in s
        return [a, t]

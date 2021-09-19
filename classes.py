import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *

class Photon:

    def __init__(self, 
                 a_launch,      # scale factor at emssion 
                 amax,          # scale factor at receiving 
                 d_source,      # The distance of the photon source in comoving distance
                 ray_direction, # +1 for travel outward, -1 for inward
                 N = int(1e3)):
        Result = self.Dc(a_launch, amax, d_source, ray_direction, N=N)
        self.a_photon = Result[0]
        self.t_photon = Result[1]
        self.d_photon = Result[2]
        self.a_source = Result[3]
        self.t_source = Result[4]
        self.d_source = Result[5]
        return

    def Dc(self,
           a_launch,        # The scale factor of when photon was launched
           amax,            # scale factor at receiving 
           d_source,        # The comoving distance of the photon source  
           ray_direction,   # +1 for travel outward, -1 for inward
           amin=1e-30,      # scale factor at emssion 
           N = int(1e3)     # Sample points
           ):
        
        a_source = np.linspace(amin, amax, N)
        a_photon = np.linspace(a_launch, amax, N)
       
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
        
        # Path of photon
        photon = solve_ivp(     fun=Cov_Dis, 
                                t_span=(a_launch, amax),
                                t_eval=a_photon,
                                y0=[d_source])
        
        # Relation between a and t
        # To avoid a=0 at the beginning of the universe, 
        # we start from current universe to work out t(a)
        tmin = solve_ivp(       fun=cosmic_time, 
                                t_span=(1, amin),
                                y0=[138e8*yr2s]).y[0][-1]

        # Note: We are using 13.8 Gyr as the age of the universe
        # This is cosmological parameter dependent!
        t_launch = solve_ivp(   fun=cosmic_time, 
                                t_span=(1, a_launch),
                                y0=[138e8*yr2s]).y[0][-1]
        
        source_t = solve_ivp(   fun=cosmic_time, 
                                t_span=(amin, amax),
                                t_eval=a_source,
                                y0=[tmin])
        
        photon_t = solve_ivp(   fun=cosmic_time, 
                                t_span=(a_launch, amax),
                                t_eval=a_photon,
                                y0=[t_launch])

        t_photon = photon_t.y[0]        # scale factor
        d_photon = a_photon*photon.y[0] # Proper distance of photon in MKS
        t_source = source_t.y[0]        # Cosmic time in s
        d_source = a_source*d_source    # Proper distance of source in MKS
        return [a_photon, t_photon, d_photon, a_source, t_source, d_source]

    


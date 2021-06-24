#====================#

# This code is made by Lin Yen-Hsing (NTHU, 2017-2021) to demonstate
# the difference between three commonly mis-understood concepts in
# cosmoloy, namely the "Particle horizon", the "Event horizon" and 
# the "Hubble radius".

# This code mainly based on 
# https://physics.stackexchange.com/questions/70339/a-cosmological-horizon-at-the-hubble-radius
# which provides a clear and short description on this topic.
# Of course, these concepts should be found in any cosmology textbook.

# If you have any question, feel free to contact me with
# julius52700@gapp.nthu.edu.tw

#====================#

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Fundamental Constants
c = 299792458           # m/s
Mpc2m = 3.08567758e22   # m
Gly2m = 9.4605e24       # m
yr2s = 31556926         # s
rho_c = 7.64e-10        # J/m^3
H0 = 67.3e3/Mpc2m       # Hubble Constant at a=1
Omega_M = 0.315         # Mass density
Omega_L = 0.685         # Dark energy density
Omega_R = 4.17e-14/rho_c# Radiation density

# Create figure
#plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(6,6))
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
#ax.set_facecolor('#303030')

class Photon:

    def __init__(self, 
                 a_launch,      # scale factor at emssion 
                 amax,          # scale factor at receiving 
                 d_source,      # The distance of the photon source in comoving distance
                 ray_direction  # +1 for travel outward, -1 for inward
                 ):
        Result = self.Dc(a_launch, amax, d_source, ray_direction)
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

# This is particle horizon demo
ray1 = Photon(1e-50, 10, 47*Gly2m, -1) 
ax.plot(ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='b') # This is photon
ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='b') # This is source galaxy
ray1 = Photon(1e-50, 10, -47*Gly2m, 1) 
ax.plot(ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='b') # This is photon
ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='b') # This is source galaxy

# This is event horizon demo
ray2 = Photon(1, 10, 16.7*Gly2m, 1) 
#ray3 = Photon(1, 10, 0, 1) # photon launch by observer
#ax.plot(ray3.d_photon/Gly2m, ray3.t_photon/yr2s/1e9, linestyle='solid',  color='r') # This is photon
ax.plot(ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='r') # This is source galaxy
ray2 = Photon(1, 10, -16.7*Gly2m, -1) 
#ray3 = Photon(1, 10, 0, -1) # photon launch by observer
#ax.plot(ray3.d_photon/Gly2m, ray3.t_photon/yr2s/1e9, linestyle='solid',  color='r') # This is photon
ax.plot(ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='r') # This is source galaxy

# This is for hubble radius
ray4 = Photon(1, 10, 14.5*Gly2m, 1) 
ax.plot(ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='goldenrod') # This is source galaxy
ray4 = Photon(1, 10, -14.5*Gly2m, -1) 
ax.plot(ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='goldenrod') # This is source galaxy

ax.axhline(13.8, 0, 1, linestyle='dashed', color='gray', label='Now') # Current age of the universe
ax.axvline(0,    0, 1, linestyle='dashed', color='gray')              # World line for observer
'''
# RHS
ax.axvline(16.7, 0, 1, linestyle='dashed', color='r',    label='Event Horizon')
ax.axvline(14.5, 0, 1, linestyle='dashed', color='green', label='Hubble Radius')
ax.axvline(47,   0, 1, linestyle='dashed', color='b',    label='Particle Horizon')

# LHS
ax.axvline(-16.7, 0, 1, linestyle='dashed', color='r',   ) 
ax.axvline(-14.5, 0, 1, linestyle='dashed', color='green')
ax.axvline(-47,   0, 1, linestyle='dashed', color='b',   ) 
'''
ax.set_xlim([-125, 125])
ax.set_ylim([0, 50 ])
#ax.set_xlabel('Proper Distance (Gly)')
#ax.set_ylabel('Cosmic Time (Gyr)')
#plt.legend()
fig.savefig('World_Line_Cosmology.png', dpi=300)
plt.show()

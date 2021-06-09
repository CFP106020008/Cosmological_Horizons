#====================#

# This code is made by Lin Yen-Hsing (NTHU, 2017-2021) to demonstate
# the difference between three commonly mis-understood concepts in
# cosmoloy, namely the "Particle horizon", the "Event horizon" and 
# the "Hubble radius".

# This code is mainly based on 
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
c = 299792458 # m/s
Mpc2m = 3.08567758e22 # m
Gly2m = 9.4605e24 # m
yr2s = 31556926 # s
rho_c = 7.64e-10 # J/m^3
H0 = 67.3e3/Mpc2m # Hubble Constant at a=1
Omega_M = 0.315 # Mass density
Omega_L = 0.685 # Dark energy density
Omega_R = 4.17e-14/rho_c

fig, ax = plt.subplots()

# Comoving distance

class Photon:

    def __init__(self, 
                 amin, # scale factor at emssion 
                 amax, # scale factor at receiving 
                 d_source, # The distance of the photon source in comoving distance
                 ray_direction # +1 for travel outward, -1 for inward
                 ):
        Result = self.Dc(amin, amax, d_source, ray_direction)
        self.a = Result[0]
        self.t = Result[2]
        self.photon = Result[1]
        self.source = Result[3]
        return

    def Dc(self,
           amin, # scale factor at emssion 
           amax, # scale factor at receiving 
           d_source, # The distance of the photon source in comoving distance
           ray_direction, # +1 for travel outward, -1 for inward
           N = int(1e3) # Sample points
           ):
        
        a = np.linspace(amin, amax, N)
        
        # The comoving distance photon travels
        def Hubble_Parameter(a):
            return H0*np.sqrt(Omega_R*a**-4 + Omega_M*a**-3 + Omega_L*a**0)
        
        def Cov_Dis(a, d): 
            H = Hubble_Parameter(a) # Hubble Parameter
            return c/(a**2*H)*ray_direction
        D = solve_ivp(fun=Cov_Dis, 
                      t_span=(amin, amax),
                      t_eval=a,
                      y0=[d_source])
        
        # Cosmic time
        def cosmic_time(a, t):
            H = Hubble_Parameter(a) # Hubble Parameter
            return 1/(a*H)
        
        # To avoid a=0 at the beginning of the universe, 
        # we start from current universe to work out t(a)
        tmin = solve_ivp(fun=cosmic_time, 
                         t_span=(1, amin),
                         y0=[138e8*yr2s])
        
        T = solve_ivp(fun=cosmic_time, 
                      t_span=(amin, amax),
                      t_eval=a,
                      y0=[tmin.y[0][-1]])
        
        a = D.t                 # scale factor
        d = a*D.y[0]            # Proper distance of photon in MKS
        t = T.y[0]              # Cosmic time in s
        source = a*d_source     # Proper distance of source in MKS
        return [a, d, t, source]

# This is particle horizon demo
ray1 = Photon(1e-50, 1, 47*Gly2m, -1) 
ax.plot(ray1.photon/Gly2m, ray1.t/yr2s/1e9, linestyle='solid',  color='b') # This is photon
ax.plot(ray1.source/Gly2m, ray1.t/yr2s/1e9, linestyle='dotted', color='b') # This is source galaxy

# This is event horizon demo
ray2 = Photon(1, 10, 16.7*Gly2m, 1) 
ray3 = Photon(1, 10, 0, 1) # photon launch by observer
ax.plot(ray3.photon/Gly2m, ray3.t/yr2s/1e9, linestyle='solid',  color='r') # This is photon
ax.plot(ray2.source/Gly2m, ray2.t/yr2s/1e9, linestyle='dotted', color='r') # This is source galaxy

# This is for hubble radius
ray4 = Photon(1, 10, 14.5*Gly2m, 1) 
ax.plot(ray4.source/Gly2m, ray4.t/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy

ax.axhline(13.8, 0, 1, linestyle='dashed', color='gray', label='Now') # Current age of the universe
ax.axvline(0, 0, 1,    linestyle='dashed', color='gray') # World line for observer
ax.axvline(16.7, 0, 1, linestyle='dashed', color='r',    label='Event Horizon')
ax.axvline(14.5, 0, 1, linestyle='dashed', color='cyan', label='Hubble Radius')
ax.axvline(47,   0, 1, linestyle='dashed', color='b',    label='Particle Horizon')
ax.set_xlim([0, 125])
ax.set_ylim([0, 50])
ax.set_xlabel('Proper Distance (Gly)')
ax.set_ylabel('Cosmic Time (Gyr)')
plt.legend()
fig.savefig('World_Line_Cosmology.png', dpi=300)
plt.show()

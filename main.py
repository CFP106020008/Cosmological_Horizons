import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

# Create figure
#plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(6,6))
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
plt.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
#ax.set_facecolor('#303030')

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

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

# Create figure
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(5,5))
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
plt.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
#ax.set_facecolor('#303030')
'''
# This is particle horizon demo
ray1 = Photon(1e-50, 10, 47*Gly2m, -1) 
ax.plot(ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='b') # This is photon
ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='b') # This is source galaxy
ray1 = Photon(1e-50, 10, -47*Gly2m, 1) 
ax.plot(ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='b') # This is photon
ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='b') # This is source galaxy
'''
# This is event horizon demo
ray2 = Photon(1, 50, 1*Gly2m, 1) 
ray3 = Photon(1, 50, 5*Gly2m, 1) 
ray4 = Photon(1, 50, 10*Gly2m, 1) 
#ray3 = Photon(1, 10, 0, 1) # photon launch by observer
#ray4 = Photon(0.35, 1, 16.7*Gly2m, -1) # photon that reach observer now
#ax.plot(ray3.d_photon/Gly2m, ray3.t_photon/yr2s/1e9, linestyle='solid',  color='y', label='Photon launched now') # This is photon
#ax.plot(ray4.d_photon/Gly2m, ray4.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
ax.plot(ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 1 Gly') # This is source galaxy
ax.plot(ray3.d_source/Gly2m, ray3.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 5 Gly') # This is source galaxy
ax.plot(ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 10 Gly') # This is source galaxy
ax.plot(-ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
ax.plot(-ray3.d_source/Gly2m, ray3.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
ax.plot(-ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
#ray2 = Photon(1, 10, -16.7*Gly2m, -1) 
#ray3 = Photon(1, 10, 0, -1) # photon launch by observer
#ax.plot(ray3.d_photon/Gly2m, ray3.t_photon/yr2s/1e9, linestyle='solid',  color='r') # This is photon
#ax.plot(ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='r') # This is source galaxy
'''
# This is for hubble radius
ray4 = Photon(1, 10, 14.5*Gly2m, 1) 
ax.plot(ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='goldenrod') # This is source galaxy
ray4 = Photon(1, 10, -14.5*Gly2m, -1) 
ax.plot(ray4.d_source/Gly2m, ray4.t_source/yr2s/1e9, linestyle='dotted', color='goldenrod') # This is source galaxy
'''
ax.axhline(13.8, 0, 1, linestyle='dashed', color='gray', label='Now') # Current age of the universe
ax.axvline(0,    0, 1, linestyle='dotted', color='gray')              # World line for observer
# RHS
#ax.axvline(16.7, 0, 1, linestyle='dashed', color='r',    label='Event Horizon')
#ax.scatter(16.7, 13.8, marker='o', color='r', label='Event Horizon')
#ax.axvline(14.5, 0, 1, linestyle='dashed', color='green', label='Hubble Radius')
#ax.axvline(47,   0, 1, linestyle='dashed', color='b',    label='Particle Horizon')

# LHS
#ax.axvline(-16.7, 0, 1, linestyle='dashed', color='r',   ) 
#ax.axvline(-14.5, 0, 1, linestyle='dashed', color='green')
#ax.axvline(-47,   0, 1, linestyle='dashed', color='b',   ) 

# Plot event horizon in every time
def EventHorizon(amin=1e-30, amax=100):
    #a = np.linspace(amin,amax,100)
    a = np.logspace(np.log10(amin),np.log10(amax),int(5e2))
    eh = []
    t = []
    for i in a:
        ray = Photon(i, 1e5, 0, 1)
        dc = ray.d_photon[-1]/ray.a_photon[-1]*i
        eh.append(dc)
        t.append(ray.t_photon[0])
    return [np.array(t), np.array(eh)]
t, eh = EventHorizon()
ax.plot(eh/Gly2m, t/yr2s/1e9, color='r', linestyle='dashed', label='Event Horizon')
ax.plot(-eh/Gly2m, t/yr2s/1e9, color='r', linestyle='dashed')

ax.set_xlim([-50, 50])
ax.set_ylim([0, 75 ])
ax.set_xlabel('Proper Distance (Gly)')
ax.set_ylabel('Cosmic Time (Gyr)')
#plt.legend()
fig.savefig('World_Line_Cosmology.png', dpi=300, transparent=True)
plt.show()

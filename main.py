import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *
from functions import *
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Create figure
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(5,5))
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
plt.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
#ax.set_facecolor('#303030')

# This is particle horizon demo
Plot_ParticleHorizonDemo(fig, ax)

# This is event horizon demo
Plot_Multiple(fig, ax, Dcs=[14.5,16.7])

# This is for hubble radius demo

# This is for plotting light cone
Plot_LightCone(fig, ax)

# Current time and position
ax.axhline(13.8, 0, 1, linestyle='dashed', color='gray', label='Now') # Current age of the universe
ax.axvline(0,    0, 1, linestyle='dotted', color='gray')              # World line for observer

# Plot event horizon in every time
Plot_EventHorizon(fig, ax)

# Plot Hubble radius
Plot_HubbleRadius(fig, ax)

ax.set_xlim([-50, 50])
ax.set_ylim([0, 30 ])

ax.set_xlabel('Proper Distance (Gly)')
ax.set_ylabel('Cosmic Time (Gyr)')
plt.legend()
plt.tight_layout()
fig.savefig('Particle_Horizon_2.png', dpi=300, transparent=True)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import solve_ivp
from constants import *
from classes import *
from functions import *
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Create figure
plt.style.use('dark_background')
axcolor = 'lightgoldenrodyellow'
fig = plt.figure(figsize=(5,7))
ax  = plt.axes([0.1,0.3,0.8,0.65])
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
ax.set_facecolor('#303030')

# Making all the Sliders
axOL = plt.axes([0.1, 0.200, 0.8, 0.02], facecolor=axcolor)
axOM = plt.axes([0.1, 0.150, 0.8, 0.02], facecolor=axcolor)
axOR = plt.axes([0.1, 0.100, 0.8, 0.02], facecolor=axcolor)
sOL   = Slider(axOL, '$\Omega_\Lambda$', 0, 1, valinit=Omega_L)
sOM   = Slider(axOM, '$\Omega_m$',       0, 1, valinit=Omega_M)
sOR   = Slider(axOR, '$\Omega_R$',       0, 1, valinit=Omega_R)

# All the stuff to plot
def Plot_All(OM, OL, OR):
    global Omega_L
    global Omega_M
    global Omega_R
    Omega_L = OL
    Omega_M = OM
    Omega_R = OR
    ax.cla()
    # This is particle horizon demo
    Plot_ParticleHorizonDemo(fig, ax)

    # This is event horizon demo
    Plot_Multiple(fig, ax, Dcs=[14.5,16.7])

    # This is for plotting light cone
    Plot_LightCone(fig, ax)

    # Plot event horizon in every time
    Plot_EventHorizon(fig, ax)

    # Plot Hubble radius
    Plot_HubbleRadius(fig, ax)
    
    ax.set_xlim([-50, 50])
    ax.set_ylim([0, 30 ])
    # Current time and position
    ax.axhline(13.8, 0, 1, linestyle='dashed', color='gray', label='Now') # Current age of the universe
    ax.axvline(0,    0, 1, linestyle='dotted', color='gray')              # World line for observer

Plot_All(Omega_M, Omega_L, Omega_R)

# Changing parameters
def updatefig(val):
    Omega_L = sOL.val
    Omega_M = sOM.val
    Omega_R = sOR.val
    Plot_All(Omega_M, Omega_L, Omega_R)
    fig.canvas.draw_idle()

sOL.on_changed(updatefig)
sOM.on_changed(updatefig)
sOR.on_changed(updatefig)

ax.set_xlabel('Proper Distance (Gly)')
ax.set_ylabel('Cosmic Time (Gyr)')
#plt.legend()
#plt.tight_layout()
fig.savefig('Particle_Horizon_2.png', dpi=300, transparent=True)
plt.show()

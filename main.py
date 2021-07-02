import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

# Create figure
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(6,4))


ParaSet1 = []

A = Photon(1e-30, 1.5) 
ax.plot(A.t/yr2s/1e9, A.a, linestyle='solid', label='$\Omega_m={}$'.format(0.5))

ax.set_xlim([0, 16])
ax.set_ylim([0, 1.5])
ax.set_xlabel('Cosmic Time (Gyr)')
ax.set_ylabel('Size of the universe (scale factor)')
plt.legend()
fig.savefig('World_Line_Cosmology.png', dpi=300, transparent=True)
plt.show()

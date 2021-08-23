import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

# Create figure
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(6,4))

ParaSet1 = []

A = Photon(1e-50, 10) 

def LogLog(A=A, fig=fig, ax=ax):
    ax.loglog(A.t/yr2s/1e9, A.a, linestyle='solid', label='$\Lambda$CDM')

    ax.vlines(9.8,0,2,colors='gray',linestyles='dashed') # Dark Energy Dominated
    ax.vlines(0.047,0,2,colors='gray',linestyles='dashed') # Dark Energy Dominated

    ax.set_xlim([1e-10, 100])
    ax.set_ylim([1e-4, 10])
    ax.set_xlabel('Cosmic Time (Gyr)')
    ax.set_ylabel('Size of the universe (scale factor)')
    ax.scatter(13.8, 1, label='Now')
    plt.legend()
    fig.savefig('Scale_Factor_loglog.png', dpi=300, transparent=True)
    ax.cla()
    #plt.show()

def Normal_Plot(A=A, fig=fig, ax=ax):
    ax.plot(A.t/yr2s/1e9, A.a, linestyle='solid', label='$\Lambda$CDM')
    #ax.vlines(9.8,0,2,colors='gray',linestyles='dashed') # Dark Energy Dominated
    #ax.vlines(0.047,0,2,colors='gray',linestyles='dashed') # Dark Energy Dominated
    ax.set_xlim([0, 20])
    ax.set_ylim([0, 2])
    ax.set_xlabel('Cosmic Time (Gyr)')
    ax.set_ylabel('Size of the universe (scale factor)')
    ax.scatter(13.8, 1, label='Now')
    plt.legend()
    fig.savefig('Scale_Factor.png', dpi=300, transparent=True)
    ax.cla()
    #plt.show()

Normal_Plot()
LogLog()

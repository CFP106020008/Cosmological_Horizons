import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

def Plot_HubbleRadius(fig, ax, amin=1e-30, amax=100):
    a = np.logspace(np.log10(amin),np.log10(amax),int(1e3))
    def Hubble_Parameter(a):
        # Assuming k = 0, standard LCDM cosmology
        return H0*np.sqrt(Omega_R*a**-4 + Omega_M*a**-3 + Omega_L*a**0)
    H = Hubble_Parameter(a)
    d = c/H/Gly2m
    t = []
    for i in a:
        ray = Photon(i, 1e5, 0, 1)
        t.append(ray.t_photon[0])
    t = np.array(t)
    ax.plot(d, t/yr2s/1e9, color='g')
    ax.plot(-d, t/yr2s/1e9, color='g', label='Hubble Radius')
    return 

def Plot_EventHorizon(fig, ax, amin=1e-30, amax=100):
    def EventHorizon(amin=1e-30, amax=100):
        #a = np.linspace(amin,amax,100)
        a = np.logspace(np.log10(amin),np.log10(amax),int(1e3))
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
    return 

def Plot_ParticleHorizonDemo(fig, ax):
    ray1 = Photon(1e-50, 1, 47*Gly2m, -1) 
    ax.plot(ray1.d_photon/Gly2m,  ray1.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    ax.plot(ray1.d_source/Gly2m,  ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    ax.plot(-ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    ax.plot(-ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    #ax.axvline(47, 0, 1, linestyle='dashed', color='cyan', label='Particle Horizon')
    points = [[0, 0, 47, -47],[0, 13.8, 13.8, 13.8]]
    ax.scatter(points[0], points[1], color='cyan', zorder=10)
    return

def Plot_Multiple(fig, ax, Dcs=[1, 5, 10]):
    for Dc in Dcs:
        ray = Photon(1e-50, 20, Dc*Gly2m, -1) 
        ax.plot(ray.d_source/Gly2m, ray.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at {} Gly'.format(Dc)) # This is source galaxy
        ax.plot(-ray.d_source/Gly2m,ray.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    '''
    ray1 = Photon(1e-50, 20, 1*Gly2m, -1) 
    ray2 = Photon(1e-50, 20, 5*Gly2m, -1) 
    ray3 = Photon(1e-50, 20, 10*Gly2m, -1) 
     
    ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 1 Gly') # This is source galaxy
    ax.plot(ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 5 Gly') # This is source galaxy
    ax.plot(ray3.d_source/Gly2m, ray3.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at 10 Gly') # This is source galaxy
    ax.plot(-ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    ax.plot(-ray2.d_source/Gly2m, ray2.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    ax.plot(-ray3.d_source/Gly2m, ray3.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    '''
    return

def Plot_LightCone(fig, ax):
    ray1 = Photon(1e-50, 100, 47*Gly2m, -1, N=int(1e4)) # For the past
    ax.fill_betweenx(ray1.t_photon/yr2s/1e9, 
                     ray1.d_photon/Gly2m, 
                     -ray1.d_photon/Gly2m,
                     color='y',
                     alpha=0.5)
    ax.plot(ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    #ax.plot(ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    ax.plot(-ray1.d_photon/Gly2m, ray1.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    #ax.plot(-ray1.d_source/Gly2m, ray1.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    
    #ray2 = Photon(1, 100, 0, 1)
    #ax.plot(ray2.d_photon/Gly2m, ray2.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    #ax.plot(-ray2.d_photon/Gly2m, ray2.t_photon/yr2s/1e9, linestyle='solid',  color='cyan') # This is photon
    #ax.fill_betweenx(ray2.t_photon/yr2s/1e9, 
    #                 ray2.d_photon/Gly2m, 
    #                 -ray2.d_photon/Gly2m,
    #                 color='y',
    #                 alpha=0.5)
     
    return

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from constants import *
from classes import *

def Plot_HubbleRadius(fig, ax, amin=1e-30, amax=100):
    print('Plotting Hubble Radius')
    a = np.logspace(np.log10(amin),np.log10(amax),int(1e3))
    
    # Solving Hubble radius
    def Hubble_Parameter(a):
        # Assuming k = 0, standard LCDM cosmology
        return H0*np.sqrt(Omega_R*a**-4 + Omega_M*a**-3 + Omega_L*a**0)
    H = Hubble_Parameter(a)
    d = c/H/Gly2m
    
    # Solving cosmic time
    def cosmic_time(a, t):
        H = Hubble_Parameter(a) # Hubble Parameter
        return 1/(a*H)
    tmin = solve_ivp(       fun=cosmic_time, 
                            t_span=(1, amin),
                            y0=[138e8*yr2s]).y[0][-1]
    source_t = solve_ivp(   fun=cosmic_time, 
                            t_span=(amin, amax),
                            t_eval=a,
                            y0=[tmin])
    t = np.array(source_t.y[0])
    
    # Plot the results
    ax.plot(d,  t/yr2s/1e9, color='g')
    ax.plot(-d, t/yr2s/1e9, color='g', label='Hubble Radius')
    return 

def Plot_EventHorizon(fig, ax, amin=1e-50, amax=100, N=int(1e3)):
    print('Plotting Event Horizon')
    '''
    def EventHorizon(amin=1e-30, amax=100):
        #a = np.linspace(amin,amax,100)
        a = np.logspace(np.log10(amin),np.log10(amax),N)
        eh = []
        t = []
        for i in a:
            ray = Photon(i, 1e5, 0, 1)
            dc = ray.d_photon[-1]/ray.a_photon[-1]*i
            eh.append(dc)
            t.append(ray.t_photon[0])
        return [np.array(t), np.array(eh)]
    t, eh = EventHorizon()
    '''
    ray_EH = Photon(1, 1e5, 0, 1)
    EH_t0 = ray_EH.d_photon[-1]/ray_EH.a_photon[-1]
    ray_PH = Photon(amin, 1, 0, 1)
    PH_t0 = ray_PH.d_photon[-1]
    ray    = Photon(amin, amax, 0, 1)
    EH = (np.ones(N)*(EH_t0 + PH_t0) - ray.d_photon/ray.a_photon)*ray.a_photon
    t = ray.t_photon
    ax.plot(EH/Gly2m,  t/yr2s/1e9, color='r', linestyle='dashed', label='Event Horizon')
    ax.plot(-EH/Gly2m, t/yr2s/1e9, color='r', linestyle='dashed')
    return 

def Plot_ParticleHorizonDemo(fig, ax):
    print('Plotting Particle Horizon Demo')
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
    print('Plotting Multiple world lines')
    for Dc in Dcs:
        ray = Photon(1e-50, 20, Dc*Gly2m, -1) 
        ax.plot(ray.d_source/Gly2m, ray.t_source/yr2s/1e9, linestyle='dotted', color='cyan', label='Galaxy at {} Gly'.format(Dc)) # This is source galaxy
        ax.plot(-ray.d_source/Gly2m,ray.t_source/yr2s/1e9, linestyle='dotted', color='cyan') # This is source galaxy
    return

def Plot_LightCone(fig, ax):
    print('Plotting Past and Future Light cones')
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
    return

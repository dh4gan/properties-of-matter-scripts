# Written 7/7/17 by dh4gan
# Calculates the acceleration profile with time of a solar sail

import numpy as np
import matplotlib.pyplot as plt

sun_luminosity = 3.86 * 10**26  # [Watt] stellar luminosity
sun_radius = 695700000  # [m]
AU = 1.496e11 # [m]

rstar = 1.0*sun_radius
lstar = 1.0*sun_luminosity
c = 299792458  # [m/sec] speed of light
year = 3.15e7
month = year/12

msail = 0.001 # kg
area = 10.0 # square metres

deltat = 0.01*month

rinit = 5.0*sun_radius

tmax = 0.083*year
nsteps = int(tmax/deltat)+1

r = np.zeros(nsteps)
v = np.zeros(nsteps)
a = np.zeros(nsteps)
t = np.zeros(nsteps)

t[0] = 0.0
v[0] = 0.0
r[0] = rinit

for i in range (1,nsteps):

    # Radial photon pressure force
    F_photon_r = lstar * area / (3 * np.pi * c * rstar**2) * \
            (1 - (1 - (rstar / r[i-1])**2)**(1.5))
            
    print t[i-1]/year, r[i-1]/sun_radius,F_photon_r 
            
    a[i] = F_photon_r/msail
    
    
    v[i] = v[i-1] + a[i]*deltat
    r[i] = r[i-1] + v[i]*deltat
    t[i] = t[i-1]+deltat
    
    
fig1, axarr = plt.subplots(2,sharex=True)

axarr[1].set_xlabel('Time (months)', fontsize=20)
axarr[0].set_ylabel('Position (AU)', fontsize=18)
#ax1.tick_params('y', colors='red')
axarr[1].set_ylabel(r'Velocity (km s$^{-1}$)',fontsize=18)
#ax2.tick_params('y', colors='blue')

axarr[0].plot(t*12.0/year,r/AU, color='red', linewidth=2, )
axarr[1].plot(t*12.0/year,v/1000, color='blue', linewidth=2)
#axarr[2].plot(t/year,a/1000, color='green', linewidth=2)

label=r'Sail mass: '+str(msail)+' kg \nSail Area: '+str(area)+r' metres sq.'
axarr[0].annotate(label, xy=(0.6*12*tmax/year, 10.0),xycoords='data',size=13)

plt.show()
fig1.savefig('solar_sail_acceleration.png')
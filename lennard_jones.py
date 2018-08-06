import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#
# Script to integrate 1D motion in the Lennard-Jones potential
# We do this using odeint, which is a first order ODE solver
# We integrate the velocity and the position
# Our input is the matrix  Y=[v,x]
#

def dydt(y,t,sigma,epsilon):
    '''Given y=[v,x], compute dx/dt = v; and dv/dt=-grad(potential)'''
    v=y[0]
    x=y[1]
    dydt = [force(x,sigma,epsilon),v]

    return dydt

def lennard_jones(x,sigma,epsilon):
    '''Compute the 1D Lennard-Jones potential'''

    V = 4*epsilon*(np.power(sigma/x,12) - np.power(sigma/x,6))
    return V

def force(x,sigma,epsilon):
    '''Compute the force from the Lennard_Jones potential (-grad potential)'''

    F = 4*epsilon*(12*np.power(sigma/x,13) - 6*np.power(sigma/x,7))
    return F


#
# Set up Lennard-Jones potential constants
#

sigma = 1.5
epsilon = 100.0

# Potential minimum is at 2^{1/6}*sigma

#
# Setup initial conditions
#
v0 = 0 # initial speed
x0 = 3.0 # initial position

# Initial position = potential minimum
#x0 = np.power(2,1.0/6.0)*sigma

tmin = 0.0 # beginning time
tmax = 10.0 # end time
t = np.linspace(tmin,tmax, num=500) # Range of points between beginning and end time

# Compute potential at all x points
domain = np.linspace(0.1,5*sigma, num=100) # Range of position points
potential = lennard_jones(domain,sigma,epsilon)

# set up and perform integration
y0 = [v0,x0]
intgr = odeint(dydt, y0,t, args=(sigma,epsilon))

# extract velocity and position from integrator output
v = intgr[:,0]
x = intgr[:,1]

# Plot results
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
ax1.plot(x,t, label='position')
ax1.legend()

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_ylim((-2.0*epsilon, 10.0*epsilon))
ax2.set_xlabel('Position')
ax2.set_ylabel('Potential')
ax2.plot(domain,potential)


plt.show()










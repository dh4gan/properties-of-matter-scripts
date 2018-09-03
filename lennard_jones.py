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
    '''Compute the force from the Lennard-Jones potential (-grad potential)'''

    F = 4.0*epsilon*(12.0*np.power(sigma/x,13) - 6.0*np.power(sigma/x,7))
    return F

def rep_force(x,sigma,epsilon):
    '''Compute the repulsive force from the Lennard-Jones potential'''

    F = 4.0*epsilon*12.0*np.power(sigma/x,13)
    return F

def att_force(x,sigma,epsilon):
    '''Compute the force from the Lennard-Jones potential (-grad potential)'''
    
    F = -4*epsilon*6*np.power(sigma/x,7)
    return F


#
# Set up Lennard-Jones potential constants
#

sigma = 1.5
epsilon = 100.0

# Equilibrium separation = force is zero
equilibriumsep = np.power(2,1.0/6.0)*sigma

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

# Compute potential and force at all x points in the domain
domain = np.linspace(0.1,5*sigma, num=200) # Range of position points
potential_domain = lennard_jones(domain,sigma,epsilon)
force_domain = force(domain,sigma,epsilon)
rep_force_domain = rep_force(domain,sigma,epsilon)
att_force_domain = att_force(domain,sigma,epsilon)


# set up and perform integration
y0 = [v0,x0]
intgr = odeint(dydt, y0,t, args=(sigma,epsilon))

# extract velocity and position from integrator output
v = intgr[:,0]
x = intgr[:,1]

# Plot Position vs Time diagram for particle in potential
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
ax1.plot(x,t, label='position')
ax1.legend()


# Plot Lennard-Jones Potential
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_ylim((-2.0*epsilon, 10.0*epsilon))
ax2.set_xlabel('Position',fontsize=18)
ax2.set_ylabel('Lennard-Jones Potential',fontsize=18)
ax2.axhline(-epsilon, color='orange', linestyle='dashed', label='$\epsilon$')
ax2.axhline(0.0,color='black', linewidth=0.5)
ax2.axvline(equilibriumsep, color='green', linestyle='dotted', label='d')

ax2.plot(domain,potential_domain)
ax2.legend(fontsize=18)

fig2.savefig('LJ_potential.png')

# Plot Lennard-Jones Force
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.set_ylim((-4.0*epsilon,10.0*epsilon))
ax3.set_xlabel('Position',fontsize=18)
ax3.set_ylabel('Lennard-Jones Force',fontsize=18)
ax3.plot(domain,force_domain, label='Total Force')
ax3.axvline(equilibriumsep, color='green', linestyle='dotted', label='d')
ax3.legend(fontsize=18)

fig3.savefig('LJ_force.png')

# Plot attractive/repulsive components of Lennard-Jones Force
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
ax4.set_ylim((-1.5*epsilon,1.5*epsilon))
ax4.set_xlim(0.0,3.0*sigma)
ax4.set_xlabel('Position',fontsize=18)
ax4.set_ylabel('Repulsive/Attractive Force',fontsize=18)
ax4.plot(domain,rep_force_domain[rep_force_domain>0.0],label='Repulsive')
ax4.axhline(0.0, color='black', linestyle='dotted')
ax4.plot(domain,att_force_domain,label='Attractive')
ax4.legend(fontsize=18, loc='lower left')

fig4.savefig('LJ_forcecomponents.png')

plt.show()










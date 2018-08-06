import numpy as np
import matplotlib.pyplot as plt

#
# This Python code generates random walk trajectories and plots them
# Trajectories move (+-1, +-1) in 2D space, with equal probability of +- for
# each co-ordinate.
# 


def get_new_position(x):
    '''Given a co-ordinate x, make a random jump (+-1)'''

    jump = 1
    rand = np.random.random()
    if(rand>0.5):
        jump = -1

    return x+jump


def do_random_walk(nsteps,x0,y0):
    '''Generate a set of (x,y) points on a random walk '''
    x = np.zeros(nsteps)
    y = np.zeros(nsteps)

    x[0] = x0
    y[0] = y0

    for i in range(1,nsteps):
        x[i] = get_new_position(x[i-1])
        y[i] = get_new_position(y[i-1])

    return x,y


def do_n_random_walks(nsteps, nwalks, x0,y0):
    ''' Generate N sets of (x,y) points on a random walk'''

    x = np.zeros((nwalks,nsteps))
    y = np.zeros((nwalks,nsteps))

    # Generate nwalks random walks, and store in x,y arrays
    for i in range(nwalks):
        xi,yi = do_random_walk(nsteps,x0,y0)
        x[i,:] = xi
        y[i,:] = yi


    # Find maximum, minimum x and y
    xmax = np.amax(x)
    xmin = np.amin(x)

    ymax = np.amax(y)
    ymin = np.amin(y)

    return x,y, xmin,xmax,ymin,ymax


#
# Setup conditions
#

nsteps = 10000
nwalks = 10

x0 = 0
y0 = 0


#
# Generate random walks
#
x,y, xmin, xmax,ymin,ymax = do_n_random_walks(nsteps,nwalks,x0,y0)



#
# Plot as a figure
#
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)

for i in range(nwalks):
    ax1.plot(x[i,:],y[i,:])

plt.show()






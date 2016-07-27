# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 15:55:37 2016

@author: clay_budin
"""
# from: https://github.com/barbagroup/CFDPython/blob/master/lessons/18_Burgers_equation.ipynb

#%matplotlib inline
import numpy
from matplotlib import pyplot
from matplotlib import animation
#from JSAnimation.IPython_display import display_animation



#Basic initial condition parameters
#defining grid size, time steps, CFL condition, etc...
nx = 81
nt = 70
sigma = 1
dx = 4.0/nx
dt = sigma*dx

#Define a quick function to set up the initial square wave condition
def u_ic():
    u = numpy.ones(nx)
    u[(nx-1)/2:]=0
    return u

#Define two lambda functions to help with common operations encountered
#with the Euler equations
#utoE = lambda u: (u/2)**2
utoA = lambda u: u**2		# this is A*u @ (n,i)
utoE = lambda u: (u**2)/2
#utoA = lambda u: u

def laxfriedrichs(u, nt, dt, dx):
    #initialize our results array with dimensions nt by nx
    un = numpy.zeros((nt,len(u)))
    #copy the initial u array into each row of our new array
    un[:,:] = u.copy()

    '''
    Now, for each timestep, we're going to calculate u^n+1,
    then set the value of u equal to u^n+1 so we can calculate
    the next iteration.  For every timestep, the entire vector
    u^n is saved in a single row of our results array un.
    '''
    for i in range(1,nt):
        E = utoE(u)
        un[i,1:-1] = .5*(u[2:]+u[:-2]) - dt/(2*dx)*(E[2:]-E[:-2])
        un[i,0] = 1
        u = un[i].copy()

    return un


u = u_ic()      #make sure that u is set to our expected initial conditions
un = laxfriedrichs(u,nt,dt,dx)

print str(un.shape)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

def animate(data):
    x = numpy.linspace(0,4,nx)
    y = data
    line.set_data(x,y)
    return line,

# CFL = 1.0
#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='loop')


# CFL = .5
dt = .5*dx
u = u_ic()  ##Reset our initial conditions (the square wave)
un = laxfriedrichs(u,nt,dt,dx)


#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')



def laxwendroff(u, nt, dt, dx):
    un = numpy.zeros((nt,len(u)))
    un[:] = u.copy()

    for i in range(1,nt):
        E = utoE(u)
        un[i,1:-1] = u[1:-1] - dt/(2*dx) * (E[2:]-E[:-2]) + dt**2/(4*dx**2) *\
        ((u[2:]+u[1:-1])*(E[2:]-E[1:-1]) -\
        (u[1:-1]+u[:-2])*(E[1:-1]-E[:-2]))
        un[i,0]=1

        u = un[i].copy()

    return un


u = u_ic()
sigma = 1
dt = sigma*dx
un = laxwendroff(u,nt,dt,dx)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

# CFL = 1.0
#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')


# CFL = .5
u = u_ic()
sigma = .5
dt = sigma*dx
un = laxwendroff(u,nt,dt,dx)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')


def maccormack(u, nt, dt, dx):
    un = numpy.zeros((nt,len(u)))
    ustar = numpy.empty_like(u)
    un[:] = u.copy()
    ustar = u.copy()

    for i in range(1,nt):
        E = utoE(u)
        ustar[:-1] = u[:-1] - dt/dx * (E[1:]-E[:-1])
        Estar = utoE(ustar)
        un[i,1:] = .5 * (u[1:]+ustar[1:] - dt/dx * (Estar[1:] - Estar[:-1]))
        u = un[i].copy()

    return un

u = u_ic()
sigma = 1
dt = sigma*dx

un = maccormack(u,nt,dt,dx)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')


u = u_ic()
sigma = 0.5
dt = sigma*dx

un = maccormack(u,nt,dt,dx)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')




from scipy import linalg

def beamwarming(u, nt, dt, dx):
    ##Tridiagonal setup##
    a = numpy.zeros_like(u)
    b = numpy.ones_like(u)
    c = numpy.zeros_like(u)
    d = numpy.zeros_like(u)

    un = numpy.zeros((nt,len(u)))
    un[:]=u.copy()

    for n in range(1,nt):
        u[0] = 1
        E = utoE(u)
        au = utoA(u)

        a[0] = -dt/(4*dx)*u[0]
        a[1:] = -dt/(4*dx)*u[0:-1]
        a[-1] = -dt/(4*dx)*u[-1]

        #b is all ones

        c[:-1] = dt/(4*dx)*u[1:]

        print str(-.5*dt/dx*(E[2:]-E[0:-2])+dt/(4*dx)*(au[2:]-au[:-2]))
        #print str(dt/(4*dx)*(au[2:]-au[:-2]))
        d[1:-1] = u[1:-1]-.5*dt/dx*(E[2:]-E[0:-2])+dt/(4*dx)*(au[2:]-au[:-2])

        ###subtract a[0]*LHS B.C to 'fix' thomas algorithm
        d[0] = u[0] - .5*dt/dx*(E[1]-E[0])+dt/(4*dx)*(au[1]-au[0]) - a[0]

        ab = numpy.matrix([c,b,a])
        u = linalg.solve_banded((1,1), ab, d)
        u[0]=1
        un[n] = u.copy()
    return un


u = u_ic()
sigma = .5
dt = sigma*dx
nt=60
nt = 40

un = beamwarming(u,nt,dt,dx)

fig = pyplot.figure();
ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
line, = ax.plot([],[],lw=2);

anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
#display_animation(anim, default_mode='once')


'''


def dampit(u,eps,dt,dx):
	d = u[2]-.5*dt/dx*(u[3]**2/2-u[1]**2/2)+dt/(4*dx)*(u[3]**2-u[1]**2)\
		-eps*(u[4]-4*u[3]+6*u[2]-4*u[1]+u[0])
	return d

def beamwarming_damp(u, nt, dt, dx):
    ##Tridiagonal setup##
    a = numpy.zeros_like(u)
    b = numpy.ones_like(u)
    c = numpy.zeros_like(u)
    d = numpy.zeros_like(u)

    un = numpy.zeros((nt,len(u)))
    un[:] = u.copy()

    eps = .125

    for n in range(1,nt):
        u[0] = 1
        E = utoE(u)
        au = utoA(u)

        a[0] = -dt/(4*dx)*u[0]
        a[1:] = -dt/(4*dx)*u[0:-1]
        a[-1] = -dt/(4*dx)*u[-1]

        #b is all ones

        c[:-1] = dt/(4*dx)*u[1:]

        ###Calculate the damping factor for MOST of our u_vector
        d[2:-2] = u[2:-2]-.5*dt/dx*(E[3:-1]-E[1:-3])+dt/(4*dx)\
                            *(au[3:-1]-au[1:-3])\
                    -eps*(u[4:]-4*u[3:-1]+6*u[2:-2]-4*u[1:-3]+u[:-4])


        ###Calculate the damping factor for d[0] and d[1]
        damp = numpy.concatenate((np.ones(2), u[:3]))
        d[0] = dampit(damp,eps,dt,dx)
        damp = numpy.concatenate((np.ones(1), u[:4]))
        d[1] = dampit(damp,eps,dt,dx)

        ###subtract a[0]*LHS B.C to 'fix' thomas algorithm
        d[0] = d[0] - u[0] * a[0]

        ab = numpy.matrix([c,b,a])
        u = linalg.solve_banded((1,1), ab, d)
        u[0]=1
        un[n] = u.copy()
    return un


u = u_ic()
sigma = 0.5
dt = sigma*dx
nt = 120
un = beamwarming_damp(u,nt,dt,dx)

#fig = pyplot.figure();
#ax = pyplot.axes(xlim=(0,4),ylim=(-.5,2));
#line, = ax.plot([],[],lw=2);

#anim = animation.FuncAnimation(fig, animate, frames=un, interval=50)
display_animation(anim, default_mode='once')




from IPython.core.display import HTML
def css_styling():
    styles = open("../styles/custom.css", "r").read()
    return HTML(styles)
css_styling()

'''

pyplot.show()





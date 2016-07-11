# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 1: 1D linear convection
	du/dt + c*du/dx = 0
linear because c is constant

Finite Difference approximations to du/dx
Forward diff:  du/dx ~= (u(x+1) - u(x)) / delta_x        O(delta_x)
Backward diff: du/dx ~= (u(x) - u(x-1)) / delta_x        O(delta_x)
Central diff:  du/dx ~= (u(x+1) - u(x-1)) / (2 * delta_x)    O(delta_x^2)

Forward diff in time
Backward diff in space

Domain: [0,2]
Range: [1,2]
Initial Conditions: square wave, half-sine, full inverted cosine
Boundary Conditions: u = 1 @ x = 0,2

c*dt/dx seems to be important to stability
	if > 1, sim blows up, eg if 2: u(t+1,i) = 2*u(t,i-1) - u(t,i)
	if == 1, sim is perfect: u(t+1,i) = u(t,i-1)
	if < 1, sim works but seems to damp high freqs, eg if .5 get averaging bwtn u(t,i-1) and u(t,i)

NOTES:
	must use global declaration when changing a global, not accessing (reading)
	don't need global when changing elements of list
	need to use float constants where floats are needed - won't promote?
	yprev = ydata doesn't copy list, just pointer to list



"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1)
ax.grid()

# simulation constants
nx = 201
nt = 150
c = 1.0 #0.5 #1.0
dt = 0.01
dx = 2.0 / (nx-1.0)

# simlation grids - 1D
xdata, ydata, yprev = [], [], []


# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt) + " c*dt/dx = " + str(c*dt/dx)
	ct = 0

	ax.set_xlim(-.1, 2.1)
	ax.set_ylim(0.0,2.5)
	del xdata[:]
	del ydata[:]
	for i in range(nx):
		x = i*dx
		xdata.append(x)
		y = 1.0

		# initialize range [.5,1] with various test fns
		if x >= .5 and x <= 1.0: y = 2.0		# square wave
		#if x >= .5 and x <= 1.0: y = 1.0 + np.sin((x-.5)*np.pi*2.0)	# half-sine
#		if x >= .5 and x <= 1.0:				# full inverted cosine
#			#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#			y = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0))
#			#print "y = " + str(y)
		#if x >= .5 and x <= 1.0: y = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	# full-sine

		ydata.append(y)
		yprev.append(y)
	line.set_data(xdata, ydata)
	return line,


def data_gen ():
	global ct
	while ct <= nt:
		#print "data_gen() ct =", ct
		ct += 1

		# copy contents of ydata to yprev
		#yprev = ydata   # NO - only sets pointer
		#for i in range(nx): yprev[i] = ydata[i]
		yprev[:] = ydata[:]

		for i in range(1,nx-1):
			ydata[i] = yprev[i] - (c*dt/dx)*(yprev[i]-yprev[i-1])
			#ydata[i] = yprev[i] - (c*dt/(2.0*dx))*(yprev[i+1]-yprev[i-1])		# central diff - not stable for any c (?)
			#ydata[i] = yprev[i] + (c*dt/dx)*(yprev[i+1]-yprev[i])				# change sign, use forward diff - wave moves to left
		yield


def run (data):
	#print "run()"
	line.set_data(xdata, ydata)
	return line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

plt.show()




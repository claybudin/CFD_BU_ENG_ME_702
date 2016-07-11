# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 2: 1D non-linear convection
	du/dt + u*du/dx = 0
non-linear because u multiplies du/dx
also called Inviscid Burger's eq
Can generate discontinuous solutions from smooth Initial Conditions (shocks)

Forward diff in time (forward Euler)
Backward diff in space

Domain: [0,2]
Range: [1,2]
Initial Conditions: square wave, half-sine, full inverted cosine
Boundary Conditions: u = 1 @ x = 0,2

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1)
ax.grid()

# simulation constants
nx = 101
nt = 150
dt = 0.01
dx = 2.0 / (nx-1.0)

# simlation grids - 1D
xdata, ydata, yprev = [], [], []


# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt)
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
		#if x >= .5 and x <= 1.0: y = 2.0		# square wave
		#if x >= .5 and x <= 1.0: y = 1.0 + np.sin((x-.5)*np.pi*2.0)	# half-sine
#		if x >= .5 and x <= 1.0:				# full inverted cosine
#			#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#			y = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0))
#			#print "y = " + str(y)
		#y -= .5
		if x >= .5 and x <= 1.0: y = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	# full-sine

		ydata.append(y)
		yprev.append(y)
	line.set_data(xdata, ydata)
	return line,


def data_gen ():
	global ct
	while ct <= nt:
		#print "data_gen() ct =", ct
		ct += 1
		for i in range(nx): yprev[i] = ydata[i]
		for i in range(1,nx-1):
			ydata[i] = yprev[i] - (yprev[i]*dt/dx)*(yprev[i]-yprev[i-1])
		yield


def run (data):
	#print "run()"
	line.set_data(xdata, ydata)
	return line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

plt.show()




# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 4: 1D Burgers's Equation
	du/dt + u*du/dx - v*d^2u/dx^2 = 0

Combines steps 2 and 3 - convection and diffusion
Has analytical solution (for some ICs?)

Domain: [0, 2*pi]
Range: ~[.9,7.1]
Initial Conditions: u(x) = 4 - 2*v*(dp/dx / p)
	where p = exp(-x^2/4v) + exp(-(x-2*pi)^2/4v)
	dp/dx = -x/2v*exp(-x^2/4v) + -(x-2*pi)/2v*exp(-(x-2*pi)^2/4v)
Boundary Conditions: u(0) = u(2*pi) - periodic BC

Analytical Solution for IC:
	u(x,t) = 4 - 2*v*(dp/dx / p)
	where p = exp(-(x-4t)^2/4v(t+1)) + exp(-(x-4t-2*pi)^2/4v(t+1))

TO DO:
	+enforce periodic BC
	+plot analytic soln next to numerical soln
	create guthub repo and check everything in
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1)
lineAna, = ax.plot([], [], lw=1)

ax.grid()
ax.set_xlim(-.1, 6.4)
ax.set_ylim(.7, 7.3)
ax.set_title("Step 4: Burgers's Equation\nBlue - Numerical, Green - Analytical")

# simulation constants
nx = 201
nt = 1000
dt = 0.002
dx = (2.0*np.pi) / nx   # range is 0 to 2pi-dx to enforce periodic BC u(2pi) = u(0)
v = .02

# simlation grids - 1D
xdata, ydata, yprev, yana = [], [], [], []


# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt) + " v*dt/(dx*dx) = " + str(v*dt/(dx*dx))
	ct = 0

	del xdata[:]
	del ydata[:]
	del yprev[:]
	del yana[:]

	ymin = 1000000000.0
	ymax = -100000000.0
	for i in range(nx):
		x = i*dx
		#print x
		xdata.append(x)

		p1 = np.exp(-(x*x)/(4.0*v))
		p2 = np.exp(-(x-2*np.pi)*(x-2*np.pi)/(4.0*v))
		y = 4.0 - 2.0*v*(((-x/(2.0*v))*p1+(-(x-2.0*np.pi)/(2.0*v))*p2) / (p1+p2))

		if y < ymin: ymin = y
		if y > ymax: ymax = y

		ydata.append(y)
		yprev.append(y)
		yana.append(y)
	print "ymin = " + str(ymin) + " ymax = " + str(ymax)
	line.set_data(xdata, ydata)
	lineAna.set_data(xdata, yana)
	#return line,


def data_gen ():
	global ct
	while ct <= nt:
		#print "data_gen() ct =", ct
		ct += 1

		for i in range(nx): yprev[i] = ydata[i]

		for i in range(nx):
			# enforce BC by wrapping indices
			im1 = i - 1
			if (im1 < 0): im1 = nx-1
			ip1 = i + 1
			if (ip1 >= nx): ip1 = 0
			ydata[i] = yprev[i] - (yprev[i]*dt/dx)*(yprev[i]-yprev[im1]) + (v*dt/(dx*dx))*(yprev[ip1]-2.0*yprev[i]+yprev[im1])

		# analytical soln
		t = ct * dt
		for i in range(nx):
			x = i*dx
			p1 = np.exp(-((x-4*t)*(x-4*t))/(4.0*v*(t+1)))
			p2 = np.exp(-(x-4*t-2*np.pi)*(x-4*t-2*np.pi)/(4.0*v*(t+1)))
			y = 4.0 - 2.0*v*(((-(x-4*t)/(2.0*v*(t+1)))*p1+(-(x-4*t-2.0*np.pi)/(2.0*v*(t+1)))*p2) / (p1+p2))
			yana[i] = y

		yield


def run (data):
	#print "run()"
	line.set_data(xdata, ydata)
	lineAna.set_data(xdata, yana)
	#return line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

plt.show()




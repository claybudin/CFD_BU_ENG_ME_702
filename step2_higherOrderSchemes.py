# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 2: 1D non-linear convection (Inviscid Burgers's eqn)
	du/dt + u*du/dx = 0
	du/dt + d/dx(u^2/2) = 0   - "conservative" form
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
ax.grid()
ax.set_xlim(-.1, 4.1)
ax.set_ylim(-.1, 1.1)


# plots
line_bd, = ax.plot([], [], 'b', lw=1)
line_fd, = ax.plot([], [], 'r', lw=1)
line_lf, = ax.plot([], [], 'g', lw=1)
line_lw, = ax.plot([], [], 'y', lw=1)



# simulation constants
nx = 201
nt = 200
dt = 0.02 #0.01
dx = 4.0 / (nx-1)


# simlation grids - 1D
xdata = []
ydbd = []
ydfd = []
ydlf = []
ydlw = []


# global var - current time step
ct = 0


def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt) + " dt/dx = " + str(dt/dx)
	ct = 0

	del xdata[:]
	del ydbd[:]
	del ydfd[:]
	del ydlf[:]
	del ydlw[:]

	for i in range(nx):
		x = i*dx
		xdata.append(x)

		y = 0.0
		if x <= 2.0: y = 1.0

		ydbd.append(y)
		ydfd.append(y)
		ydlf.append(y)
		ydlw.append(y)

	line_bd.set_data(xdata, ydbd)
	line_fd.set_data(xdata, ydfd)
	line_lf.set_data(xdata, ydlf)
	line_lw.set_data(xdata, ydlw)



def data_gen ():
	global ct
	while ct < nt:
		#print "data_gen() ct =", ct
		ct += 1

		# backwards difference scheme
		# nothing happens here with these ICs
		ypbd = ydbd[:]
		for i in range(1,nx-1):
			ydbd[i] = ypbd[i] - (ypbd[i]*dt/dx)*(ypbd[i]-ypbd[i-1])

		# forwards difference scheme
		# unstable
#		ypfd = ydfd[:]
#		for i in range(1,nx-1):
#			ydfd[i] = ypfd[i] - (ypfd[i]*dt/dx)*(ypfd[i+1]-ypfd[i])

		# Lax-Friedrichs: u(n+1,i) = 1/2*(u(n,i+1)+u(n,i-1)) - u(n,i)*dt/(2*dx) * (u(n,i+1)-u(n,i-1))
		yplf = ydlf[:]
		for i in range(1,nx-1):
			ydlf[i] = .5*(yplf[i+1]+yplf[i-1]) - (yplf[i]*dt/(2*dx))*(yplf[i+1]-yplf[i-1])


		# Lax-Wendroff: u(n+1,i) = u(n,i) - sigma/2*(u(n,i+1)-u(n,i-1)) + sigma^2/2*(u(n,i+1)-2*u(n,i)+u(n,i-1))
		# do I need to re-do derivations because substitution of d/dx for d/dt has changed?
		yplw = ydlw[:]
		for i in range(1,nx-1):
			# this does LW same as linear case, just substituting u for c
			# doesn't move at all
			#ydlw[i] = yplw[i] - (yplw[i]*dt/(2*dx))*(yplw[i+1]-yplw[i-1]) + (yplw[i]*dt/dx)*(yplw[i]*dt/dx)/2*(yplw[i+1]-2*yplw[i]+yplw[i-1])
			# rederived LW from Taylor series using u*du/dx instead of c*du/dx
			# doesn't move either
			ydlw[i] = yplw[i] - (yplw[i]*dt/(2*dx))*(yplw[i+1]-yplw[i-1]) + \
					  (dt/dx)*(dt/dx)/2*((yplw[i]*yplw[i]*(yplw[i+1]-2*yplw[i]+yplw[i-1])) + .5*yplw[i]*(yplw[i+1]-yplw[i-1]))

		yield


def run (data):
	#print "run()"
	line_bd.set_ydata(ydbd)
	#line_fd.set_ydata(ydfd)
	line_lf.set_ydata(ydlf)
	line_lw.set_ydata(ydlw)


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)
plt.show()




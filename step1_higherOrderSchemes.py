# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 1: 1D linear convection
	du/dt + c*du/dx = 0
linear because c is constant

Schemes:
Upwind - Forward diff in time, Backward diff in space
Leapfrog - Centered diff in time, CD in space
	u(n+1,i) = u(n-1,i) - c*dt/dx * (u(n,i+1)-u(n,i-1)) - centered diff using point from 2 time steps back
Lax-Friedrichs - FD in time, CD in space using avg for u(n,i)
	u(n+1,i) = 1/2*(u(n,i+1)+u(n,i-1)) - c*dt/(2*dx) * (u(n,i+1)-u(n,i-1))
Lax-Wendroff - FD in time, CD in space
	u(n+1,i) = u(n,i) - sigma/2*(u(n,i+1)-u(n,i-1)) + sigma^2/2*(u(n,i+1)-2*u(n,i)+u(n,i-1))
	adds second derivive term for more accuracy
	I think Lax-Wendroff applies specifically to Convection, since it uses a substitution based on convection eq (?)
Richtmyer/Lax-Wendroff - multi-step (Predictor-Corrector)
	Variant 1: Step 1 - Lax-Friedrich to time n+1/2, Step 2 - Leap frog using values from Step 1
		Step 1:  u(n+1/2, i) = 1/2*(u(n,i+1)+u(n,i-1)) - c*dt/(4*dx)*(u(n,i+1)-u(n,i-1))
		Step 2:  u(n+1,i) = u(n,i) - c*dt/(2*dx)*(u(n+1/2, i+1) - u(n+1/2, i-1))
	Variant 2: Step 1 - Lax-Friedrich to time n+1/2 and space i+1/2, Step 2 - Leap frog using values from Step 1
		Step 1:  u(n+1/2, i+1/2) = 1/2*(u(n,i+1)+u(n,i)) - c*dt/(2*dx)*(u(n,i+1)-u(n,i))
		Step 2:  u(n+1,i) = u(n,i) - c*dt/dx*(u(n+1/2, i+1/2) - u(n+1/2, i-1/2))
MacCormack - multi-step - creates intermediate values u*(x)
	Step 1 - FD: u*(i) = u(n,i) - c*dt/dx*(u(n,i+1)-u(n,i))
	Step 2 - BD: u(n+1,i) = .5*(u(n,i)+u*(i) - c*dt/dx*(u*(i)-u*(i-1)))
	Can also alternate FD and BD in alternate steps
Beam-Warming - 2nd order BD, gets u(n+1,i) using u(n,i), u(n,i-1) and u(n,i-2)

Domain: [0,2]
Range: [0,1]
Initial Conditions: u(0) = 1, 0 everywhere else
Boundary Conditions: u(0) = 1

c*dt/dx (= sigma) seems to be important to stability - CFL number
	if > 1, sim blows up, eg if 2: u(t+1,i) = 2*u(t,i-1) - u(t,i)
	if == 1, sim is perfect: u(t+1,i) = u(t,i-1)
	if < 1, sim works but seems to damp high freqs, eg if .5 get averaging bwtn u(t,i-1) and u(t,i) - numerical diffusion



"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig, ax = plt.subplots()
ax.grid()
ax.set_xlim(-.1, 2.1)
ax.set_ylim(-.1, 1.4)

line_uw, = ax.plot([], [], "b", lw=1)
line_lf, = ax.plot([], [], "r", lw=1)
line_lxfr, = ax.plot([], [], "g", lw=1)
line_lw, = ax.plot([], [], "c", lw=1)
line_rlw1, = ax.plot([], [], "y", lw=1)
line_rlw2, = ax.plot([], [], "m", lw=1)
line_mc1, = ax.plot([], [], "y", lw=1)
line_mc2, = ax.plot([], [], "g", lw=1)
line_bw, = ax.plot([], [], "r", lw=1)



# simulation constants
nx = 201
nt = 300
c = 0.5 #1.0
dt = 0.01
dx = 2.0 / (nx-1.0)

# simlation grids - 1D
xdata = []
yduw = []				# upwind scheme - backwards difference
ydlf, yplf = [], []		# leapfrog scheme - save previous data since we need u(t-2)
ydlxfr = []
ydlw = []
ydrlw1 = []
ydrlw2 = []
ydmc1 = []
ydmc2 = []
ydbw = []


# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt) + " c*dt/dx = " + str(c*dt/dx)
	ct = 0

	del xdata[:]
	del yduw[:]
	del ydlf[:]
	del yplf[:]
	del ydlxfr[:]
	del ydlw[:]
	del ydrlw1[:]
	del ydrlw2[:]
	del ydmc1[:]
	del ydmc2[:]
	del ydbw[:]

	for i in range(nx):
		x = i*dx
		xdata.append(x)
		y = 0.0

		# IC u(0) = 1
		if i == 0: y = 1.0
		if i == 1: y = 1.0		# need this for Beam-Warming scheme since it starts at u[2]

		yduw.append(y)
		ydlf.append(y)
		yplf.append(y)
		ydlxfr.append(y)
		ydlw.append(y)
		ydrlw1.append(y)
		ydrlw2.append(y)
		ydmc1.append(y)
		ydmc2.append(y)
		ydbw.append(y)

	line_uw.set_data(xdata, yduw)
	line_lf.set_data(xdata, ydlf)
	line_lxfr.set_data(xdata, ydlxfr)
	line_lw.set_data(xdata, ydlw)
	line_rlw1.set_data(xdata, ydrlw1)
	line_rlw2.set_data(xdata, ydrlw2)
	line_mc1.set_data(xdata, ydmc1)
	line_mc2.set_data(xdata, ydmc2)
	line_bw.set_data(xdata, ydbw)




def data_gen ():
	global ct
	while ct < nt:
		#print "data_gen() ct =", ct
		ct += 1

		# for useful visual comparissons among these schemes, see:
		#	http://www.thevisualroom.com/heavy_side_and_sinusoidal_input.html

		# Upwind scheme (backwards diff)
		# copy contents of yduw to ypuw
		ypuw = []
		ypuw[:] = yduw[:]

		for i in range(1,nx-1):
			yduw[i] = ypuw[i] - (c*dt/dx)*(ypuw[i]-ypuw[i-1])
			#yduw[i] = ypuw[i] - (c*dt/(2.0*dx))*(ypuw[i+1]-ypuw[i-1])		# central diff - not stable for any c


		# leapfrog scheme
		# stable for sigma <= 1.0
		# better than upwind on leading edge, but trailing edge has oscillations
		ypplf = []
		ypplf[:] = yplf[:]
		yplf[:] = ydlf[:]
		if ct == 1:
			# need a startup scheme - upwind
			#print "Startup"
			for i in range(1,nx-1):
				ydlf[i] = yplf[i] - (c*dt/dx)*(yplf[i]-yplf[i-1])
		else:
			# all other times use leapfrog update step
			for i in range(1,nx-1):
				ydlf[i] = ypplf[i] - (c*dt/dx)*(yplf[i+1]-yplf[i-1])


		# Lax-Friedrichs
		# seems like worse numerical diffusion than with Upwind
		yplxfr = []
		yplxfr[:] = ydlxfr[:]
		for i in range(1,nx-1):
			ydlxfr[i] = .5*(yplxfr[i+1]+yplxfr[i-1]) - (c*dt/(2*dx))*(yplxfr[i+1]-yplxfr[i-1])

		# Lax-Wendroff
		# improved on leading edge like leapfrog, oscillations on trailing edge but better than leapfrog
		yplw = []
		yplw[:] = ydlw[:]
		for i in range(1,nx-1):
			ydlw[i] = yplw[i] - (c*dt/(2*dx))*(yplw[i+1]-yplw[i-1]) + (c*dt/dx)*(c*dt/dx)*.5*(yplw[i+1]-2*yplw[i]+yplw[i-1])


		# Richtmyer/Lax-Wendroff variant 1
		# seems less accurate than just Lax-Wendroff
		# step 1
		yhrlw1 = []
		yhrlw1.append(ydrlw1[0])
		for i in range(1,nx-1):
			yhrlw1.append(.5*(ydrlw1[i+1]+ydrlw1[i-1]) - (c*dt/(4*dx))*(ydrlw1[i+1]-ydrlw1[i-1]))
		yhrlw1.append(ydrlw1[nx-1])
		# step 2
		# can just reuse ydrlw1 since we only need to reference point at i from previous time
		for i in range(1,nx-1):
			ydrlw1[i] = ydrlw1[i] - (c*dt/(2*dx))*(yhrlw1[i+1]-yhrlw1[i-1])

		# Richtmyer/Lax-Wendroff variant 2
		# seems identical to Lax-Wendroff - true for linear PDEs
		# step 1
		yhrlw2 = []
		for i in range(nx-1):	# yhrlw2[i] hold values for i+1/2, need to compute yhrlw2[0]
			yhrlw2.append(.5*(ydrlw2[i+1]+ydrlw2[i]) - (c*dt/(2*dx))*(ydrlw2[i+1]-ydrlw2[i]))
		# step 2
		for i in range(1,nx-1):
			ydrlw2[i] = ydrlw2[i] - (c*dt/dx)*(yhrlw2[i]-yhrlw2[i-1])

		# MacCormack
		# also equivalent to Lax-Wendoff for linear PDEs
		# step 1
		ymcus1 = []
		for i in range(nx-1):	# compute u*(0)
			ymcus1.append(ydmc1[i] - (c*dt/dx)*(ydmc1[i+1]-ydmc1[i]))
		# step 2
		for i in range(1,nx-1):
			ydmc1[i] = .5*(ydmc1[i] + ymcus1[i] - (c*dt/dx)*(ymcus1[i]-ymcus1[i-1]))

		# MacCormack variation - alternate order of FD and BD in steps on alternate iterations
		# also equivalent to Lax-Wendoff for linear PDEs
		ymcus2 = []
		if ct % 2 == 1:
			for i in range(nx-1):	# compute u*(0)
				ymcus2.append(ydmc2[i] - (c*dt/dx)*(ydmc2[i+1]-ydmc2[i]))
			for i in range(1,nx-1):
				ydmc2[i] = .5*(ydmc2[i] + ymcus2[i] - (c*dt/dx)*(ymcus2[i]-ymcus2[i-1]))
		else:
			ymcus2.append(0)			# fill u*(0) with dummy value - not used
			for i in range(1,nx):	# compute u*(nx-1)
				ymcus2.append(ydmc2[i] - (c*dt/dx)*(ydmc2[i]-ydmc2[i-1]))
			for i in range(1,nx-1):
				ydmc2[i] = .5*(ydmc2[i] + ymcus2[i] - (c*dt/dx)*(ymcus2[i+1]-ymcus2[i]))


		# Beam-Warming - similar to LW but uses 2nd order BDs
		# oscillations lead rather than trail
		ypbw = ydbw[:]
		for i in range(2,nx):
			ydbw[i] = ypbw[i] - (c*dt/(2*dx))*(3*ypbw[i]-4*ypbw[i-1]+ypbw[i-2]) + (c*dt/dx)*(c*dt/dx)*.5*(ypbw[i]-2*ypbw[i-1]+ypbw[i-2])

		# don't need to enfore BC since leftmost element of ydata arrays not touched by iterations

		yield


def run (data):
	#print "run()"
	line_uw.set_ydata(yduw)
	#line_lf.set_ydata(ydlf)
	#line_lxfr.set_ydata(ydlxfr)
	line_lw.set_ydata(ydlw)
	#line_rlw1.set_ydata(ydrlw1)
	line_rlw2.set_ydata(ydrlw2)
	line_mc1.set_ydata(ydmc1)
	line_mc2.set_ydata(ydmc2)
	line_bw.set_ydata(ydbw)



ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

plt.show()




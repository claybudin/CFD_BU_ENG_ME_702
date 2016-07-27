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
ax.set_ylim(-.1, 1.5)


# plots
line_bd, = ax.plot([], [], 'b', lw=1)
line_bdc, = ax.plot([], [], 'r', lw=1)
line_fd, = ax.plot([], [], 'g', lw=1)
line_lf, = ax.plot([], [], 'c', lw=1)
line_lfc, = ax.plot([], [], 'm', lw=1)
line_lw, = ax.plot([], [], 'y', lw=1)
line_mc, = ax.plot([], [], 'g', lw=1)
line_bw, = ax.plot([], [], 'b', lw=1)
line_bw2, = ax.plot([], [], 'r', lw=1)
line_bw3, = ax.plot([], [], 'g', lw=1)



# simulation constants
nx = 201
nt = 200
dt = 0.02 #0.01
dx = 4.0 / (nx-1)


# simlation grids - 1D
xdata = []
ydbd = []
ydbdc = []
ydfd = []
ydlf = []
ydlfc = []
ydlw = []
ydmc = []
ydbw = []
ydbw2 = []
ydbw3 = []


# global var - current time step
ct = 0


def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt) + " dt/dx = " + str(dt/dx)
	ct = 0

	del xdata[:]
	del ydbd[:]
	del ydbdc[:]
	del ydfd[:]
	del ydlf[:]
	del ydlfc[:]
	del ydlw[:]
	del ydmc[:]

	# because of the tridiagonal solver, ydbw is no longer a simple array that can be deleted like the others
	global ydbw, ydbw2, ydbw3
	del ydbw
	del ydbw2
	del ydbw3
	ydbw = []
	ydbw2 = []
	ydbw3 = []


	for i in xrange(nx):
		x = i*dx
		xdata.append(x)

		y = 0.0
		if x <= 2.0: y = 1.0

		ydbd.append(y)
		ydbdc.append(y)
		ydfd.append(y)
		ydlf.append(y)
		ydlfc.append(y)
		ydlw.append(y)
		ydmc.append(y)
		ydbw.append(y)
		ydbw2.append(y)
		ydbw3.append(y)


	line_bd.set_data(xdata, ydbd)
	line_bdc.set_data(xdata, ydbdc)
	#line_fd.set_data(xdata, ydfd)
	line_lf.set_data(xdata, ydlf)
	line_lfc.set_data(xdata, ydlfc)
	line_lw.set_data(xdata, ydlw)
	line_mc.set_data(xdata, ydmc)
	line_bw.set_data(xdata, ydbw)
	line_bw2.set_data(xdata, ydbw2)
	line_bw3.set_data(xdata, ydbw3)


# https://gist.github.com/ofan666/1875903
## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver (a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = ac
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    del bc, cc, dc  # delete variables from memory

    return xc



def data_gen ():
	global ct
	while ct < nt:
		#print "data_gen() ct =", ct
		ct += 1

		# backwards difference scheme
		# nothing happens here with these ICs
		ypbd = ydbd[:]
		for i in xrange(1,nx-1):
			ydbd[i] = ypbd[i] - (ypbd[i]*dt/dx)*(ypbd[i]-ypbd[i-1])

		# backwards difference scheme - conservative
		# 	du/dt + df/dx = 0,  f = 1/2*u^2, df/dx = u*du/dx (by the chain rule)
		# approximate df/dx with BD: (f(i)-f(i-1)) / dx, du/dt with FD and solve for u(n+1,i)
		# this works - seems better than LF or LW for these ICs
		ypbdc = ydbdc[:]
		for i in xrange(1,nx-1):
			ydbdc[i] = ypbdc[i] - (dt/(2*dx))*(ypbdc[i]*ypbdc[i]-ypbdc[i-1]*ypbdc[i-1])

		# forwards difference scheme
		# unstable
#		ypfd = ydfd[:]
#		for i in xrange(1,nx-1):
#			ydfd[i] = ypfd[i] - (ypfd[i]*dt/dx)*(ypfd[i+1]-ypfd[i])

		# see if FD works better using conservative formulation
		# still unstable
#		ypfd = ydfd[:]
#		for i in xrange(1,nx-1):
#			ydfd[i] = ypfd[i] - (dt/(2*dx))*(ypfd[i+1]*ypfd[i+1]-ypfd[i]*ypfd[i])

		# try CD using conservative formulation
		# still unstable
#		ypfd = ydfd[:]
#		for i in xrange(1,nx-1):
#			ydfd[i] = ypfd[i] - (dt/(4*dx))*(ypfd[i+1]*ypfd[i+1]-ypfd[i-1]*ypfd[i-1])



		# Lax-Friedrichs: u(n+1,i) = 1/2*(u(n,i+1)+u(n,i-1)) - u(n,i)*dt/(2*dx) * (u(n,i+1)-u(n,i-1))
		# works, some diffusion
		yplf = ydlf[:]
		for i in xrange(1,nx-1):
			ydlf[i] = .5*(yplf[i+1]+yplf[i-1]) - (yplf[i]*dt/(2*dx))*(yplf[i+1]-yplf[i-1])

		# conservative version
		yplfc = ydlfc[:]
		for i in xrange(1,nx-1):
			ydlfc[i] = .5*(yplfc[i+1]+yplfc[i-1]) - (dt/(4*dx))*(yplfc[i+1]*yplfc[i+1]-yplfc[i-1]*yplfc[i-1])


		# Lax-Wendroff
		# works using conservative formulation - some trailing oscillations
		yplw = ydlw[:]
		for i in xrange(1,nx-1):
			# the standard derivation seems to use the conservative form
			# from several sources, eg:
			# http://www.bcamath.org/projects/NUMERIWAVES/Burgers_Equation_M_Landajuela.pdf
			# https://people.sc.fsu.edu/~jburkardt/m_src/fd1d_burgers_lax/fd1d_burgers_lax.m
			# note that PDF file lists first term on RHS as u(n,i+1) which fails, while MatLab program lists it as u(n,i) which works
			# this appears to be confirmed in http://www.mat.univie.ac.at/~obertsch/literatur/burgers_equation.pdf
			ydlw[i] = yplw[i] - dt/(4*dx)*(yplw[i+1]*yplw[i+1]-yplw[i-1]*yplw[i-1]) + \
						dt*dt/(8*dx*dx)*((yplw[i]+yplw[i+1])*(yplw[i+1]*yplw[i+1]-yplw[i]*yplw[i]) - \
										(yplw[i]+yplw[i-1])*(yplw[i]*yplw[i]-yplw[i-1]*yplw[i-1]))


		# MacCormack - from lbarba's notes in "assignmentBurger.pdf"
		yusmc = ydmc[:]
		for i in xrange(1,nx-1):
			yusmc[i] = ydmc[i] - dt/(2*dx)*(ydmc[i+1]*ydmc[i+1]-ydmc[i]*ydmc[i])
		for i in xrange(1,nx-1):
			ydmc[i] = .5*(ydmc[i]+yusmc[i]-dt/(2*dx)*(yusmc[i]*yusmc[i]-yusmc[i-1]*yusmc[i-1]))


		# implicit Beam-Warming - solve tri-diagonal matrix
		# set up arrays
		global ydbw
		a = np.zeros(nx)
		b = np.ones(nx)
		c = np.zeros(nx)
		d = np.zeros(nx)
		d[0] = ydbw[0]
		d[nx-1] = ydbw[nx-1]
		sig = dt/(4*dx)
		for i in xrange(1,nx-1):
			a[i] = -sig*ydbw[i-1]
			c[i] = sig*ydbw[i+1]
			d[i] = ydbw[i]		# second 2 terms are equal and opposite, so cancel (???)
		ydbw = TDMAsolver(a,b,c,d)

		# implicit Beam-Warming with explicit damping - solve tri-diagonal matrix
		# set up arrays
		global ydbw2
		epsilon = .1		# damping factor
		a = np.zeros(nx)
		b = np.ones(nx)
		c = np.zeros(nx)
		d = np.zeros(nx)
		d[0] = ydbw2[0]
		d[nx-1] = ydbw2[nx-1]
		sig = dt/(4*dx)
		for i in xrange(1,nx-1):
			a[i] = -sig*ydbw2[i-1]
			c[i] = sig*ydbw2[i+1]
			d[i] = ydbw2[i]		# second 2 terms are equal and opposite in Inviscid Burgers eqn, so cancel (???)

			# damping - not clear what to do at edges
			if i >= 2 and i <= nx-3:
				d[i] -= epsilon*(ydbw2[i+2]-4*ydbw2[i+1]+6*ydbw2[i]-4*ydbw2[i-1]+ydbw2[i-2])

		ydbw2 = TDMAsolver(a,b,c,d)

		# implicit Beam-Warming with implicit damping - solve tri-diagonal matrix
		# set up arrays
		global ydbw3
		epsilon = .1		# damping factor
		a = np.zeros(nx)
		b = np.ones(nx)
		c = np.zeros(nx)
		d = np.zeros(nx)
		d[0] = ydbw3[0]
		d[nx-1] = ydbw3[nx-1]
		sig = dt/(4*dx)
		for i in xrange(1,nx-1):
			a[i] = -sig*ydbw3[i-1] - epsilon*ydbw3[i-1]
			b[i] = 1.0 + epsilon*2*ydbw3[i]
			c[i] = sig*ydbw3[i+1] - epsilon*ydbw3[i+1]
			d[i] = ydbw3[i]		# second 2 terms are equal and opposite, so cancel (???)
		ydbw3 = TDMAsolver(a,b,c,d)


		yield


def run (data):
	#print "run()"
#	line_bd.set_ydata(ydbd)
#	line_bdc.set_ydata(ydbdc)
#	#line_fd.set_ydata(ydfd)
#	line_lf.set_ydata(ydlf)
#	line_lfc.set_ydata(ydlfc)
#	line_lw.set_ydata(ydlw)
#	line_mc.set_ydata(ydmc)
#	line_bw.set_ydata(ydbw)
	line_bw2.set_ydata(ydbw2)
	line_bw3.set_ydata(ydbw3)


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)
plt.show()




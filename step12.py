# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 12: Navier-Stokes Equation - Channel Flow
	du/dt + u*du/dx + v*du/dy = -(1/rho)*dp/dx + nu*(d^2u/dx^2 + d^2u/dy^2) + F
	dv/dt + u*dv/dx + v*dv/dy = -(1/rho)*dp/dy + nu*(d^2v/dx^2 + d^2v/dy^2)
	d^2p/dx^2 + d^2p/dy^2 = -rho*((du/dx)^2+2*(du/dy*dv/dx)+(dv/dy)^2) + rho*d/dt(du/dx+dv/dy)

(u,v) are 2D velocity of the fluid
p is the pressure
rho is the density - constant
nu is the viscosity - constant

Add force F = 1 to u equation

The first 2 eqs are N-S, the 3rd is a Poisson eq on the pressure p to enforce div(velocity) = 0 (incompressible flow)
Proceed by alternately solving Poisson eq for p and N-S eqs for u and v
Note rhs of Poisson eq (b) now depends on u and v - not a constant

Forward diff in time and psuedo-time (for Poisson eq)
Centered diff in space for diffusion term in NS and pressure eq
Upwind (Backwards diff) for convection terms in NS

Note: the final term of pressure eq contains d/dt of (du/dx+dv/dy) - it looks like this is discretized
 numerically just by dividing by dt.  Should it be a forward difference?
Yes, but (du/dx+dv/dy) is div(velocity) which we are setting to 0 in the new time point, so a forward
 difference reduces to 1/dt of the spatial derivatives at the old time point (negated)

Domain: [0,2] x [0,2]
Range:
Initial Conditions: u,v,p = 0 everywhere
Boundary Conditions: u,v,p periodic at x = 0,2; u,v = 0 @ y = 0,2; dp/dy = 0 @ y = 0,2


Plot U and V as quiver:
fig = pyplot.figure(figsize = (11,7), dpi=100)
plt.quiver(X,Y,U1,V1)

"""

import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import os
import shutil



# graphing vars
fig = plt.figure(figsize=(18.0,6.0), subplotpars=matplotlib.figure.SubplotParams(left=.05,right=.95))
p1 = fig.add_subplot(131);
p1.set_title("U")
p2 = fig.add_subplot(132);
p2.set_title("V")
p3 = fig.add_subplot(133);
p3.set_title("P")



# frame dir for movie
frmDir = "tmp1"

# physical constants
rho = 1.0		# fluid density
nu = .04 #0.0025 #0.01		# fluid viscosity - scales as 1/dx
Fu = 1	# force term on U

# simulation constants
nx = 51 #201		# num grid X points
ny = 51 #201		# num grid Y points
nt = 300 #201		# number of time steps in sim
nit = 250 #100	# pseudo-time iteration steps per time step (Poisson eq) - MUST BE EVEN - this should probably scale with dx
dt = .004 #.002 #0.01   # time step - scales as 1/dx
dx = 2.0 / (nx-1.0)
dy = 2.0 / (ny-1.0)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
U1 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
U2 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
V1 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
V2 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
P1 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
P2 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
B  = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]		# rhs of pressure eq
ping1to2 = True

# default color map is rainbow
im1 = p1.imshow(U1, vmin=-.1, vmax=1.0)
im2 = p2.imshow(V1, vmin=-.1, vmax=1.0)		#vmin=0.0, vmax=.1)
im3 = p3.imshow(P1, vmin=-1.0, vmax=1.0)

# global var - current time step
ct = 0



def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dy = " + str(dy) + " dt = " + str(dt)
	ct = 0





def data_gen ():
	global ct
	while ct < nt:
		#print "data_gen() ct =", ct
		ct += 1

		uo, un = [], []
		vo, vn = [], []
		if ping1to2:
			uo = U1
			un = U2
			vo = V1
			vn = V2
		else:
			uo = U2
			un = U1
			vo = V2
			vn = V1


		# set up B, the rhs for the pressure eq
		for xi in xrange(nx):
			# enforce periodic BC by wrapping indices
			xim1 = xi - 1
			if (xim1 < 0): xim1 = nx-1
			xip1 = xi + 1
			if (xip1 >= nx): xip1 = 0

			for yi in xrange(1,ny-1):
				dudx = (uo[yi][xip1]-uo[yi][xim1]) / (2.0*dx)
				dudy = (uo[yi+1][xi]-uo[yi-1][xi]) / (2.0*dy)
				dvdx = (vo[yi][xip1]-vo[yi][xim1]) / (2.0*dx)
				dvdy = (vo[yi+1][xi]-vo[yi-1][xi]) / (2.0*dy)
				B[yi][xi] = -rho*(dudx*dudx + 2.0*(dudy*dvdx) + dvdy*dvdy) + (rho/dt)*(dudx+dvdy)

		# update pressure by iterating Poisson eq
		# NOTE: number of iterations must be even to get us back to starting point
		po, pn = [], []
		po = P1
		pn = P2
		for it in xrange(nit):
			for xi in xrange(nx):
				# enforce periodic BC by wrapping indices
				xim1 = xi - 1
				if (xim1 < 0): xim1 = nx-1
				xip1 = xi + 1
				if (xip1 >= nx): xip1 = 0

				for yi in xrange(1,ny-1):
					pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xip1]+po[yi][xim1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]) - dx*dx*dy*dy*B[yi][xi])

			# enforce BC: dp/dy = 0 @ y = 0,2
			#for xi in xrange(1,nx-1):
			for xi in xrange(nx):
				pn[0][xi] = pn[1][xi]
				pn[ny-1][xi] = pn[ny-2][xi]

			# swap old and new buffers
			ptmp = pn
			pn = po
			po = ptmp

		# set pn to newest buffer
		pn = P1

		# do one step in N-S for u and v using new pressure
		for xi in xrange(nx):
			# enforce periodic BC by wrapping indices
			xim1 = xi - 1
			if (xim1 < 0): xim1 = nx-1
			xip1 = xi + 1
			if (xip1 >= nx): xip1 = 0

			for yi in xrange(1,ny-1):
				un[yi][xi] = Fu*dt + uo[yi][xi] - \
								(dt/dx)*uo[yi][xi]*(uo[yi][xi]-uo[yi][xim1]) - (dt/dy)*vo[yi][xi]*(uo[yi][xi]-uo[yi-1][xi]) - \
								(dt/(2.0*rho*dx))*(pn[yi][xip1]-pn[yi][xim1]) + \
								nu*dt*((uo[yi][xip1]-2*uo[yi][xi]+uo[yi][xim1])/(dx*dx) + (uo[yi+1][xi]-2*uo[yi][xi]+uo[yi-1][xi])/(dy*dy))
				vn[yi][xi] = vo[yi][xi] - \
								(dt/dx)*uo[yi][xi]*(vo[yi][xi]-vo[yi][xim1]) - (dt/dy)*vo[yi][xi]*(vo[yi][xi]-vo[yi-1][xi]) - \
								(dt/(2.0*rho*dy))*(pn[yi+1][xi]-pn[yi-1][xi]) + \
								nu*dt*((vo[yi][xip1]-2*vo[yi][xi]+vo[yi][xim1])/(dx*dx) + (vo[yi+1][xi]-2*vo[yi][xi]+vo[yi-1][xi])/(dy*dy))

		yield


def run (data):
	global ping1to2

	#print "run()"
	if ping1to2:
		im1.set_array(U2)
		im2.set_array(V2)
	else:
		im1.set_array(U1)
		im2.set_array(V1)

	# P1 is always the newest pressure grid
	im3.set_array(P1)
	#p1.set_title("P - %d" % ct)

	# save the current figure out as a frame for our movie
	# can build movie on command line with:
	#	ffmpeg -i tmp3/frm%04d.png -r 15 -vcodec mpeg4 -y step10.mp4
	fig.savefig("tmp1/frm%04d.png" % ct, dpi='figure')
	print "Frame " + str(ct)

	ping1to2 = not ping1to2


# delete frmDir if it exists and create it again
if os.path.exists(frmDir):
    shutil.rmtree(frmDir)
os.makedirs(frmDir)

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=1, repeat=False, init_func=init)





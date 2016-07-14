# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 11: Navier-Stokes Equation - Cavity Flow
	du/dt + u*du/dx + v*du/dy = -(1/rho)*dp/dx + nu*(d^2u/dx^2 + d^2u/dy^2)
	dv/dt + u*dv/dx + v*dv/dy = -(1/rho)*dp/dy + nu*(d^2v/dx^2 + d^2v/dy^2)
	d^2p/dx^2 + d^2p/dy^2 = -rho*((du/dx)^2+2*(du/dy*dv/dx)+(dv/dy)^2) + rho*d/dt(du/dx+dv/dy)

(u,v) are 2D velocity of the fluid
p is the pressure
rho is the density - constant
nu is the viscosity - constant

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
Boundary Conditions: u = 1 @ y = 2, u,v = 0 @ x = 0,2 & y = 0, p = 0 @ y = 2, dp/dy = 0 @ y = 0, dp/dx = 0 @ x = 0,2

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
nu = 0.1		# fluid viscosity

# simulation constants
nx = 20		# num grid X points
ny = 20		# num grid Y points
nt = 300		# number of time steps in sim
nit = 100	# pseudo-time iteration steps per time step (Poisson eq) - MUST BE EVEN
dt = 0.01
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
im1 = p1.imshow(U1, vmin=0.0, vmax=1.0)
im2 = p2.imshow(V1, vmin=0.0, vmax=1.0)
im3 = p3.imshow(P1, vmin=0.0, vmax=1.0)

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
		#print "B"
		for yi in xrange(1,ny-1):
			for xi in xrange(1,nx-1):
				#dudx = (uo[yi][xi+1]-uo[yi][xi-1]) / (2.0*dx)
				#dudy = (uo[yi+1][xi]-uo[yi-1][xi]) / (2.0*dy)
				#dvdx = (vo[yi][xi+1]-vo[yi][xi-1]) / (2.0*dx)
				#dvdy = (vo[yi+1][xi]-vo[yi-1][xi]) / (2.0*dy)
				#B[yi][xi] = -rho*(dudx*dudx + 2.0*(dudy*dvdx) + dvdy*dvdy) + (rho/dt)*(dudx+dvdy)
				B[yi][xi] = rho*( ((uo[yi][xi+1]-uo[yi][xi-1])/2.0/dx + (vo[yi+1][xi]-vo[yi-1][xi])/2.0/dx ) / dt \
								+ pow((uo[yi][xi+1]-uo[yi][xi-1])/2.0/dx, 2.0) \
								+ 2.0*((uo[yi+1][xi]-uo[yi-1][xi])/2.0/dy * (vo[yi][xi+1]-vo[yi][xi-1])/2.0/dx) \
								+ pow((vo[yi+1][xi]-vo[yi-1][xi])/2.0/dy, 2.0) )

		# update pressure by iterating Poisson eq
		# NOTE: number of iterations must be even to get us back to starting point
		#print "P"
		po, pn = [], []
		po = P1
		pn = P2
		for it in xrange(nit):
			for yi in xrange(1,ny-1):
				for xi in xrange(1,nx-1):
					#pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xi+1]+po[yi][xi-1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]) - dx*dx*dy*dy*B[yi][xi])
					pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xi+1]+po[yi][xi-1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]) - dx*dx*dy*dy*B[yi][xi])

			# enforce BCs: dp/dy = 0 @ y = 0, dp/dx = 0 @ x = 0,2
			for xi in xrange(nx):
				pn[0][xi] = pn[1][xi]
				pn[ny-1][xi] = pn[ny-2][xi]
			for yi in xrange(ny):
				pn[yi][0] = pn[yi][1]
				pn[yi][nx-1] = pn[yi][nx-2]

			# swap old and new buffers
			ptmp = pn
			pn = po
			po = ptmp

		# set pn to newest buffer
		pn = P1

		# do one step in N-S for u and v using new pressure
		#print "UV"
		for yi in xrange(1,ny-1):
			for xi in xrange(1,nx-1):
#				un[yi][xi] = uo[yi][xi] - \
#								(dt/dx)*uo[yi][xi]*(uo[yi][xi]-uo[yi][xi-1]) - (dt/dy)*vo[yi][xi]*(uo[yi][xi]-uo[yi-1][xi]) - \
#								(dt/(2.0*rho*dx))*(pn[yi][xi+1]-pn[yi][xi-1]) + \
#								nu*dt*((uo[yi][xi+1]-2*uo[yi][xi]+uo[yi][xi-1])/(dx*dx) + (uo[yi+1][xi]-2*uo[yi][xi]+uo[yi-1][xi])/(dy*dy))
#				vn[yi][xi] = vo[yi][xi] - \
#								(dt/dx)*uo[yi][xi]*(vo[yi][xi]-vo[yi][xi-1]) - (dt/dy)*vo[yi][xi]*(vo[yi][xi]-vo[yi-1][xi]) - \
#								(dt/(2.0*rho*dy))*(pn[yi+1][xi]-pn[yi-1][xi]) + \
#								nu*dt*((vo[yi][xi+1]-2*vo[yi][xi]+vo[yi][xi-1])/(dx*dx) + (vo[yi+1][xi]-2*vo[yi][xi]+vo[yi-1][xi])/(dy*dy))
				#print "xi = " + str(xi) + " yi = " + str(yi)
				un[yi][xi] = uo[yi][xi] \
							- uo[yi][xi]*dt/dx*(uo[yi][xi]-uo[yi][xi-1])  \
							- vo[yi][xi]*dt/dy*(uo[yi][xi]-uo[yi-1][xi]) \
							- 1.0/rho*(pn[yi][xi+1]-pn[yi][xi-1])*dt/2.0/dx \
							+ nu*dt/dx/dx*(uo[yi][xi+1]-2.0*uo[yi][xi]+uo[yi][xi-1]) \
							+ nu*dt/dy/dy*(uo[yi+1][xi]-2.0*uo[yi][xi]+uo[yi-1][xi])
				vn[yi][xi] = vo[yi][xi] \
							- uo[yi][xi]*dt/dx*(vo[yi][xi]-vo[yi][xi-1]) \
							- vo[yi][xi]*dt/dy*(vo[yi][xi]-vo[yi-1][xi]) \
							- 1.0/rho*(pn[yi+1][xi]-pn[yi-1][xi])*dt/2.0/dy \
							+ nu*dt/dx/dx*(vo[yi][xi+1]-2.0*vo[yi][xi]+vo[yi][xi-1]) \
							+ nu*dt/dy/dy*(vo[yi+1][xi]-2.0*vo[yi][xi]+vo[yi-1][xi])

		# enforce BCs: u = 1 @ y = 2, u,v = 0 @ x = 0,2 & y = 0
		# NO - set at start and time-iteration doesn't touch edges, so shouldn't have to
		# BC: u = 1 @ y = 2
		#print "BC"
		for xi in xrange(nx):
			un[0][xi] = 0.0
			un[ny-1][xi] = 1.0
			vn[0][xi] = 0.0
			vn[ny-1][xi] = 0.0
		for yi in xrange(ny):
			un[yi][0] = 0.0
			un[yi][nx-1] = 0.0
			vn[yi][0] = 0.0
			vn[yi][nx-1] = 0.0

		yield


def run (data):
	global ping1to2

	#print "run()"
	if ping1to2:
		im1.set_array(U2)
		im2.set_array(V2)
		im3.set_array(P2)
	else:
		im1.set_array(U1)
		im2.set_array(V1)
		im3.set_array(P1)
	#p1.set_title("P - %d" % ct)

	# save the current figure out as a frame for our movie
	# can build movie on command line with:
	#	ffmpeg -i tmp3/frm%04d.png -r 15 -vcodec mpeg4 -y step10.mp4
	fig.savefig(frmDir + "/frm%04d.png" % ct, dpi='figure')
	print "Frame " + str(ct)

	ping1to2 = not ping1to2


# delete frmDir if it exists and create it again
if os.path.exists(frmDir):
    shutil.rmtree(frmDir)
os.makedirs(frmDir)

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=1, repeat=False, init_func=init)





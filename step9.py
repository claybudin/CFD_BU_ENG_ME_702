# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 9: 2D Laplace Equation
	d^2p/dx^2 + d^2p/dy^2 = 0

Solving for pressure, p(x,y) -> R (scalar)

No time dependency (steady-state solution) but we will use an interative pseudo-time approach
to solving

Forward diff in psuedo-time
Centered diff in space

Domain: [0,2] x [0,1]
Range: [0,1]
Initial Conditions: p = 0 everywhere (except at boundary condition below?)
Boundary Conditions: p = 0 @ x = 0, p = y @ x = 2, dp/dy = 0 @ y = 0,1
  these BCs seem to conflict at x = 2 since dp/dy there = 1

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig = plt.figure(figsize=[16.0,6.0])
#print str(fig.get_figwidth()) + " x " + str(fig.get_figheight())
p1 = fig.add_subplot(121);
p2 = fig.add_subplot(122);

p1.set_title("P")
p2.set_title("P - Analytic Solution")



# simulation constants
nx = 201
ny = 101
nt = 100000	# pseudo-time iteration steps - wow, need a lot of these
dx = 2.0 / (nx-1.0)
dy = 1.0 / (ny-1.0)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
P1 = [ [ 0.0 for xi in range(nx) ] for yi in range(ny) ]
P2 = [ [ 0.0 for xi in range(nx) ] for yi in range(ny) ]
PA = [ [ 0.0 for xi in range(nx) ] for yi in range(ny) ]
ping1to2 = True

# default color map is rainbow
#im1 = p1.imshow(P1, vmin=0.0, vmax=2.0)
im1 = p1.imshow(P1, vmin=0.0, vmax=1.0)
im2 = p2.imshow(PA, vmin=0.0, vmax=1.0)




# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dy = " + str(dy)
	ct = 0

	pmin = 1000000000.0
	pmax = -1000000000.0
	for yi in range(ny):
		for xi in range(nx):
			x = xi*dx
			y = yi*dy

			# initialize arrays
			p = 0.0
			# boundary condition p = y @ x = 2
			if x == 2.0: p = y

			# analytic solution
			sum = 0.0
			for n in range(10):
				nn = 2.0*n+1.0
				sum += (np.sinh(nn*np.pi*x)*np.cos(nn*np.pi*y)) / ((nn*np.pi)*(nn*np.pi)*np.sinh(2*np.pi*nn))
			pa = x/4.0 - 4.0*sum

			P1[yi][xi] = p
			P2[yi][xi] = p
			PA[yi][xi] = pa

			if pa < pmin: pmin = pa
			if pa > pmax: pmax = pa

			#if yi == 10: U1[yi][xi] = x    # visualize color range

	im2.set_array(PA)
	print "pmin = " + str(pmin) + " pmax = " + str(pmax)


def data_gen ():
	global ct
	while ct < nt:
		#print "data_gen() ct =", ct
		ct += 1

		po, pn = [], []
		if ping1to2:
			po = P1
			pn = P2
		else:
			po = P2
			pn = P1

		# iterate
		# wouldn't be hard to run this on multiple CPUs - all would need access to po and write to different places in pn
		for yi in range(1,ny-1):
			for xi in range(1,nx-1):
				pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xi+1]+po[yi][xi-1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]))

		# enforce boundary conditions dp/dy = 0 @ y = 0,1
		# not sure why we aren't changing x values at indices 0 and nx-1 - guess to enfore BCs on x
		for xi in range(1,nx-1):
			pn[0][xi] = pn[1][xi]
			pn[ny-1][xi] = pn[ny-2][xi]

		yield


def run (data):
	global ping1to2

	#print "run()"
	#print "ping1to2 = " + str(ping1to2)
	if ct % 100 == 0:
		if ping1to2:
			im1.set_array(P2)
		else:
			im1.set_array(P1)
		p1.set_title("P - %d" % ct)

	ping1to2 = not ping1to2


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=1, repeat=False, init_func=init)

#plt.show()




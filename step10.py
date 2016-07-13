# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 10: 2D Poisson Equation
	d^2p/dx^2 + d^2p/dy^2 = b

Solving for pressure, p(x,y) -> R (scalar)
b is a fixed scalar function defined over domain

No time dependency (steady-state solution) but we will use an interative pseudo-time approach
to solving

Forward diff in psuedo-time
Centered diff in space

Domain: [0,2] x [0,1]
Range: [-100,100]
Initial Conditions: p = 0
Boundary Conditions: p = 0 @ x = 0,2, y = 0,1

Source (function b):
	b = 100 @ (nx/4, ny/4)
	b = -100 @ (3nx/4, 3ny/4)
	b = 0 everywhere else
made B spikes +/- 10000 because otherwise they were having almost no effect
despite large values of B, they seem to have a tiny effect on the graph, due probably to the fact
 that B is scaled by dx*dx*dy*dy

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig = plt.figure(figsize=[8.0,6.0])
p1 = fig.add_subplot(111);
p1.set_title("P")



# simulation constants
nx = 201
ny = 101
nt = 30000	# pseudo-time iteration steps - wow, need a lot of these
dx = 2.0 / (nx-1.0)
dy = 1.0 / (ny-1.0)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
P1 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
P2 = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
B = [ [ 0.0 for xi in xrange(nx) ] for yi in xrange(ny) ]
ping1to2 = True

# default color map is rainbow
im1 = p1.imshow(P1, vmin=-1.0, vmax=1.0)	# doesn't seem to get anywhere near +/-100 of B spikes


# global var - current time step
ct = 0



def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dy = " + str(dy)
	ct = 0

	B[ny/4][nx/4] = 10000.0
	B[3*ny/4][3*nx/4] = -10000.0


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
		for yi in xrange(1,ny-1):
			for xi in xrange(1,nx-1):
				pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xi+1]+po[yi][xi-1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]) - dx*dx*dy*dy*B[yi][xi])

		yield


def run (data):
	global ping1to2

	#print "run()"
	if ct == 1 or ct % 100 == 0:
		if ping1to2:
			im1.set_array(P2)
		else:
			im1.set_array(P1)
		p1.set_title("P - %d" % ct)

		# save the current figure out as a frame for our movie
		# can build movie on command line with:
		#	ffmpeg -i tmp3/frm%04d.png -r 15 -vcodec mpeg4 -y step10.mp4
		fig.savefig("tmp3/frm%04d.png" % (ct/100), dpi='figure')


	ping1to2 = not ping1to2


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=1, repeat=False, init_func=init)





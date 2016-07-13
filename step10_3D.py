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
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

import time



# graphing vars
fig = plt.figure(figsize=[8.0,6.0])
p1 = fig.add_subplot(111, projection='3d');
p1.set_title("P - 0")




# simulation constants
nx = 201
ny = 101
nt = 30000	# pseudo-time iteration steps - wow, need a lot of these
dx = 2.0 / (nx-1.0)
dy = 1.0 / (ny-1.0)


# for plotting
xs = np.linspace(0.0, 2.0, nx)
ys = np.linspace(0.0, 1.0, ny)
X, Y = np.meshgrid(xs, ys)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
P1 = X + Y
P2 = X + Y
B = X + Y
ping1to2 = True


# initial conditions
for yi in xrange(ny):
	for xi in xrange(nx):
		x = xi*dx
		y = yi*dy

		# initialize arrays
		p = 0.0

		# initialize source
		b = 0.0
		if xi == nx/4 and yi == ny/4: b = 10000.0
		if xi == 3*nx/4 and yi == 3*ny/4: b = -10000.0

		P1[yi][xi] = p
		P2[yi][xi] = p
		B[yi][xi] = b


# plot in 3D - based on wire3d_animation_demo.py
wframe = None
tstart = time.time()

frmNum = 0
for ct in xrange(nt):

	oldcol = wframe

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


	if ct % 100 == 0:
		#wframe = p1.plot_wireframe(X, Y, un, rstride=5, cstride=5)	# default strides are 1
		wframe = p1.plot_surface(X, Y, pn, rstride=5, cstride=5, cmap=cm.magma, vmin=-1.0, vmax=1.0, linewidth=0)	# default strides are 10
		#wframe = p1.scatter(X, Y, un, s=1, c='b')	# NOTE: doesn't have strides

		# for some reason, this moves down through surface plot as frames progress
		p1.set_title("P - %d" % ct)

		# Remove old line collection before drawing
		if oldcol is not None:
			p1.collections.remove(oldcol)

		# save the current figure out as a frame for our movie
		# can build movie on command line with:
		#	ffmpeg -i tmp4/frm%04d.png -r 15 -vcodec mpeg4 -y step9_3Dsurf.mp4
		fig.savefig("tmp4/frm%04d.png" % frmNum, dpi='figure')
		frmNum += 1

	plt.pause(.001)

	ping1to2 = not ping1to2


print('Done. FPS: %f' % (nt / (time.time() - tstart)))




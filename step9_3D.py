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
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

import time



# graphing vars
fig = plt.figure(figsize=[16.0,6.0])
#print str(fig.get_figwidth()) + " x " + str(fig.get_figheight())
p1 = fig.add_subplot(121, projection='3d');
p2 = fig.add_subplot(122, projection='3d');

p1.set_title("P - 0")
p2.set_title("P - Analytic")




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
PA = X + Y
ping1to2 = True


# initial conditions
for yi in xrange(ny):
	for xi in xrange(nx):
		x = xi*dx
		y = yi*dy

		# initialize arrays
		p = 0.0
		# boundary condition p = y @ x = 2
		if x == 2.0: p = y

		# analytic solution
		sum = 0.0
		for n in xrange(10):
			nn = 2.0*n+1.0
			sum += (np.sinh(nn*np.pi*x)*np.cos(nn*np.pi*y)) / ((nn*np.pi)*(nn*np.pi)*np.sinh(2*np.pi*nn))
		pa = x/4.0 - 4.0*sum

		P1[yi][xi] = p
		P2[yi][xi] = p
		PA[yi][xi] = pa


#ax.scatter(X, Y, U1, s=1, c='b')	# doesn't have strides
#ax.plot_wireframe(X, Y, U1, rstride=10, cstride=10)

p2.plot_surface(X, Y, PA, rstride=5, cstride=5, cmap=cm.magma, vmin=0.0, vmax=1.0, linewidth=0)


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
			pn[yi][xi] = (1.0/(2.0*(dx*dx+dy*dy))) * (dy*dy*(po[yi][xi+1]+po[yi][xi-1]) + dx*dx*(po[yi+1][xi]+po[yi-1][xi]))

	# enforce boundary conditions dp/dy = 0 @ y = 0,1
	# not sure why we aren't changing x values at indices 0 and nx-1 - guess to enfore BCs on x
	for xi in xrange(1,nx-1):
		pn[0][xi] = pn[1][xi]
		pn[ny-1][xi] = pn[ny-2][xi]

	if ct % 100 == 0:
		#wframe = p1.plot_wireframe(X, Y, un, rstride=5, cstride=5)	# default strides are 1
		wframe = p1.plot_surface(X, Y, pn, rstride=5, cstride=5, cmap=cm.magma, vmin=0.0, vmax=1.0, linewidth=0)	# default strides are 10
		#wframe = p1.scatter(X, Y, un, s=1, c='b')	# NOTE: doesn't have strides

		p1.set_title("P - %d" % ct)

		# Remove old line collection before drawing
		if oldcol is not None:
			p1.collections.remove(oldcol)

		# save the current figure out as a frame for our movie
		# can build movie on command line with:
		#	ffmpeg -i tmp1/frm%04d.png -vframes 301 -r 15 -vcodec mpeg4 -y step9_3Dsurf.mp4
		fig.savefig("tmp1/frm%04d.png" % frmNum, dpi='figure')
		frmNum += 1

	plt.pause(.001)

	ping1to2 = not ping1to2


print('Done. FPS: %f' % (nt / (time.time() - tstart)))




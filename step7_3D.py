# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 7: 2D diffusion
	du/dt - visc*(d^2u/dx^2 + d^2u/dy^2) = 0

Forward diff in time
Backward diff in space

Domain: [0,2]
Range: [1,2]
Initial Conditions: square wave, half-sine, full inverted cosine
Boundary Conditions: u = 1 @ x = 0,2, y = 0,2

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as anim

import time
import sys



# graphing vars
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# movie generations
#matplotlib.use("Agg")
#FFMpegWriter = anim.writers['ffmpeg']
matplotlib.verbose.level = 'debug'
matplotlib.verbose.fileo = sys.stdout
metadata = dict(title='"Movie Test"', artist='Matplotlib', comment='"Movie support!"')
#writer = FFMpegWriter(fps=15, metadata=metadata)
writer = anim.FFMpegFileWriter(fps=15, metadata=metadata)



# simulation constants
nx = 201
ny = 201
nt = 150
visc = .0025
dt = 0.01
dx = 2.0 / (nx-1.0)
dy = 2.0 / (ny-1.0)


# for plotting
xs = np.linspace(0.0, 2.0, nx)
ys = np.linspace(0.0, 2.0, ny)
X, Y = np.meshgrid(xs, ys)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
U1 = X + Y
U2 = X + Y
ping1to2 = True


# initial conditions
for yi in range(ny):
	for xi in range(nx):
		x = xi*dx
		y = yi*dy

		# initialize range [.5,1] with various test fns
		v = 1.0
		if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 2.0		# square wave
		#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 1.0 + np.sin((x-.5)*np.pi*2.0)*np.sin((y-.5)*np.pi*2.0)	# half-sine
#		if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0:				# full inverted cosine
#			#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#			v = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0)) * .5*(1.0 - np.cos((y-.5)*np.pi*4.0))
#			#print "y = " + str(y)
		#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	 * .5*np.sin((y-.5)*np.pi*4.0)      # full-sine

		U1[yi][xi] = v
		U2[yi][xi] = v

		#if yi == 2: U1[yi][xi] = x    # visualize color range


#ax.scatter(X, Y, U1, s=1, c='b')	# doesn't have strides
#ax.plot_wireframe(X, Y, U1, rstride=10, cstride=10)


# plot in 3D - based on wire3d_animation_demo.py
wframe = None
tstart = time.time()

# something's broken in the writer - this will run through the animation sequance and save the
#  files but will crash when it goes to build the movie
# fortunately it prints out the command line it is trying to use so it can be run manually and
#  then the temp image frame files can be deleted
# might need to set something in the matplotlibrc or something
# it seems to be complaining about bad file handles
#
with writer.saving(fig, "step7_3D.mp4", fig.dpi):

	for ct in range(nt):

		oldcol = wframe

		o, n = [[]], [[]]
		if ping1to2:
			o = U1
			n = U2
		else:
			o = U2
			n = U1

		for yi in range(1,ny-1):
			for xi in range(1,nx-1):
				n[yi][xi] = o[yi][xi] + (visc*dt/(dx*dx))*(o[yi][xi+1]-2*o[yi][xi]+o[yi][xi-1]) + \
									    (visc*dt/(dy*dy))*(o[yi+1][xi]-2*o[yi][xi]+o[yi-1][xi])

		wframe = ax.plot_wireframe(X, Y, n, rstride=5, cstride=5)	# default strides are 1
		#wframe = ax.plot_surface(X, Y, n, rstride=5, cstride=5)	# default strides are 10
		#wframe = ax.scatter(X, Y, n, s=1, c='b')	# NOTE: doesn't have strides

		# Remove old line collection before drawing
		if oldcol is not None:
			ax.collections.remove(oldcol)

		plt.pause(.001)

		ping1to2 = not ping1to2

		writer.grab_frame()

print('Done. FPS: %f' % (nt / (time.time() - tstart)))




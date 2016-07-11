# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 5: 2D linear convection
	du/dt + c*du/dx + c*du/dy = 0
linear because c is constant

Forward diff in time
Backward diff in space

Domain: [0,2]
Range: [1,2]
Initial Conditions: square wave, half-sine, full inverted cosine
Boundary Conditions: u = 1 @ x = 0,2, y = 0,2

stability depends on c*dt/dx and c*dt/dy
looks like sum should be <= 1 to remain stable

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

import time



# graphing vars
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')



# simulation constants
nx = 201
ny = 201
nt = 150
c = 0.5 #0.5 #1.0
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
			n[yi][xi] = o[yi][xi] - (c*dt/dx)*(o[yi][xi]-o[yi][xi-1]) - (c*dt/dy)*(o[yi][xi]-o[yi-1][xi])
			#ydata[i] = yprev[i] - (c*dt/(2.0*dx))*(yprev[i+1]-yprev[i-1])		# central diff - not stable for any c (?)
			#ydata[i] = yprev[i] + (c*dt/dx)*(yprev[i+1]-yprev[i])				# change sign, use forward diff - wave moves to left

	wframe = ax.plot_wireframe(X, Y, n, rstride=5, cstride=5)	# default strides are 1
	#wframe = ax.plot_surface(X, Y, n, rstride=5, cstride=5)	# default strides are 10
	#wframe = ax.scatter(X, Y, n, s=1, c='b')	# NOTE: doesn't have strides

	# Remove old line collection before drawing
	if oldcol is not None:
		ax.collections.remove(oldcol)

	plt.pause(.001)

	ping1to2 = not ping1to2

print('Done. FPS: %f' % (nt / (time.time() - tstart)))




# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 6: 2D non-linear convection
	du/dt + u*du/dx + v*du/dy = 0
	dv/dt + u*dv/dx + v*dv/dy = 0

Now 2 output variables (or one 2D solution function): u and v

Forward diff in time
Backward diff in space

Domain: [0,2]
Range: [1,2]
Initial Conditions: square wave, half-sine, full inverted cosine
Boundary Conditions: u,v = 1 @ x = 0,2, y = 0,2


"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

import time



# graphing vars
fig = plt.figure(figsize=[16.0,6.0])
#print str(fig.get_figwidth()) + " x " + str(fig.get_figheight())
p1 = fig.add_subplot(121, projection='3d');
p2 = fig.add_subplot(122, projection='3d');

p1.set_title("U")
p2.set_title("V")




# simulation constants
nx = 201
ny = 201
nt = 150
c = 0.5 #0.5 #1.0
dt = 0.0025
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
V1 = X + Y
V2 = X + Y
ping1to2 = True


# initial conditions
for yi in range(ny):
	for xi in range(nx):
		x = xi*dx
		y = yi*dy

		# initialize range [.5,1] with various test fns
		u = 1.0
		if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 2.0		# square wave
		#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 1.0 + np.sin((x-.5)*np.pi*2.0)*np.sin((y-.5)*np.pi*2.0)	# half-sine
#		if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0:				# full inverted cosine
#			#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#			u = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0)) * .5*(1.0 - np.cos((y-.5)*np.pi*4.0))
#			#print "y = " + str(y)
		#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	 * .5*np.sin((y-.5)*np.pi*4.0)      # full-sine

		#v = u
		v = 1.0
		if x >= .75 and x <= 1.25 and y >= .75 and y <= 1.25: v = 2.0		# square wave

		U1[yi][xi] = u
		U2[yi][xi] = u
		V1[yi][xi] = v
		V2[yi][xi] = v

		#if yi == 2: U1[yi][xi] = x    # visualize color range


#ax.scatter(X, Y, U1, s=1, c='b')	# doesn't have strides
#ax.plot_wireframe(X, Y, U1, rstride=10, cstride=10)


# plot in 3D - based on wire3d_animation_demo.py
wframe1 = None
wframe2 = None
tstart = time.time()

for ct in range(nt):

	oldcol1 = wframe1
	oldcol2 = wframe2

	uo, un = [[]], [[]]
	vo, vn = [[]], [[]]
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

	for yi in range(1,ny-1):
		for xi in range(1,nx-1):
			un[yi][xi] = uo[yi][xi] - (uo[yi][xi]*dt/dx)*(uo[yi][xi]-uo[yi][xi-1]) - (vo[yi][xi]*dt/dy)*(uo[yi][xi]-uo[yi-1][xi])
			vn[yi][xi] = vo[yi][xi] - (uo[yi][xi]*dt/dx)*(vo[yi][xi]-vo[yi][xi-1]) - (vo[yi][xi]*dt/dy)*(vo[yi][xi]-vo[yi-1][xi])
			#ydata[i] = yprev[i] - (c*dt/(2.0*dx))*(yprev[i+1]-yprev[i-1])		# central diff - not stable for any c (?)
			#ydata[i] = yprev[i] + (c*dt/dx)*(yprev[i+1]-yprev[i])				# change sign, use forward diff - wave moves to left

	wframe1 = p1.plot_wireframe(X, Y, un, rstride=5, cstride=5)	# default strides are 1
	#wframe1 = p1.plot_surface(X, Y, un, rstride=5, cstride=5)	# default strides are 10
	#wframe1 = p1.scatter(X, Y, un, s=1, c='b')	# NOTE: doesn't have strides

	wframe2 = p2.plot_wireframe(X, Y, vn, rstride=5, cstride=5)	# default strides are 1
	#wframe2 = p2.plot_surface(X, Y, vn, rstride=5, cstride=5)	# default strides are 10
	#wframe2 = p2.scatter(X, Y, vn, s=1, c='b')	# NOTE: doesn't have strides

	# Remove old line collection before drawing
	if oldcol1 is not None:
		p1.collections.remove(oldcol1)
	if oldcol2 is not None:
		p2.collections.remove(oldcol2)

	plt.pause(.001)

	ping1to2 = not ping1to2

print('Done. FPS: %f' % (nt / (time.time() - tstart)))




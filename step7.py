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
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# graphing vars
fig = plt.figure()



# simulation constants
nx = 201
ny = 201
nt = 150
visc = .0025
dt = 0.01
dx = 2.0 / (nx-1.0)
dy = 2.0 / (ny-1.0)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
U1 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
U2 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
ping1to2 = True

# default color map is rainbow
im = plt.imshow(U1, vmin=0.0, vmax=2.0)




# global var - current time step
ct = 0

def init ():
	global ct
	if ct == 0: print "init() dx = " + str(dx) + " dt = " + str(dt)
	ct = 0

	for yi in range(ny):
		for xi in range(nx):
			x = xi*dx
			y = yi*dy

			# initialize range [.5,1] with various test fns
			v = 1.0
			if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 2.0		# square wave
			#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 1.0 + np.sin((x-.5)*np.pi*2.0)*np.sin((y-.5)*np.pi*2.0)	# half-sine
#			if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0:				# full inverted cosine
#				#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#				v = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0)) * .5*(1.0 - np.cos((y-.5)*np.pi*4.0))
#				#print "y = " + str(y)
			#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: v = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	 * .5*np.sin((y-.5)*np.pi*4.0)      # full-sine

			U1[yi][xi] = v

			#if yi == 10: U1[yi][xi] = x    # visualize color range


def data_gen ():
	global ct
	while ct <= nt:
		#print "data_gen() ct =", ct
		ct += 1

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
		yield


def run (data):
	global ping1to2

	#print "run()"
	#print "ping1to2 = " + str(ping1to2)
	if ping1to2:
		im.set_array(U2)
	else:
		im.set_array(U1)

	ping1to2 = not ping1to2


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

plt.show()




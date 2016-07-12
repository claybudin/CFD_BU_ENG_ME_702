# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 8: 2D Burgers's Equation
	du/dt + u*du/dx + v*du/dy - visc*(d^2u/dx^2 + d^2u/dy^2) = 0
	dv/dt + u*dv/dx + v*dv/dy - visc*(d^2v/dx^2 + d^2v/dy^2) = 0

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
import matplotlib.animation as animation

# graphing vars
fig = plt.figure(figsize=[16.0,6.0])
#print str(fig.get_figwidth()) + " x " + str(fig.get_figheight())
p1 = fig.add_subplot(121);
p2 = fig.add_subplot(122);

p1.set_title("U")
p2.set_title("V")



# simulation constants
nx = 201
ny = 201
nt = 200
visc = 0.0025
dt = .002 #.0015 #0.0025 #0.01
dx = 2.0 / (nx-1.0)
dy = 2.0 / (ny-1.0)

# simlation grids - 2D - ping-pong between them
# access as A[yi][xi]
U1 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
U2 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
V1 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
V2 = [ [ 1.0 for xi in range(nx) ] for yi in range(ny) ]
ping1to2 = True

# default color map is rainbow
im1 = p1.imshow(U1, vmin=0.0, vmax=2.0)
im2 = p2.imshow(V1, vmin=0.0, vmax=2.0)




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
			u = 1.0
			if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 2.0		# square wave
			#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 1.0 + np.sin((x-.5)*np.pi*2.0)*np.sin((y-.5)*np.pi*2.0)	# half-sine
#			if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0:				# full inverted cosine
#				#print "t = " + str(np.cos((x-.5)*np.pi*4.0))
#				u = 1.0 + .5*(1.0 - np.cos((x-.5)*np.pi*4.0)) * .5*(1.0 - np.cos((y-.5)*np.pi*4.0))
#				#print "y = " + str(y)
			#if x >= .5 and x <= 1.0 and y >= .5 and y <= 1.0: u = 1.0 + .5*np.sin((x-.5)*np.pi*4.0)	 * .5*np.sin((y-.5)*np.pi*4.0)      # full-sine

			#v = u
			v = 1.0
			if x >= .75 and x <= 1.25 and y >= .75 and y <= 1.25: v = 2.0		# square wave - offset from u

			U1[yi][xi] = u
			U2[yi][xi] = u
			V1[yi][xi] = v
			V2[yi][xi] = v


			#if yi == 10: U1[yi][xi] = x    # visualize color range


def data_gen ():
	global ct
	while ct <= nt:
		#print "data_gen() ct =", ct
		ct += 1

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
				un[yi][xi] = uo[yi][xi] - (uo[yi][xi]*dt/dx)*(uo[yi][xi]-uo[yi][xi-1]) - (vo[yi][xi]*dt/dy)*(uo[yi][xi]-uo[yi-1][xi]) + \
								 (visc*dt/(dx*dx))*(uo[yi][xi+1]-2*uo[yi][xi]+uo[yi][xi-1]) + (visc*dt/(dy*dy))*(uo[yi+1][xi]-2*uo[yi][xi]+uo[yi-1][xi])
				vn[yi][xi] = vo[yi][xi] - (uo[yi][xi]*dt/dx)*(vo[yi][xi]-vo[yi][xi-1]) - (vo[yi][xi]*dt/dy)*(vo[yi][xi]-vo[yi-1][xi]) + \
								 (visc*dt/(dx*dx))*(vo[yi][xi+1]-2*vo[yi][xi]+vo[yi][xi-1]) + (visc*dt/(dy*dy))*(vo[yi+1][xi]-2*vo[yi][xi]+vo[yi-1][xi])
		yield


def run (data):
	global ping1to2

	#print "run()"
	#print "ping1to2 = " + str(ping1to2)
	if ping1to2:
		im1.set_array(U2)
		im2.set_array(V2)
	else:
		im1.set_array(U1)
		im2.set_array(V1)

	ping1to2 = not ping1to2


ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=True, init_func=init)

#plt.show()




# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:34:06 2016

@author: clay_budin
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)


#x = [ 0,0,1 ]
#y = [ 0,0,0 ]
#z = [ 0,1,0 ]
##ax.plot(x,y,z, 'bo')
#ax.plot(x,y,z, 'r+')
#

npts = 25 #101
X, Y, Z = [], [], []

for yi in range(npts):
	for xi in range(npts):
		y = yi / (npts-1.0)
		x = xi / (npts-1.0)

		#z = 0
		z = np.cos(x*np.pi*2) * np.sin(y*np.pi*2)

		X.append(x)
		Y.append(y)
		Z.append(z)


#ax.plot(X,Y,Z, 'bo')
#ax.plot(X,Y,Z, 'b^')
#ax.plot(X,Y,Z, 'b-')
#ax.plot(X,Y,Z, 'b.')



#t = 0.0
#def gen ():
#	#print "gen()"
#	global t
#	while True:
#		t += .1
#		i = 0
#		for yi in range(npts):
#			for xi in range(npts):
#				y = yi / (npts-1.0)
#				x = xi / (npts-1.0)
#
#				z = np.cos(x*np.pi*2+t) * np.sin(y*np.pi*2+t)
#				Z[i] = z
#				i += 1
#
#		yield
#
#def run (data):
#	#print "run()"
#	ax.plot(X,Y,Z, 'b.')
#
#
#a = animation.FuncAnimation(fig, run, gen, interval=10, blit=False)
#

#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X,Y,Z)

# slows down over time - something is being accumulated - see wire3d_animation_demo.py
t = 0.0
def run (data):
	#print "run()"
	global t
	t += .1
	i = 0
	for yi in range(npts):
		for xi in range(npts):
			y = yi / (npts-1.0)
			x = xi / (npts-1.0)

			z = np.cos(x*np.pi*2+t) * np.sin(y*np.pi*2+t)
			Z[i] = z
			i += 1

	ax.plot(X,Y,Z, 'b.')


a = animation.FuncAnimation(fig, run, interval=10)

#plt.show()



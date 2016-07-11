#!/usr/bin/env python
"""
An animated image
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def f (x, y):
	#return np.sin(x) + np.cos(y)
	return x+y

# x is horizontal (regular) array, y is vertical 1xn array
x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)

fig = plt.figure()
#im = plt.imshow(f(x, y), cmap=plt.get_cmap('viridis'), animated=True)
#im = plt.imshow(f(x, y))


#im = plt.imshow([i for i in range(256)], shape=[256,]))

X = [ [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	  [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	  [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	  [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	  [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	  [ 255,255,255,255 ],
	  [ 0,0,0,0 ],
	]

#im = plt.imshow(X)

#im = plt.imshow([ [ xi*yi*.01 for xi in range(1,11) ] for yi in range(1,11) ] )
#im = plt.imshow([ [ xi*yi*.0000325 for xi in range(640) ] for yi in range(480) ] )

A = [ [ xi*yi*.00003 for xi in range(320) ] for yi in range(240) ]
A[50][10] = 1.0		# A[y][x] y goes top-to-bottom
im = plt.imshow(A)


def updatefig (*args):
	print "updatefig()"

	global x, y
	x += np.pi / 15.
	y += np.pi / 20.
	#im.set_array(f(x, y))

	#return im,
	yield

def run (data):
	#print "run()"
	im.set_array(f(x, y))


#ani = animation.FuncAnimation(fig, run, updatefig, interval=50, blit=False)

plt.show()






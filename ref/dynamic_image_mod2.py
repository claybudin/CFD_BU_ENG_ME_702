#!/usr/bin/env python
"""
An animated image
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import random


fig = plt.figure()
#im = plt.imshow(f(x, y), cmap=plt.get_cmap('viridis'), animated=True)
#im = plt.imshow(f(x, y))


#im = plt.imshow([i for i in range(256)], shape=[256,]))

# 2D - image representation
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
	#print "updatefig()"

	#A[10][200] = 1.0
	#A[random.randint(0,239)][random.randint(0,319)] = 1.0
	#yield

	for yi in range(240):
		for xi in range(320):
			#print "xi = " + str(xi) + " yi = " + str(yi)
			A[yi][xi] = 1.0
			yield

def run (data):
	#print "run()"
	im.set_array(A)


ani = animation.FuncAnimation(fig, run, updatefig, interval=10, blit=False)

plt.show()






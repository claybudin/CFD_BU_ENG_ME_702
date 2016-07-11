from __future__ import print_function
"""
A very simple 'animation' of a 3D plot
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import time


def generate(X, Y, phi):
	R = 1 - np.sqrt(X**2 + Y**2)
	return np.cos(2 * np.pi * X + phi) * R
	#return 0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = np.linspace(-1, 1, 50)
ys = np.linspace(-1, 1, 50)
#xs = np.linspace(-1, 1, 200)
#ys = np.linspace(-1, 1, 200)

X, Y = np.meshgrid(xs, ys)
Z = generate(X, Y, 0.0)

wframe = None
tstart = time.time()
#for phi in np.linspace(0, 360 / 2 / np.pi, 100):
for phi in np.linspace(0, 10 * 360 / 2 / np.pi, 10 * 100):

	oldcol = wframe

	Z = generate(X, Y, phi)

	#wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
	#wframe = ax.plot_wireframe(X, Y, Z)	# default strides are 1
	#wframe = ax.plot_surface(X, Y, Z, rstride=2, cstride=2)	# default strides are 10
	#wframe = ax.plot(X, Y, Z, 'b.')	# not working
	#wframe = ax.plot3d(X, Y, Z, 'b.')	# not working
	wframe = ax.scatter(X, Y, Z, s=1, c='b')	# doesn't have strides

	# Remove old line collection before drawing
	if oldcol is not None:
		ax.collections.remove(oldcol)

	plt.pause(.001)

print('FPS: %f' % (100 / (time.time() - tstart)))




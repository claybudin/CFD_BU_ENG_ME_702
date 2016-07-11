# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1)
ax.grid()
xdata, ydata = [], []


def data_gen (t=0):
#cnt = 0
#while cnt < 1000:
#		cnt += 1
#		t += 0.1
#		yield t, np.sin(2*np.pi*t) * np.exp(-t/10.)

	while True:
		print "data_gen() t = ", t
		t += .1
		for i in range(100):
			x = i * .1
			y = np.sin(2*np.pi*x+t) * np.exp(-x/10.0)
			xdata[i] = x
			ydata[i] = y
		yield


def init ():
	print "init()"
	ax.set_ylim(-1.1, 1.1)
	ax.set_xlim(0, 10)
	del xdata[:]
	del ydata[:]
	for i in range(100):
		xdata.append(0)
		ydata.append(0)
	line.set_data(xdata, ydata)
	return line,


def run (data):
	print "run()"
	# update the data
#	t, y = data
#	xdata.append(t)
#	ydata.append(y)
#	xmin, xmax = ax.get_xlim()
#
#	if t >= xmax:
#		ax.set_xlim(xmin, 2*xmax)
#		ax.figure.canvas.draw()
	line.set_data(xdata, ydata)

	return line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,
				repeat=False, init_func=init)

plt.show()

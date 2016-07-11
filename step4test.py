# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:05:12 2016

@author: clay_budin

Step 4: 1D Burgers's Equation
	du/dt + u*du/dx - v*d^2u/dx^2 = 0

Plot analytic fn outside original bounds to see if it is periodic
It doesn't look like it

"""

import numpy as np
import matplotlib.pyplot as plt

# simulation constants
nx = 601
dt = 0.002
dx = (2.0*np.pi) / nx   # range is 0 to 2pi-dx to enforce periodic BC u(2pi) = u(0)
v = .02

# simlation grids - 1D
xdata, yana = [], []


# global var - current time step
ct = 0

def data_gen ():
	global ct

	# analytical soln
	t = ct * dt
	for i in range(nx):
		x = i*dx*3.0
		x = x - nx*dx
		xdata.append(x)
		p1 = np.exp(-((x-4*t)*(x-4*t))/(4.0*v*(t+1)))
		p2 = np.exp(-(x-4*t-2*np.pi)*(x-4*t-2*np.pi)/(4.0*v*(t+1)))
		y = 4.0 - 2.0*v*(((-(x-4*t)/(2.0*v*(t+1)))*p1+(-(x-4*t-2.0*np.pi)/(2.0*v*(t+1)))*p2) / (p1+p2))
		yana.append(y)


data_gen()
plt.plot(xdata, yana, 'r')

ct = 1000
data_gen()
plt.plot(xdata, yana)





# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
line1, = ax.plot([], [], lw=1)
line2, = ax.plot([], [], lw=1)

ax.grid()
ax.set_xlim(0, 1)
ax.set_ylim(0,2)


x1 = [ 0, 1 ]
y1 = [ 2, 1 ]
line1.set_data(x1, y1)

x2 = [ 0, 1 ]
y2 = [ 1, 2 ]
line2.set_data(x2, y2)

plt.show()





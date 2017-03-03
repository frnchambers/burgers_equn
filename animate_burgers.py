#!/usr/bin/python

# Double pendulum formula translated from the C code at
# http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c


import numpy as np
import matplotlib

import sys as sys


## ----------------------------- whether to save the profile or not ----------------------------- ##
save = False
name='data/animation.mp4'
for arg in sys.argv:
    if( arg=="save" ):
        save = True
        name = sys.argv[2]

if( save==True ):
    # if working remotely need to do this for saving profiles
    # ---> must be before importing matplotlib.pyplot or pylab!
    matplotlib.use('Agg')

## ---------------------------------------------------------------------------------------------- ##


import matplotlib.pyplot as plt
import matplotlib.animation as animation


ls = [ '-', '--' ]
pt = [ 'o' ]
cl = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']


print "# reading data from files...."

grid = np.loadtxt("data/grid.dat", unpack=True)
N_grid  = len(grid)
print "# .... have grid data"


time, nu = np.loadtxt("data/solution.dat", comments=["#","-----","d:"], usecols=(1,2), unpack=True)
N_steps = len(time)
print "# .... have time data"

solun = np.loadtxt("data/solution.dat", comments=["#","-----","t:"], usecols=(1,2,3), unpack=True)
u = np.reshape(solun[0], (N_steps,N_grid))
du = np.reshape(solun[1], (N_steps,N_grid))
print "# .... have solution data"


lim_x = [np.amin(grid), np.amax(grid)]
lim_u = [np.amin(u), np.amax(u)]
lim_du = [np.amin(du), np.amax(du)]


fig = plt.figure( figsize=(8, 10) )

ax_u = fig.add_subplot( 211 )
ax_du = fig.add_subplot( 212 )

#####################################################################################
################################ Set up plots domain ################################


# Common for all plots
# ax_u.set_xscale('log');  ax_u.set_yscale('log')
# ax_du.set_xscale('log'); ax_du.set_yscale('linear')

# xtics = [1.0e5,1.0e6,1.0e7,3.0e8,1.0e9,1.0e10,1.0e11,1.0e12,1.0e13,1.0e14,1.0e15,1.0e16,1.0e8]
# ax_u.xaxis.set_ticks(xtics)
# ax_du.xaxis.set_ticks(xtics)


ax_u.set_xlim( lim_x )
ax_u.set_ylim( lim_u )
ax_u.set_ylabel("u ()")

ax_du.set_xlim( lim_x )
ax_du.set_ylim( lim_u )
ax_du.set_ylabel("du/dx ()")
ax_du.set_xlabel("x ()")


# ax_T.grid();  ax_T.margins(0, 0)
# ax_DT.grid(); ax_DT.margins(0, 0)

line_u,  = ax_u.plot([], [], marker='o', color =cl[1], lw =1.5)
line_du, = ax_du.plot([], [], linestyle ='-', color =cl[4], lw =1.5)

    
time_template = 'time = %.2e s'
time_text = ax_u.text(0.05, 0.9, '', transform=ax_u.transAxes)

def init():

    result = []
    
    line_u.set_data([], [])
    line_du.set_data([], [])
    time_text.set_text('')

    result.append(line_u)
    result.append(line_du)
    result.append(time_text)

    return result


def animate(i):
    result = []

    line_u.set_data( grid, u[i])
    line_du.set_data(grid, du[i])
    time_text.set_text(time_template % (time[i]))

    result.append(line_u)
    result.append(line_du)
    result.append(time_text)
        
    return result


step_delay=5

# FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)
ani = animation.FuncAnimation(fig, animate, np.arange(0, len(time)), init_func=init, interval=step_delay, blit=True )


if( save==True ):
    ani.save(name, fps=30, codec="libx264")
else:
    plt.show()

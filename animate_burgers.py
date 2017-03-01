#!/usr/bin/python

# Double pendulum formula translated from the C code at
# http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c


import numpy as np
import matplotlib

import sys as sys


## ----------------------------- whether to save the profile or not ----------------------------- ##
save = False
name='back_evolve.mp4'
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

data_dir = "data/"

t      = np.genfromtxt(data_dir+"EVOLVE_TIMESTEP.dat", unpack=True)
y_grid = np.genfromtxt(data_dir+"EVOLVE_Y_GRIDPOINTS.dat", unpack=True)
print "# get evolution grids"

T_grid = np.genfromtxt(data_dir+"EVOLVE_TEMP.dat", unpack=True).transpose()
print "# .... got T"
DT_grid = np.genfromtxt(data_dir+"EVOLVE_DLGTDLGY.dat", unpack=True).transpose()
print "# .... got DT"


fig = plt.figure( figsize=(10, 10) )

ax_u = fig.add_subplot( 211 )
ax_DT = fig.add_subplot( 212 )


#####################################################################################
################################ Set up plots domain ################################


# Common for all plots
ax_T.set_xscale('log');  ax_T.set_yscale('log')
ax_DT.set_xscale('log'); ax_DT.set_yscale('linear')



maxx =np.amax(y_grid);   minx =np.amin(y_grid)
maxyT=np.amax(T_grid);   minyT=np.amin(T_grid)
maxyDT=np.max(DT_grid); minyDT=np.min(DT_grid)
# maxyr=np.amax(r_grid);   minyr=np.amin(r_grid)
# maxyCP=np.amax(CP_grid); minyCP=np.amin(CP_grid)
# maxyK=np.amax(K_grid);   minyK=np.amin(K_grid)



# Default is full intagration region
xlimits    =[minx, maxx]
ylimits_T  =[minyT, maxyT]
ylimits_DT =[minyDT, maxyDT]
# ylimits_r  =[minyr, maxyr]
# ylimits_CP =[minyCP, maxyCP]
# ylimits_K  =[minyK, maxyK]


xtics = [1.0e5,1.0e6,1.0e7,3.0e8,1.0e9,1.0e10,1.0e11,1.0e12,1.0e13,1.0e14,1.0e15,1.0e16,1.0e8]
ax_T.xaxis.set_ticks(xtics)
ax_DT.xaxis.set_ticks(xtics)
# ax_r.xaxis.set_ticks(xtics)
# ax_CP.xaxis.set_ticks(xtics)
# ax_K.xaxis.set_ticks(xtics)


ax_T.set_xlim( xlimits )
ax_T.set_ylim( ylimits_T )
ax_T.set_ylabel("Temperature (K)")
ax_T.margins(0, 0)

ax_DT.set_xlim( xlimits )
ax_DT.set_ylim( ylimits_DT )
ax_DT.set_ylabel("dlogT/dlogy ()")
ax_DT.set_xlabel("Column Depth = P/g (g cm^-2)")

# ax_r.set_xlim( xlimits )
# ax_r.set_ylim( ylimits_r )
# ax_r.set_ylabel("Density (g cm^-3)")
# ax_r.margins(0, 0)

# ax_CP.set_xlim( xlimits )
# ax_CP.set_ylim( ylimits_CP )
# ax_CP.set_ylabel("Heat Capacity (g cm^-3)")
# ax_CP.margins(0, 0)

# ax_K.set_xlim( xlimits )
# ax_K.set_ylim( ylimits_K )
# ax_K.set_ylabel("Opacity (cm^2 g^-1)")
# ax_K.margins(0, 0)


ax_T.grid();  ax_T.margins(0, 0)
ax_DT.grid(); ax_DT.margins(0, 0)
# ax_r.grid();  ax_r.margins(0, 0)
# ax_CP.grid(); ax_CP.margins(0, 0)
# ax_K.grid();  ax_K.margins(0, 0)

line_T,  = ax_T.plot([], [], marker='o', color =cl[1], lw =1.5)
line_DT, = ax_DT.plot([], [], linestyle ='-', color =cl[4], lw =1.5)
# line_r,  = ax_r.plot([], [], linestyle  ='-', color =cl[5], lw =1.5)
# line_CP, = ax_CP.plot([], [], linestyle ='-', color =cl[5], lw =1.5)
# line_K,  = ax_K.plot([], [], linestyle  ='-', color =cl[4], lw =1.5)

    
time_template = 'time = %.2e s'
time_text = ax_T.text(0.05, 0.9, '', transform=ax_T.transAxes)

def init():

    result = []
    
    line_T.set_data([], [])
    line_DT.set_data([], [])
    result.append(line_T)
    result.append(line_DT)

    time_text.set_text('')
    result.append(time_text)

    return result


def animate(i):

    xcoords = y_grid
    result = []

    line_T.set_data(xcoords, T_grid[i])
    line_DT.set_data(xcoords, DT_grid[i])
    result.append(line_T)
    result.append(line_DT)

    time_text.set_text(time_template % (t[i]))
    result.append(time_text)
        
    return result


step_delay=1

# FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)
ani = animation.FuncAnimation(fig, animate, np.arange(1, len(t)), init_func=init, interval=step_delay, blit=True )


if( save==True ):
    ani.save(name, fps=30, codec="libx264")
else:
    plt.show()

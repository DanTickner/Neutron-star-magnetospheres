'''
Smooth_increase_and_decrease_1D_multiple_v1.py

Plot a function in one dimension which smoothly increases toward a set value,
and then smoothly decreases to zero.
This will be used for implementing smoothing functions for magnetospheric twists,
avoiding potential jump discontinuities at the coordinate boundaries.

This function is two exponential decays, rising from zero to some value A and then back to zero.
Plot multiple graphs with differing values of B and C.
'''

import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
x_min     = -4  # Minimum value of x to plot.
x_max     = 4   # Maximum value of x to plot.
n_points  = 100 # Number of gridpoints to plot.
x_0       = 3   # Value of x at which the function should hit zero.
A         = 2   # Value of f(x) far from the ramping regions.
B         = 3   # Speed at which the function ramps up (larger B is faster, less smooth ramp-up).



#----- Define functions -----
def f( x ):
	if x >= 0:
		return A * ( 1.0 - np.exp( - B * x ) )
	else:
		return 0


#----- Build array of function values -----
x_list = np.linspace( x_min, x_max, n_points )

f_list = np.zeros( n_points ) # Can't do f = function(x) because that doesn't work with if-statements.

for n in range( n_points ):
	f_list[n] = f( x_list[n] )




#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( x_list, f_list )
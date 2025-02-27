'''
Smooth_increase_and_decrease_1D_multiple_v2.py

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
x_1       = -2  # Minimum value of x for which function is defined.
x_2       = 3   # Maximum value of x for which function is defined.
n_points  = 1000 # Number of gridpoints to plot.
A         = 2   # Value of f(x) far from the ramping regions.
B         = 3   # Speed at which the function ramps up (larger B is faster, less smooth ramp-up).
C         = 3   # Speed at which the function ramps up (larger B is faster, less smooth ramp-up).

fig_filepath = "../Writeup/Figures/Numerical-Method/20241025_Smooth_increase_and_decrease_exponential_1D_multiple"



#----- Define functions -----
def f( x ):
	if x >= x_1 and x <= x_2:
		return A * ( 1.0 - np.exp( - B * ( x - x_1 ) ) ) * ( 1.0 - np.exp( - C * ( -x + x_2 ) ) )
	else:
		return 0


#----- Build array of function values -----
x_list = np.linspace( x_min, x_max, n_points )

f_list_a = np.zeros( n_points ) # Can't do f = function(x) because that doesn't work with if-statements.
f_list_b = np.zeros( n_points )
f_list_c = np.zeros( n_points )

for n in range( n_points ):
	f_list_a[n] = f( x_list[n] )

B = 6
C = 6

for n in range( n_points ):
	f_list_b[n] = f( x_list[n] )
	
B = 100
C = 100

for n in range( n_points ):
	f_list_c[n] = f( x_list[n] )




#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( x_list, f_list_a, color="k"   , ls="-" , label=r"$B=C=3$" )
ax.plot( x_list, f_list_b, color="b"   , ls="--", label=r"$B=C=6$" )
ax.plot( x_list, f_list_c, color="grey", ls=":" , label=r"$B,C\to\infty$" )

ax.set_xlim( x_min, x_max )
ax.set_ylim( -0.05, A*1.05 )

xticks_values = [ x_1     , 0     , x_2      ]
xticks_labels = [ r"$x_1$", r"$0$", r"$x_2$" ]

yticks_values = [ 0     , A      ]
yticks_labels = [ r"$0$", r"$A$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.legend( loc="best" )

fig.tight_layout()

plt.savefig( fig_filepath )
print( f"Figure saved:\t{ fig_filepath }" )
'''
Smooth_increase_and_decrease_1D_linear.py

Plot a function in one dimension which smoothly increases toward a set value,
and then smoothly decreases to zero.
This will be used for implementing smoothing functions for magnetospheric twists,
avoiding potential jump discontinuities at the coordinate boundaries.

This function is two exponential decays, rising from zero to some value A and then back to zero.
Plot only one graph (one set of values B,C) and fully label it.
'''

import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
x_min     = -4   # Minimum value of x to plot.
x_max     = 4    # Maximum value of x to plot.
x_L       = -3   # Lower limit of x for which function is nonzero.
x_U       = 3    # Upper limit of x for which function is nonzero.
n_points  = 1000 # Number of gridpoints to plot.
A         = 2    # Desired ramped-up function value.

x_A_minus = -2   # Smallest value of x at which f(x)=A.
x_A_plus  = 1.5  # Largest  value of x at which f(x)=A.

fig_filepath = "../Writeup/Figures/Numerical-Method/20250217_Smooth_increase_and_decrease_linear_1D"



#----- Define functions -----
def f( x ):
	if x >= x_L and x <= x_A_minus:
		return ( A / ( x_A_minus - x_L ) ) * ( x - x_L )
	elif x >= x_A_minus and x <= x_A_plus:
		return A
	elif x >= x_A_plus and x <= x_U:
		return ( A / ( x_A_plus - x_U ) ) * ( x - x_U )
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

ax.plot( x_list, f_list, color="k"   , ls="-" )

ax.set_xlim( x_min, x_max )
ax.set_ylim( -0.05, A*1.05 )

xticks_values = [ x_L            , x_A_minus     , x_A_plus     , x_U             ]
xticks_labels = [ r"$x_{\rm{L}}$", r"$x_{A-}$", r"$x_{A+}$", r"$x_{\rm{U}}$" ]

yticks_values = [ 0     , A      ]
yticks_labels = [ r"$0$", r"$A$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.plot( [ x_A_minus, x_A_minus ] , [ -1 , f(x_A_minus) ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_A_plus , x_A_plus  ] , [ -1 , f(x_A_plus ) ], color="grey", ls="--", lw=0.5 )

ax.axhline( A, color="grey", ls="--", lw=1 )

fig.tight_layout()

plt.savefig( fig_filepath )
print( f"Figure saved:\t{ fig_filepath }" )
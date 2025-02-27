'''
Smooth_increase_and_decrease_1D_exponential_single.py

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
p         = 0.75 # Control ramp-up speed: fraction of maximum value to reach at x_p.
x_p       = -2.4 # Control ramp-up speed: value of x at which f(x_p)=pA. Need x_L < x_p x_p < < x_U.
q         = 0.9  # Control ramp-down speed: fraction of maximum value to reach at x_q.
x_q       = 2.3  # Control ramp-down speed: value of x at which f(x_q)=qA. Need x_L < x_p < x_q < x_U.

fig_filepath = "../Writeup/Figures/Numerical-Method/20241025_Smooth_increase_and_decrease_exponential_1D_single"



#----- Define functions -----
B = - np.log( 1.0 - p ) / (   x_p - x_L )
C = - np.log( 1.0 - q ) / ( - x_q + x_U )

def f( x ):
	if x >= x_L and x <= x_U:
		return A * ( 1.0 - np.exp( - B * ( x - x_L ) ) ) * ( 1.0 - np.exp( - C * ( -x + x_U ) ) )
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

xticks_values = [ x_L            , x_p     , x_q     , x_U             ]
xticks_labels = [ r"$x_{\rm{L}}$", r"$x_p$", r"$x_q$", r"$x_{\rm{U}}$" ]

yticks_values = [ 0     , p*A    , q*A    , A      ]
yticks_labels = [ r"$0$", r"$pA$", r"$qA$", r"$A$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.plot( [ x_min, x_p ] , [ p*A, p*A    ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p  , x_p ] , [ -1 , f(x_p) ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min, x_q ] , [ q*A, q*A    ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q  , x_q ] , [ -1 , f(x_q) ], color="grey", ls="--", lw=0.5 )

ax.axhline( A, color="grey", ls="--", lw=1 )

fig.tight_layout()

plt.savefig( fig_filepath )
print( f"Figure saved:\t{ fig_filepath }" )
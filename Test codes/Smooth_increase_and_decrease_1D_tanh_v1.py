'''
Smooth_increase_and_decrease_1D_tanh_v1.py

Plot a function in one dimension which smoothly increases toward a set value,
and then smoothly decreases to zero.
This will be used for implementing smoothing functions for magnetospheric twists,
avoiding potential jump discontinuities at the coordinate boundaries.

This function is a hyperbolic tangent rising from zero to some value A.
'''

import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
x_min     = -1  # Minimum value of x to plot.
x_max     = 3   # Maximum value of x to plot.
n_points  = 1000 # Number of gridpoints to plot.

A         = 2    # Value of f(x) far from the ramping regions.
p         = 0.2 # Fraction of maximum value to obtain at x = x_p_minus.
x_p_minus = 1.1  # Value of x at which f(x) = pA.
x_p_plus  = 1.9  # Value of x at which f(x) = (1-p)A.

x_p = 0.5 * ( x_p_minus + x_p_plus )
B = ( 2.0 / ( x_p_minus - x_p_plus ) ) * np.arctanh( 2.0 * p - 1 )
print( B )




#----- Define functions -----
def f( x ):
	return 0.5 * A * ( np.tanh( B * ( x - x_p ) ) + 1 )


print( f( x_p_minus ) )
print( f( x_p_plus  ) )
print( x_p )
print( f( x_p ) )

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

xticks_values = [ x_p_minus  , x_p     , x_p_plus    ]
xticks_labels = [ r"$x_{p-}$", r"$x_p$", r"$x_{p+}$" ]

yticks_values = [ 0     , p*A     , 0.5*A          , (1-p)*A    , A      ]
yticks_labels = [ r"$0$", r"$pA$", r"$\frac{A}{2}$", r"$(1-p)A$", r"$A$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.plot( x_list, f_list )

ax.plot( [ x_min    , x_p_minus ] , [ p*A    , p*A          ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p_minus, x_p_minus ] , [ 0      , f(x_p_minus) ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_p       ] , [ 0.5*A  , 0.5*A        ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p      , x_p       ] , [ 0      , f(x_p)       ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_p_plus  ] , [ (1-p)*A, (1-p)*A      ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p_plus , x_p_plus  ] , [ 0      , f(x_p_plus)  ], color="grey", ls="--", lw=0.5 )

ax.set_xlim( x_min, x_max )
ax.set_ylim( -0.02, A*1.01 )

fig.tight_layout()
'''
Smooth_increase_and_decrease_1D_tanh_v2.py

Plot a function in one dimension which smoothly increases toward a set value,
and then smoothly decreases to zero.
This will be used for implementing smoothing functions for magnetospheric twists,
avoiding potential jump discontinuities at the coordinate boundaries.

This function is a hyperbolic tangent falling from some value A to zero.
'''

import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
x_min     = -1  # Minimum value of x to plot.
x_max     = 3   # Maximum value of x to plot.
n_points  = 1000 # Number of gridpoints to plot.

A         = 2    # Value of f(x) far from the ramping regions.
q         = 0.2 # Fraction of maximum value to obtain at x = x_q_minus.
x_q_minus = 1.1  # Value of x at which f(x) = (1-q)A.
x_q_plus  = 1.9  # Value of x at which f(x) = qA.

x_q = 0.5 * ( x_q_minus + x_q_plus )
B = ( 2.0 / ( x_q_minus - x_q_plus ) ) * np.arctanh( 2.0 * q - 1 )
print( B )




#----- Define functions -----
def f( x ):
	return 0.5 * A * ( np.tanh( B * ( -x + x_q ) ) + 1 )

print( f( x_q_minus ) )
print( f( x_q_plus  ) )
print( x_q )
print( f( x_q ) )

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

xticks_values = [ x_q_minus  , x_q     , x_q_plus    ]
xticks_labels = [ r"$x_{q-}$", r"$x_q$", r"$x_{q+}$" ]

yticks_values = [ 0     , (1-q)*A     , 0.5*A            , q*A    , A      ]
yticks_labels = [ r"$0$", r"$(1-q)A$"  , r"$\frac{A}{2}$", r"$qA$", r"$A$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.plot( x_list, f_list )

ax.plot( [ x_min    , x_q_plus  ] , [ q*A    , q*A          ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q_plus , x_q_plus  ] , [ 0      , f(x_q_plus)  ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_q       ] , [ 0.5*A  , 0.5*A        ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q      , x_q       ] , [ 0      , f(x_q)       ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_q_minus ] , [ (1-q)*A, (1-q)*A      ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q_minus, x_q_minus ] , [ 0      , f(x_q_minus) ], color="grey", ls="--", lw=0.5 )

ax.set_xlim( x_min, x_max )
ax.set_ylim( -0.02, A*1.01 )

fig.tight_layout()
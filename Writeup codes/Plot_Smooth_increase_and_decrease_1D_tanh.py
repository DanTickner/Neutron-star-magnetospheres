'''
Smooth_increase_and_decrease_1D_tanh_v3.py

Plot a function in one dimension which smoothly increases toward a set value,
and then smoothly decreases to zero.
This will be used for implementing smoothing functions for magnetospheric twists,
avoiding potential jump discontinuities at the coordinate boundaries.

This function is two hyperbolic tangents, rising from zero to some value A and then back to zero.
'''

import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
fig_filename = "../Writeup/Figures/Numerical-Method/20241024_Smooth_increase_and_decrease_tanh_1D"

x_min     = -0.5 # Minimum value of x to plot.
x_max     = 8    # Maximum value of x to plot.
n_points  = 1000 # Number of gridpoints to plot.

x_L       = 0    # Smallest value of x for which the smoothing function is nonzero.
x_U       = 7    # Largest  value of x for which the smoothing function is nonzero.

A         = 2    # Value of f(x) far from the ramping regions.
p         = 0.2  # Fraction of maximum value to obtain at x = x_p_minus.
x_p_minus = 1.1  # Value of x at which f(x) = pA.
x_p_plus  = 1.9  # Value of x at which f(x) = (1-p)A.

q         = 0.1  # Fraction of maximum value to obtain at x = x_q_minus.
x_q_minus = 4.8  # Value of x at which f(x) = (1-q)A.
x_q_plus  = 6    # Value of x at which f(x) = qA.

x_p = 0.5 * ( x_p_minus + x_p_plus )
B = ( 2.0 / ( x_p_minus - x_p_plus ) ) * np.arctanh( 2.0 * p - 1 )
print( B )

x_q = 0.5 * ( x_q_minus + x_q_plus )
C = ( 2.0 / ( x_q_minus - x_q_plus ) ) * np.arctanh( 2.0 * q - 1 )
print( C )




#----- Define functions -----
def f( x ):
	if x > x_L and x < x_U:
		return 0.25 * A * ( np.tanh( B * ( x - x_p ) ) + 1 ) * ( np.tanh( C * ( -x + x_q ) ) + 1 )
	else:
		return 0


print( f( x_p_minus ) )
print( f( x_p_plus  ) )
print( x_p )
print( f( x_p ) )

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

fig = plt.figure( figsize=(8,5) )
ax  = fig.add_subplot( 111 )

xticks_values = [ x_L           , x_p_minus  , x_p     , x_p_plus    ] + [ x_q_minus  , x_q     , x_q_plus   , x_U             ]
xticks_labels = [ r"$x_{\rm{L}}$", r"$x_{p-}$", r"$x_p$", r"$x_{p+}$" ] + [ r"$x_{q-}$", r"$x_q$", r"$x_{q+}$", r"$x_{\rm{U}}$" ]

yticks_values = [ 0     , p*A     , 0.5*A          , (1-p)*A    , A      ] + [ (1-q)*A     , q*A     ]
yticks_labels = [ r"$0$", r"$pA$", r"$\frac{A}{2}$", r"$(1-p)A$", r"$A$" ] + [ r"$(1-q)A$" , r"$qA$" ]

ax.set_xticks( xticks_values, xticks_labels )
ax.set_yticks( yticks_values, yticks_labels )

ax.plot( x_list, f_list )


ax.plot( [ x_min    , x_p_minus ] , [ p*A    , p*A          ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p_minus, x_p_minus ] , [ -1     , f(x_p_minus) ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_p_plus  ] , [ (1-p)*A, (1-p)*A      ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p_plus , x_p_plus  ] , [ -1     , f(x_p_plus)  ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_q_plus  ] , [ q*A    , q*A          ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q_plus , x_q_plus  ] , [ -1     , f(x_q_plus)  ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min    , x_q_minus ] , [ (1-q)*A, (1-q)*A      ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q_minus, x_q_minus ] , [ -1     , f(x_q_minus) ], color="grey", ls="--", lw=0.5 )

ax.plot( [ x_min, x_p, x_q ] , [ 0.5*A, 0.5*A , 0.5*A ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_p  , x_p      ] , [ -1   , f(x_p)        ], color="grey", ls="--", lw=0.5 )
ax.plot( [ x_q  , x_q      ] , [ -1   , f(x_q)        ], color="grey", ls="--", lw=0.5 )

ax.axhline( A, color="grey", ls="--", lw=1 )

ax.set_xlim( x_min, x_max )
ax.set_ylim( -0.02, A*1.01 )

fig.tight_layout()

plt.savefig( fig_filename )
print( f"\nFigure saved:\t{ fig_filename }" )
'''
Plot_Parabola_From_3_Points_and_Derivatives.py

It is always possible to construct a parabola from three gridpoints.
ax^2+bx+c
given (x0,y0), (x1,y1) and (x2,y2).
However, the derivatives computed from the parabola may not match the original gridpoints.

Here, we test a formula for a, b and c given the gridpoints.
Use three random gridpoints.
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
import random


#----- Define variables -----
x_min = 0.5
x_max = 2.5

y_min = 0
y_max = 10

# Ensure the x-values are in ascending order so that we can pick out the central value.
xlist = [ x_min + ( x_max - x_min ) * random.random() for i in range(3) ]
xlist.sort()

x0 = xlist[0]
x1 = xlist[1]
x2 = xlist[2]

n_points = 1000

fig_filename = "Plot_Parabola_From_3_Points_and_Derivatives_20241223"




#----- Define functions -----
def f( x ):
	# Change this function to taste.
	return np.sin( x )

def f_dx( x ):
	# Change this function to taste.
	return np.cos( x )

def f_dx2( x ):
	# Change this function to taste.
	return -np.sin( x )

y0 = f( x0 )
y1 = f( x1 )
y2 = f( x2 )



#----- Calculate parabola parameters -----
#print( x0**2, x1**2, x2**2 )	# So can check det_A against online matrix determinant calculator.

det_A = x0**2 * ( x1 - x2 ) + x0 * ( x2**2 - x1**2 ) + x1 * x2 * ( x1 - x2 )

a = ( ( x1 - x2 ) * y0 + ( x2 - x0 ) * y1 + ( x0 - x1 ) * y2 ) / det_A
b = ( ( x2**2 - x1**2 ) * y0 + ( x0**2 - x2**2 ) * y1 + ( x1**2 - x0**2 ) * y2 ) / det_A
c = ( x1 * x2 * ( x1 - x2 ) * y0 + x0 * x2 * ( x2 - x0 ) * y1 + x0 * x1 * ( x0 - x1 ) * y2 ) / det_A

print( f"a = { a }" )
print( f"b = { b }" )
print( f"c = { c }" )
print()

f_dx_exact          = f_dx( x1 )
f_dx_from_parabola  = 2.0 * a * x1 + b
f_dx2_exact         = f_dx2( x1 )
f_dx2_from_parabola = 2.0 * a

print( f"f'(x_1) exact         :\t{ f_dx_exact         }" )
print( f"f'(x_1) from parabola :\t{ f_dx_from_parabola }" )
print( f"f''(x_1) exact        :\t{ f_dx2_exact        }" )
print( f"f''(x_1) from parabola:\t{ f_dx2_from_parabola}" )


#----- Build list of parabola values -----
x_list = np.linspace( x_min, x_max, n_points )
y_list_exact         = np.zeros( n_points )
y_list_from_parabola = np.zeros( n_points )

for i, x in enumerate( x_list ):
	y_list_exact        [i] = f( x )
	y_list_from_parabola[i] = a * x*x + b * x + c




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Plot graph -----
fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.scatter( [x0,x1,x2], [y0,y1,y2], color="k" )
ax.plot( x_list, y_list_exact        , color="r", ls="-", label="Function (Exact)" )
ax.plot( x_list, y_list_from_parabola, color="b", ls="--", label="Function (Parabola)" )

#----- Plot gradient lines -----
# Create a new gridpoint ( x_b, y_b ).
# Calculate y_b from the gradient at x_1.
# Use the straight-line formula y-y1 = m( x-x1 ), so y = m( x-x1 ) + y1.
gradient_line_x_b = x1 + 0.5
gradient_line_y_b_exact         = f_dx_exact         * ( gradient_line_x_b - x1 ) + y1
gradient_line_y_b_from_parabola = f_dx_from_parabola * ( gradient_line_x_b - x1 ) + y1

ax.plot( [ x1, gradient_line_x_b ], [ y1, gradient_line_y_b_exact         ], ls=":", color="red"  , label="Tangent (Exact)" )
ax.plot( [ x1, gradient_line_x_b ], [ y1, gradient_line_y_b_from_parabola ], ls="-.", color="grey", label="Tangent (Parabola)" )

ax.legend( loc="best" )

fig.tight_layout()
plt.savefig( fig_filename )
print( f"Figure saved:\t{ fig_filename }" )
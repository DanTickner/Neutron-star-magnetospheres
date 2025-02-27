'''
Plot_Parabola_From_3_Points_Random.py

It is always possible to construct a parabola from three gridpoints.
ax^2+bx+c
given (x0,y0), (x1,y1) and (x2,y2).

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

x0 = x_min + ( x_max - x_min ) * random.random()
y0 = y_min + ( y_max - y_min ) * random.random()

x1 = x_min + ( x_max - x_min ) * random.random()
y1 = y_min + ( y_max - y_min ) * random.random()

x2 = x_min + ( x_max - x_min ) * random.random()
y2 = y_min + ( y_max - y_min ) * random.random()

n_points = 1000

fig_filename = "Plot_Parabola_From_3_Points_Random_20241223"



#----- Calculate parabola parameters -----
#print( x0**2, x1**2, x2**2 )	# So can check det_A against online matrix determinant calculator.

det_A = x0**2 * ( x1 - x2 ) + x0 * ( x2**2 - x1**2 ) + x1 * x2 * ( x1 - x2 )

a = ( ( x1 - x2 ) * y0 + ( x2 - x0 ) * y1 + ( x0 - x1 ) * y2 ) / det_A
b = ( ( x2**2 - x1**2 ) * y0 + ( x0**2 - x2**2 ) * y1 + ( x1**2 - x0**2 ) * y2 ) / det_A
c = ( x1 * x2 * ( x1 - x2 ) * y0 + x0 * x2 * ( x2 - x0 ) * y1 + x0 * x1 * ( x0 - x1 ) * y2 ) / det_A

print( f"a = { a }" )
print( f"b = { b }" )
print( f"c = { c }" )




#----- Build list of parabola values -----
x_list = np.linspace( x_min, x_max, n_points )
y_list = np.zeros( n_points )

for i, x in enumerate( x_list ):
	y_list[i] = a * x*x + b * x + c




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Plot graph -----
fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.scatter( [x0,x1,x2], [y0,y1,y2], color="k" )
ax.plot( x_list, y_list, color="b", ls="--" )

fig.tight_layout()
plt.savefig( fig_filename )
print( f"Figure saved:\t{ fig_filename }" )
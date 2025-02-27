'''
Plot_Finite_Difference_Visualisation.py

Plot a small number of discrete gridpoints and colour/label them according to which method you can use.
Don't actually calculate the derivatives since this has been done already in other codes.
This code is purely for the production of an image for a report or presentation.
'''

import matplotlib.pyplot as plt
import numpy             as np


#----- Define function to plot -----
def f(x):
	return x**(-3.0)


#----- Define variables -----
x_min = 1.0
x_max = 2.0
n     = 10				# Number of points to plot

x = np.linspace( x_min, x_max, n )
y = f( x )



#----- Setup graphs -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":20 })

i_forward_backward = [ 0, -1 ]
#i_1_offset         = [ 1, -2 ]
#i_2_offset         = [ 2, -3 ]
i_symmetric        = [ 1, 2, 3, 4, 5, 6, 7, 8 ]

x_forward_backward = [ x[i] for i in i_forward_backward ]
y_forward_backward = [ y[i] for i in i_forward_backward ]
#x_1_offset         = [ x[i] for i in i_1_offset         ]
#y_1_offset         = [ y[i] for i in i_1_offset         ]
#x_2_offset         = [ x[i] for i in i_2_offset         ]
#y_2_offset         = [ y[i] for i in i_2_offset         ]
x_symmetric        = [ x[i] for i in i_symmetric        ]
y_symmetric        = [ y[i] for i in i_symmetric        ]


#----- Plot graph with one step -----
fig = plt.figure( figsize=(7,5) )
ax  = fig.add_subplot( 111 )

ax.scatter( x_forward_backward, y_forward_backward, color="b"         , label="Forward and backward" )
#ax.scatter( x_1_offset        , y_1_offset        , color="orange"    , label="1-offset"             )
#ax.scatter( x_2_offset        , y_2_offset        , color="lightgreen", label="2-offset"             )
ax.scatter( x_symmetric       , y_symmetric       , color="k"         , label="Symmetric"            )

ax.legend( loc="upper right" )
ax.set_title( r"Order $h^2$" )

xticks_values = [ x_min, x_max ]
xticks_labels = [ r"$x_0$", r"$x_{N-1}$" ]
yticks_values = [ f(x_min), f(x_max) ]
yticks_labels = [ f"$f(x_0)$", r"$f(x_{N-1})$" ]
plt.xticks( xticks_values, xticks_labels )
plt.yticks( yticks_values, yticks_labels )

fig.tight_layout()
plt.savefig( "../Figures/Finite_Difference_Visualisation_h_to_the_2" )
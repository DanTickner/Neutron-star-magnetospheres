import matplotlib.pyplot as plt
import numpy             as np
import math

# This code aims to develop a guessed form of atan2( sin(x), cos(x) ) in terms of x.

#----- Define guessed atan2 functions -----
def atan2_sin_cos( x ):
	# Guessed atan2( y, -x )
	
	#return 2.0 * ( ( ( x - np.pi ) / 2.0 ) % np.pi ) - np.pi
	#return ( x + np.pi ) % (2.0*np.pi) - np.pi
	return x


def atan2_build1( x ):
	# One option that the correct result could be.
	return - x

def atan2_build2( x ):
	# One option that the correct result could be.
	return - x + 2.0*np.pi


#----- Define variables -----
x_min    = -4
x_max    = 4
n_points = 1000


#----- Build lists of datapoints -----
x_plot = np.linspace( x_min, x_max, n_points )

atan2_plot        = [ math.atan2   ( np.sin(x), np.cos(x) ) for x in x_plot ]

atan2_plot_guess  = [ atan2_sin_cos( x                    ) for x in x_plot ]

atan2_plot_build1 = [ atan2_build1 ( x                    ) for x in x_plot ]

atan2_plot_build2 = [ atan2_build2 ( x                    ) for x in x_plot ]


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 16 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( x_plot, atan2_plot       , ls="-"  , label="Exact"  )
#ax.plot( x_plot, atan2_plot_build1, ls="--" , label="Build1" )
#ax.plot( x_plot, atan2_plot_build2, ls="-." , label="Build2" )
ax.plot( x_plot, atan2_plot_guess , ls=":"  , label="Guess"  )

ax.legend( loc="best" )
ax.set_xlabel( r"$x$" )
ax.set_ylabel( r"${\rm atan2}(x,y)$" )
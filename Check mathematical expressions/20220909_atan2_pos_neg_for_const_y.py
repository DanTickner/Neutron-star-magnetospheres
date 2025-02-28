import matplotlib.pyplot as plt
import numpy             as np
import math

# This code aims to develop a guessed form of atan2(y,-x) in terms of atan2(y,x).

#----- Define guessed atan2 functions -----
def atan2_pos_neg( y, x ):
	# Guessed atan2( y, -x )
	
	if   y >= 0:
		return - math.atan2( y, x ) + math.pi
	elif y < 0:
		return - math.atan2( y, x ) - math.pi
	else:
		return 3.5	# Value slightly larger than pi should be easily visible as an error.


#----- Define variables -----
x_min    = -3
x_max    = 3
y_neg    = -3
y_nil    = 0
y_pos    = 3
n_points = 1000


#----- Build lists of datapoints -----
x_plot = np.linspace( x_min, x_max, n_points )

atan2_plot_y_pos        = [  math.atan2   ( y_pos, -x )       for x in x_plot ]
atan2_plot_y_nil        = [  math.atan2   ( y_nil, -x )       for x in x_plot ]
atan2_plot_y_neg        = [  math.atan2   ( y_neg, -x )       for x in x_plot ]

atan2_plot_y_pos_guess  = [  atan2_pos_neg( y_pos,  x )       for x in x_plot ]
atan2_plot_y_nil_guess  = [  atan2_pos_neg( y_nil,  x )       for x in x_plot ]
atan2_plot_y_neg_guess  = [  atan2_pos_neg( y_neg,  x )       for x in x_plot ]

atan2_plot_y_pos_build1 = [ -math.atan2   ( y_pos,  x )+np.pi for x in x_plot ]
atan2_plot_y_nil_build1 = [ -math.atan2   ( y_nil,  x )+np.pi for x in x_plot ]
atan2_plot_y_neg_build1 = [ -math.atan2   ( y_neg,  x )+np.pi for x in x_plot ]

atan2_plot_y_pos_build2 = [ -math.atan2   ( y_pos,  x )-np.pi for x in x_plot ]
atan2_plot_y_nil_build2 = [ -math.atan2   ( y_nil,  x )-np.pi for x in x_plot ]
atan2_plot_y_neg_build2 = [ -math.atan2   ( y_neg,  x )-np.pi for x in x_plot ]

#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 16 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( x_plot, atan2_plot_y_pos       , color="green" , ls="-" , label=r"$y>0$, exact"  )
ax.plot( x_plot, atan2_plot_y_nil       , color="blue"  , ls="-" , label=r"$y=0$, exact"  )
ax.plot( x_plot, atan2_plot_y_neg       , color="orange", ls="-" , label=r"$y<0$, exact"  )

ax.plot( x_plot, atan2_plot_y_pos_guess , color="green" , ls=":" , label=r"$y>0$, guess"  )
ax.plot( x_plot, atan2_plot_y_nil_guess , color="blue"  , ls=":" , label=r"$y=0$, guess"  )
ax.plot( x_plot, atan2_plot_y_neg_guess , color="orange", ls=":" , label=r"$y<0$, guess"  )

#ax.plot( x_plot, atan2_plot_y_pos_build1, color="green" , ls="--", label=r"$y>0$, build1" )
#ax.plot( x_plot, atan2_plot_y_nil_build1, color="blue"  , ls="--", label=r"$y=0$, build1" )
#ax.plot( x_plot, atan2_plot_y_neg_build1, color="orange", ls="--", label=r"$y<0$, build1" )

#ax.plot( x_plot, atan2_plot_y_pos_build2, color="green" , ls="-.", label=r"$y>0$, build2" )
#ax.plot( x_plot, atan2_plot_y_nil_build2, color="blue"  , ls="-.", label=r"$y=0$, build2" )
#ax.plot( x_plot, atan2_plot_y_neg_build2, color="orange", ls="-.", label=r"$y<0$, build2" )

ax.legend( loc="best" )
ax.set_xlabel( r"$x$" )
ax.set_ylabel( r"${\rm atan2}(x,y)$" )
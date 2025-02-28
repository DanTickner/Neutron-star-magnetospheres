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
y_min    = -3
y_max    = 3
x_neg    = -3
x_nil    = 0
x_pos    = 3
n_points = 1000


#----- Build lists of datapoints -----
y_plot = np.linspace( y_min, y_max, n_points )

atan2_plot_x_pos        = [  math.atan2   ( y, -x_pos )       for y in y_plot ]
atan2_plot_x_nil        = [  math.atan2   ( y, -x_nil )       for y in y_plot ]
atan2_plot_x_neg        = [  math.atan2   ( y, -x_neg )       for y in y_plot ]

atan2_plot_x_pos_guess  = [  atan2_pos_neg( y,  x_pos )       for y in y_plot ]
atan2_plot_x_nil_guess  = [  atan2_pos_neg( y,  x_nil )       for y in y_plot ]
atan2_plot_x_neg_guess  = [  atan2_pos_neg( y,  x_neg )       for y in y_plot ]

atan2_plot_x_pos_build1 = [ -math.atan2   ( y,  x_pos )+np.pi for y in y_plot ]
atan2_plot_x_nil_build1 = [ -math.atan2   ( y,  x_nil )+np.pi for y in y_plot ]
atan2_plot_x_neg_build1 = [ -math.atan2   ( y,  x_neg )+np.pi for y in y_plot ]

atan2_plot_x_pos_build2 = [ -math.atan2   ( y,  x_pos )-np.pi for y in y_plot ]
atan2_plot_x_nil_build2 = [ -math.atan2   ( y,  x_nil )-np.pi for y in y_plot ]
atan2_plot_x_neg_build2 = [ -math.atan2   ( y,  x_neg )-np.pi for y in y_plot ]

#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 16 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

#ax.plot( y_plot, atan2_plot_x_pos       , color="green" , ls="-" , label=r"$x>0$, exact"  )
#ax.plot( y_plot, atan2_plot_x_nil       , color="blue"  , ls="-" , label=r"$x=0$, exact"  )
ax.plot( y_plot, atan2_plot_x_neg       , color="orange", ls="-" , label=r"$x<0$, exact"  )

#ax.plot( y_plot, atan2_plot_x_pos_guess , color="green" , ls=":" , label=r"$x>0$, guess"  )
#ax.plot( y_plot, atan2_plot_x_nil_guess , color="blue"  , ls=":" , label=r"$x=0$, guess"  )
#ax.plot( y_plot, atan2_plot_x_neg_guess , color="orange", ls=":" , label=r"$x<0$, guess"  )

#ax.plot( y_plot, atan2_plot_x_pos_build1, color="green" , ls="--", label=r"$x>0$, build1" )
#ax.plot( y_plot, atan2_plot_x_nil_build1, color="blue"  , ls="--", label=r"$x=0$, build1" )
#ax.plot( y_plot, atan2_plot_x_neg_build1, color="orange", ls="--", label=r"$x<0$, build1" )

#ax.plot( y_plot, atan2_plot_x_pos_build2, color="green" , ls="-.", label=r"$yx0$, build2" )
#ax.plot( y_plot, atan2_plot_x_nil_build2, color="blue"  , ls="-.", label=r"$x=0$, build2" )
ax.plot( y_plot, atan2_plot_x_neg_build2, color="orange", ls="-.", label=r"$x<0$, build2" )

ax.legend( loc="best" )
ax.set_xlabel( r"$y$" )
ax.set_ylabel( r"${\rm atan2}(x,y)$" )
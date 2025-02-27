'''
Plot_Gridpoint_Separation_Visualisation.py

Plot two arbitrary gridpoints on a 2D grid.
Do it in Cartesian coordinates but label it in polar coordinates (r,theta).

https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Arc.html
'''

import matplotlib.pyplot  as plt
import numpy              as np
import matplotlib.patches as mpatches # mpatches.Arc()



#----- Define variables -----
x_min = 0
x_max = 2.0
y_min = 0
y_max = 2.0

x_1 = 0.8
y_1 = 0.3

x_2 = 1.7
y_2 = 1.2



#----- Setup graphs -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":20 })


#----- Plot graph with points and joining line -----
fig = plt.figure( figsize=(7,5) )
ax  = fig.add_subplot( 111 )

#- Radius vector to points from origin -
ax.plot( [ 0, x_1 ], [ 0, y_1 ], "darkgrey", zorder=0 )
ax.plot( [ 0, x_2 ], [ 0, y_2 ], "darkgrey", zorder=0 )

#- The points themselves -
ax.scatter( x_1, y_1, color="k", s=60, zorder=2 )
ax.scatter( x_2, y_2, color="k", s=60, zorder=2 )

#- Separation btween the points -
ax.plot( [ x_1, x_2 ], [y_1, y_2 ], "k", ls="--", zorder=1 )

ax.set_xlim( x_min, x_max )
ax.set_ylim( y_min, y_max )


#--- Draw arcs representing angles ---

theta_1_degrees = np.arctan2( y_1, x_1 ) * 180.0/np.pi
theta_2_degrees = np.arctan2( y_2, x_2 ) * 180.0/np.pi

'''
First argument: x,y coordinates of origin of arc.
Second argument: width of arc (length of horizontal axis).
Third argument: height of arc (length of vertical axis).
Fourth argument angle: rotation of ellipse.
Fifth argument theta1: Starting angle in degrees.
Sixth argument theta2: Ending angle in degrees.

Start my arc at the origin.
width and height are equal. Using 2.0*x_1 would mean it interesected X_1 exactly if x_1 were on the x-axis.
(Because then it's just a diameter.) Set to a number higher than 2.0 that works for that datapoint.
Rotation is zero - too confusing otherwise.
Starts at the theta-coordinate of x_1 (from the x-axis anticlockwise) in degrees.
Ends at the y-axis.
But my code plots theta=0 at the North pole, so I want my arc to come from the y-axis, not the x-axis.

But we don't want the angle-arc to go to the point because it doesn't look good, especially for far-out points.
So just choose an arbitrary number inside.

'''
pac_1 = mpatches.Arc(
	[ 0, 0 ],
	0.5,
	0.5,
	angle  = 0,
	theta1 = theta_1_degrees,
	theta2 = 90,
	color="darkgrey",
	zorder=0
)

pac_2 = mpatches.Arc(
	[ 0, 0 ],
	0.9,
	0.9,
	angle  = 0,
	theta1 = theta_2_degrees,
	theta2 = 90,
	color="darkgrey",
	zorder=0
)

ax.add_patch( pac_1 )
ax.add_patch( pac_2 )


#----- Add coordinate labels to graph -----
ztext_r1     = r"$r_1$"
ztext_r2     = r"$r_2$"
ztext_theta1 = r"$\theta_1$"
ztext_theta2 = r"$\theta_2$"
ztext_delta  = r"$\delta(r_1,\theta_1,r_2,\theta_2)$"

plt.text(
	0.75 * x_1,
	0.75 * y_1 + 0.08,
	ztext_r1,
	ha       = "center",
	va       = "center",
	rotation = 0.8 * theta_1_degrees,
	size     = 14
)

plt.text(
	0.75 * x_2,
	0.75 * y_2 + 0.08,
	ztext_r2,
	ha       = "center",
	va       = "center",
	rotation = 0.8 * theta_2_degrees,
	size     = 14
)

plt.text(
	0.1,
	0.15,
	ztext_theta1,
	ha       = "center",
	va       = "center",
	rotation = 0,
	size     = 14
)

plt.text(
	0.25,
	0.3,
	ztext_theta2,
	ha       = "center",
	va       = "center",
	rotation = 0,
	size     = 14
)

plt.text(
	0.5  * ( x_1 + x_2 ),
	0.24 * ( x_1 + x_2 ),
	ztext_delta,
	ha       = "center",
	va       = "center",
	rotation = 0.6 * ( theta_1_degrees + theta_2_degrees ),
	size     = 14
)



#----- xticks and yticks -----
xticks_values = [ x_1, x_2 ]
xticks_labels = [ r"$x_1$", r"$x_2$" ]
yticks_values = [ y_1, y_2 ]
yticks_labels = [ f"$z_1$", r"$z_2$" ]
plt.xticks( xticks_values, xticks_labels )
plt.yticks( yticks_values, yticks_labels )

fig.tight_layout()
plt.savefig( "../Figures/Gridpoint_Separation_Visualisation" )
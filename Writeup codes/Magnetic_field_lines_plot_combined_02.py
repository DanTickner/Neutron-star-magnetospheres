'''
Magnetic_field_lines_plot_combined_02.py
Plot the combined CSV file for many magnetic field lines originating from various points on the surface of the object.

https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Circle.html
https://stackoverflow.com/questions/9215658/plot-a-circle-with-pyplot

remove default xmax and just allow user to make independent zmax and xmax
'''


import csv
import matplotlib.pyplot as plt
import numpy as np# Only for np.sqrt() to get golden ratio figsize
from matplotlib.patches import Circle

csv_filename = "CSV/20230831_Dipole_combined.csv"
fig_filename = "../Figures/Magnetic dipole field lines V02"

show_legend     = False
use_custom_xmax = True
use_custom_zmax = True

custom_xmax = 5
custom_zmax = 2

fig_height = 10


#----- Read CSV -----
labels = []

x_values = []
z_values = []

with open( csv_filename, "r" ) as csv_file:
	data = csv.reader( csv_file )
	
	headers = next( data )
	
	for row in data:
		
		if not row[0] in labels:
			labels.append( row[0] )
			x_values.append( [] )
			z_values.append( [] )
		
		n = labels.index( row[0] )
		x_values[n].append( float( row[3] ) )
		z_values[n].append( float( row[4] ) )



#----- Plot on graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":30 })

fig = plt.figure( figsize=( fig_height, fig_height/np.sqrt(2)) )
ax  = fig.add_subplot( 111 )

for n in range(len( labels )):
	plt.plot( x_values[n], z_values[n], label=labels[n], color="grey", ls="-", lw=0.5 )

ax.add_patch( Circle( (0,0), radius=1, fill=True, color="k"  ) )

plt.axis( "scaled" )

if use_custom_xmax:
	ax.set_xlim( 0, custom_xmax )
if use_custom_zmax:
	ax.set_ylim( -custom_zmax, custom_zmax )	
if show_legend:
	ax.legend( loc="best" )

ax.set_xlabel( r"$x/R_*$" )
ax.set_ylabel( r"$z/R_*$" )

fig.tight_layout()
plt.savefig( fig_filename, bbox_inches="tight" )
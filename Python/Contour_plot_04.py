''''
Contour_plot_04
Contour plot a scalar as a function of x and z, converting to Cartesian coordinates.

V04: Extract R_LC from _2_history.csv
'''

import matplotlib.pyplot as plt
import csv
import numpy as np

#----- Define variables -----
csv_filename_base  = "20240411_a_benchmark_divergenceless"

csv_filename_profiles = "CSV/" + csv_filename_base + "_1_profiles.csv"
csv_filename_history  = "CSV/" + csv_filename_base + "_2_history.csv"

function = "B_p"

# Value of T at which we wish to produce the plot
T_index_chosen = 4725

fig_filename = "Figures/20231222/Contour_" + csv_filename_base + "_" + function + "_" + str( T_index_chosen )

use_absolute_values                 = False
use_filled_contour                  = True
use_custom_contours                 = False
use_scientific_notation_in_colorbar = False

title_time_unit = "code units" # "seconds" or "code units"

#custom_contours = [ 0.03*x for x in range( 21 ) ]
#custom_contours = [ -20, -1e-3, 0, 1e-3, 0.3 ]
custom_contours = [ -1.1, -0.1, 0.1, 1.9 ]
custom_colours  = [ [                                                               "red" ],
					[                                                     "orange", "red" ],
					[                                           "yellow", "orange", "red" ],
					[                             "lightgreen", "yellow", "orange", "red" ],
					[                     "blue", "lightgreen", "yellow", "orange", "red" ],
					[           "indigo", "blue", "lightgreen", "yellow", "orange", "red" ],
					[ "violet", "indigo", "blue", "lightgreen", "yellow", "orange", "red" ] ]


T_seconds    = 0
T_codeunits  = 0
T_index_list = []	# Store all T indices to print so that the user can make decisions

i_max = 0
j_max = 0

i_list = []
j_list = []


#--- Automatic conditioning ---
if function == "is_gridpoint_within_FF_region":
	use_custom_contours                 = True
	use_scientific_notation_in_colorbar = False
	custom_contours = [ -0.1, 0.1, 1.1 ]
	
if function == "were_FF_conditions_applied":
	use_custom_contours                 = True
	use_scientific_notation_in_colorbar = False
	custom_contours = [ -1.1, -0.1, 0.1, 1.1 ]


#----- Read profiles CSV file first time: Get numbers of coordinates -----
with open( folder+"/"+csv_filename_profiles, "r" ) as csv_file:
	data = csv.reader( csv_file )
	
	#----- Determine the columns to use -----
	headers_profiles = next( data )
	
	for row in data:
		
		T_index = int( row[ headers_profiles.index( "T_index" ) ] )
		i       = int( row[ headers_profiles.index( "i"       ) ] )
		j       = int( row[ headers_profiles.index( "j"       ) ] )
		
		if T_index not in T_index_list:
			T_index_list.append( T_index )
		if i not in i_list:
			i_list.append( i )
		if j not in j_list:
			j_list.append( j )
		
		if T_index == T_index_chosen:
			T_seconds   = float( row[ headers_profiles.index( "T_sec" ) ] )
			T_codeunits = float( row[ headers_profiles.index( "T"     ) ] )

n_points_r = len( i_list )
n_points_t = len( j_list )

print( f"n_points_r = { n_points_r },\tn_points_t = { n_points_t }" )
print( f"Chosen T_index:\t{ T_index_chosen }" )
print( T_index_list )
print( f"Number of time values:\t{ len( T_index_list ) }" )

if T_index_chosen not in T_index_list:
	raise ValueError( "Chosen time index not in CSV file." )

x_grid = np.zeros( shape=(n_points_r,n_points_t) )
z_grid = np.zeros( shape=(n_points_r,n_points_t) )
f_grid = np.zeros( shape=(n_points_r,n_points_t) )


#----- Read profiles CSV file second time: Populate magnetic field arrays -----
with open( folder+"/"+csv_filename_profiles, "r" ) as csv_file:
	data = csv.reader( csv_file )
	next( data )
	for row in data:
		if int( row[ headers_profiles.index( "T_index" ) ] ) == T_index_chosen:
			
			i_csv = int  ( row[ headers_profiles.index( "i"      ) ] )
			j_csv = int  ( row[ headers_profiles.index( "j"      ) ] )
			x     = float( row[ headers_profiles.index( "x"      ) ] )
			z     = float( row[ headers_profiles.index( "z"      ) ] )
			f     = float( row[ headers_profiles.index( function ) ] )
			
			if use_absolute_values:
				f = abs( f )
			
			i = i_list.index( i_csv )
			j = j_list.index( j_csv )
			
			x_grid[i][j] = x
			z_grid[i][j] = z
			f_grid[i][j] = f



#----- Read history CSV file to obtain R_LC values -----
R_LC = 0
with open( folder+"/"+csv_filename_history, "r" ) as csv_file:
	data = csv.reader( csv_file )
	headers_history = next( data )
	for row in data:
		if int( row[ headers_history.index( "T_index" ) ] ) == T_index_chosen:
			R_LC = float( row[ headers_history.index( "R_LC" ) ] )
			print( f"\nR_LC:\t{ R_LC }" )

#----- Print extremal values to the screen for quick validation -----
x_min = min([ min(xi) for xi in x_grid ])
x_max = max([ max(xi) for xi in x_grid ])

z_min = min([ min(zi) for zi in z_grid ])
z_max = max([ max(zi) for zi in z_grid ])

f_min = min([ min(fi) for fi in f_grid ])
f_max = max([ max(fi) for fi in f_grid ])

print( f"\nmin( x ):\t{ x_min }" )
print( f"max( x ):\t{ x_max }" )

print( f"\nmin( z ):\t{ z_min }" )
print( f"max( z ):\t{ z_max }" )

print( f"\nmin( { function } ):\t{ f_min }" )
print( f"max( { function } ):\t{ f_max }" )

if use_custom_contours:
	print( "\nChosen contours:" )
	print( custom_contours )
	if custom_contours[0] > f_min:
		print( "WARNING: Minimum contour is larger than minimum field value; some gridpoints will be missing from the plot." )
	if custom_contours[-1] < f_max:
		print( "WARNING: Maximum contour is larger than maximum field value; some gridpoints will be missing from the plot." )		

#----- Plot graph of scalar field -----
title_time = 0

# Eventually it would be good to have the unit in the title. Don't need yet as I always use code units, and the title string would be too long for me.
if   title_time_unit == "seconds":
	title_time = T_seconds
elif title_time_unit == "code units":
	title_time = T_codeunits
else:
	print( f"title_time_unit { title_time_unit } not understood. Setting to zero." )


plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig, ax = plt.subplots(1,1)

if use_filled_contour:
	if use_custom_contours:
		cp = ax.contourf( x_grid, z_grid, f_grid, levels=custom_contours, colors=custom_colours[len(custom_contours)-2] )
	else:
		cp = ax.contourf( x_grid, z_grid, f_grid )
else:
	if use_custom_contours:
		cp = ax.contour( x_grid, z_grid, f_grid, levels=custom_contours, colors=custom_colours[len(custom_contours)-2] )
	else:
		cp = ax.contour( x_grid, z_grid, f_grid )

if use_absolute_values:
	# For the graph title.
	function = "abs(" + function + ")"

plt.axis( "scaled" )
ax.set_xlabel( r"$x$" )
ax.set_ylabel( r"$z$" )
ax.set_title( rf"{ function } at $T={ '{0:.3g}'.format(title_time) }$" )

if R_LC > 0:
	ax.axvline( R_LC, color="grey", ls="--" )

if use_scientific_notation_in_colorbar:
	fig.colorbar( cp, format="%.0e" )
else:
	fig.colorbar( cp )

plt.savefig( folder+"/"+fig_filename, bbox_inches="tight" )
print( f"\nSaved figure:\t{ folder+'/'+fig_filename }" )
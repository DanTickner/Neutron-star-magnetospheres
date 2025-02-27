'''
Plot_Profiles.py

Plot a 2D snapshot of a quantity in the "Profiles" CSV with radius and angle, at a certain timestep.

I thought about comparing this to a second snapshot, where you can either use
the same CSV again by setting csv_filename_2_base to the same string,
or two different evolutions at the same time, by setting chosen_T_index_2 to the same value.
But you can't really compare two different plots on a heatmap,
so I think it's easier to just run the same code back-to-back and produce two separate heatmaps.
For the same reason, I won't calculate an exact value with this code.
Option to output data to a CSV in the future, which would allow all of this, but I don't want this for now.

REASON ARCHIVED: The quadrilateral approximation of the twisted region could be easily turned into one with arcs. Keep this version for reference.
'''

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../../../../../Documents/NS_code_offline_results/CSV" #../CSV"
fig_folder = "../../../../../Documents/NS_code_offline_results/Figures/20250220" #"../Figures/20241017"
csv_filename_base = "20250220_a_twisted_axisymmetric"

chosen_header  = "B_r"
chosen_T_index = 150

#fig_filename = csv_filename_base + "_" + chosen_header + "_" + str( chosen_T_index )
fig_filename = "Twisted_AS_" + chosen_header + "_" + str( chosen_T_index )

use_absolute_values = True

use_filled_contour                  = True
use_custom_contours                 = True
use_scientific_notation_in_colorbar = True
plot_force_free_region              = True

custom_contours = [ 0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000 ]
custom_colours  = [ [                                                                      "darkred" ],
					[                                                                      "red", "darkred" ],
					[                                                            "orange", "red", "darkred" ],
					[                                                  "yellow", "orange", "red", "darkred" ],
					[                                    "lightgreen", "yellow", "orange", "red", "darkred" ],
					[                           "green", "lightgreen", "yellow", "orange", "red", "darkred" ],
					[                   "blue", "green", "lightgreen", "yellow", "orange", "red", "darkred" ],
					[         "violet", "blue", "green", "lightgreen", "yellow", "orange", "red", "darkred" ],
					[ "navy", "violet", "blue", "green", "lightgreen", "yellow", "orange", "red", "darkred" ] ]

use_custom_ylims = False
ylim_min = 6
ylim_max = 7

use_custom_xlims = True
xlim_min = 0
xlim_max = 1.6

list_T_index = []
list_T       = []
list_f       = []
list_i       = []
list_j       = []
list_x       = []
list_z       = []

chosen_T    = 0 # To be determined by the code.
chosen_R_LC = 0 # To be determined by the code.

highlight_twisted_region = True # Option to add grey dotted lines defining the boundary of the twisted region.
twist_r1 = 1								# Copy these values from the output log files.
twist_r2 = 10
twist_t1 = 0.523598775598299
twist_t2 = 2.61799387799149


#--- Automatic conditioning ---
if chosen_header in [ "force_free_condition_1", "force_free_condition_2", "force_free_conditions", "is_gridpoint_within_FF_region" ]:
	use_custom_contours                 = True
	use_scientific_notation_in_colorbar = False
	custom_contours = [ -0.1, 0.1, 1.1 ]
	#custom_colours  = [ "red", "green" ]
	
if chosen_header == "were_FF_conditions_applied":
	use_absolute_values                 = False
	use_custom_contours                 = True
	use_scientific_notation_in_colorbar = False
	custom_contours = [ -1.1, -0.1, 0.1, 1.1 ]
		
#----- Read CSV first time: Get time and coordinates -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_4_profiles.csv" )

print( f"Reading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	if not chosen_header in headers:
		print( "Headers:" )
		print( headers )
		raise ValueError( f"Chosen header { chosen_header } not in CSV file headers." )
	
	# Get correct columns. Unfortunate that index "i" is repeated.
	i_T_index = headers.index( "T_index"     )
	i_i       = headers.index( "i"           )
	i_j       = headers.index( "j"           )
	i_x       = headers.index( "x"           )
	i_z       = headers.index( "z"           )
	i_f       = headers.index( chosen_header )
	
	for row in data:
		
		T_index = int  ( row[ i_T_index ] )
		i       = int  ( row[ i_i       ] )
		j       = int  ( row[ i_j       ] )
		
		
		if not T_index in list_T_index:
			list_T_index.append( T_index )
			list_T.append( 0 ) # To be filled in later.
		
		if not i in list_i:
			list_i.append( i )
		
		if not j in list_j:
			list_j.append( j )


#----- Read Time values CSV -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_1_time_values.csv" )

print( f"Reading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns
	i_T_index = headers.index( "T_index" )
	i_T       = headers.index( "T"       )
	i_R_LC    = headers.index( "R_LC"    )
	
	for row in data:
		
		T_index = int  ( row[ i_T_index ] )
		T       = float( row[ i_T       ] )
		R_LC    = float( row[ i_R_LC    ] )
		
		if T_index in list_T_index:
			list_T[ list_T_index.index( T_index ) ] = T
		
		if T_index == chosen_T_index:
			chosen_T    = T
			chosen_R_LC = R_LC



#----- Print list of times and fail code if chosen T_index not present -----
print( "T_index:" )
print( list_T_index )
print( "\nT:" )
print( list_T )

if not chosen_T_index in list_T_index:
	raise ValueError( f"Chosen_T_index { chosen_T_index } not in Profiles CSV." )

print( f"\nChosen T_index and T:\t{ chosen_T_index }\t{ chosen_T }" )

#----- Now we can build data arrays of the correct size -----
# Get number of radial and angular gridpoints that are used in the file.
n_r_used = len( list_i )
n_t_used = len( list_j )

grid_x = np.zeros( shape=(n_r_used,n_t_used) )
grid_z = np.zeros( shape=(n_r_used,n_t_used) )
grid_f = np.zeros( shape=(n_r_used,n_t_used) )


#----- Read CSV second time: Populate grids -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_4_profiles.csv" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	if not chosen_header in headers:
		print( "Headers:" )
		print( headers )
		raise ValueError( f"Chosen header { chosen_header } not in CSV file headers." )
	
	# Get correct columns. Unfortunate that index "i" is repeated.
	i_T_index = headers.index( "T_index"     )
	i_i       = headers.index( "i"           )
	i_j       = headers.index( "j"           )
	i_x       = headers.index( "x"           )
	i_z       = headers.index( "z"           )
	i_f       = headers.index( chosen_header )
	
	for row in data:
		
		T_index = int  ( row[ i_T_index ] )
		i       = int  ( row[ i_i       ] )
		j       = int  ( row[ i_j       ] )
		x       = float( row[ i_x       ] )
		z       = float( row[ i_z       ] )
		f       = float( row[ i_f       ] )
		
		if use_absolute_values:
			f = abs( f )
		
		# This refers to the element of the list of i values etc.
		# Good enough for now but need better third index system.
		index_i = list_i.index( i )
		index_j = list_j.index( j )
		
		if T_index == chosen_T_index:
			grid_x[index_i][index_j] = x
			grid_z[index_i][index_j] = z
			grid_f[index_i][index_j] = f


#----- Print extremal values to the screen for quick validation -----
x_min = min([ min(xi) for xi in grid_x ])
x_max = max([ max(xi) for xi in grid_x ])

z_min = min([ min(zi) for zi in grid_z ])
z_max = max([ max(zi) for zi in grid_z ])

f_min = min([ min(fi) for fi in grid_f ])
f_max = max([ max(fi) for fi in grid_f ])

print( f"\nmin( x ):\t{ x_min }" )
print(   f"max( x ):\t{ x_max }" )

print( f"\nmin( z ):\t{ z_min }" )
print(   f"max( z ):\t{ z_max }" )

print( f"\nmin( { chosen_header } ):\t{ f_min }" )
print(   f"max( { chosen_header } ):\t{ f_max }" )

if use_custom_contours:
	print( "\nChosen contours:" )
	print( custom_contours )
	if custom_contours[0] > f_min:
		print( "WARNING: Minimum contour is larger than minimum field value; some gridpoints will be missing from the plot." )
	if custom_contours[-1] < f_max:
		print( "WARNING: Maximum contour is larger than maximum field value; some gridpoints will be missing from the plot." )



#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig, ax = plt.subplots(1,1)

if use_filled_contour:
	if use_custom_contours:
		cp = ax.contourf( grid_x, grid_z, grid_f, levels=custom_contours, colors=custom_colours[len(custom_contours)-2] )
	else:
		cp = ax.contourf( grid_x, grid_z, grid_f )
else:
	if use_custom_contours:
		print( "hiorhe" )
		cp = ax.contour( grid_x, grid_z, grid_f, levels=custom_contours, colors=custom_colours[len(custom_contours)-2] )
	else:
		cp = ax.contour( grid_x, grid_z, grid_f )
		
plt.axis( "scaled" )
ax.set_xlabel( r"$x$" )
ax.set_ylabel( r"$z$" )
ax.set_title( rf"{ chosen_header } at $T={ '{0:.3g}'.format(chosen_T) }$" )

if chosen_R_LC > 0:
	ax.axvline( chosen_R_LC, color="black", ls="--" )

if use_scientific_notation_in_colorbar:
	fig.colorbar( cp, format="%.0e" )
else:
	fig.colorbar( cp )

if highlight_twisted_region:
	
	twist_x1 = twist_r1 * np.sin( twist_t1 )
	twist_x2 = twist_r1 * np.sin( twist_t2 )
	twist_x3 = twist_r2 * np.sin( twist_t1 )
	twist_x4 = twist_r2 * np.sin( twist_t2 )
	
	twist_z1 = twist_r1 * np.cos( twist_t1 )
	twist_z2 = twist_r1 * np.cos( twist_t2 )
	twist_z3 = twist_r2 * np.cos( twist_t1 )
	twist_z4 = twist_r2 * np.cos( twist_t2 )
	
	#ax.plot( [ twist_x1, twist_x3, twist_x4, twist_x2, twist_x1 ], [ twist_z1, twist_z3, twist_z4, twist_z2, twist_z1 ], color="black", ls=":" )
	
	arc_twist_t_list = np.linspace( twist_t1, twist_t2, 1000 )
	arc_inner_x_list = [ twist_r1 * np.sin( t ) for t in arc_twist_t_list ]
	arc_inner_z_list = [ twist_r1 * np.cos( t ) for t in arc_twist_t_list ]
	arc_outer_x_list = [ twist_r2 * np.sin( t ) for t in arc_twist_t_list ]
	arc_outer_z_list = [ twist_r2 * np.cos( t ) for t in arc_twist_t_list ]
	
	ax.plot( arc_inner_x_list, arc_inner_z_list, color="black", ls=":" )
	ax.plot( arc_outer_x_list, arc_outer_z_list, color="black", ls=":" )
	#ax.plot( [ twist_x1, twist_x3 ], [ twist_z1, twist_z3 ], color="black", ls=":" )
	#ax.plot( [ twist_x4, twist_x2 ], [ twist_z4, twist_z2 ], color="black", ls=":" )
	ax.plot( [ arc_inner_x_list[0], arc_outer_x_list[0] ], [ arc_inner_z_list[0], arc_outer_z_list[0] ], color="black", ls=":" )

fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath, bbox_inches="tight" )
print( f"\nSaved figure:\t{ fig_filepath }" )
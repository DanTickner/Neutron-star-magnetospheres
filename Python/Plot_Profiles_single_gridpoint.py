'''
Plot_Profiles_single_gridpoint.py

Plot the time-evolution of a single gridpoints from the "Profiles" CSV.
'''

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20240603"
csv_filename_base = "20240603_d_sponge_in_terms_of_time_derivative"

chosen_header  = "B_p_dT_sponge_rel_diff"
chosen_i       = 99
chosen_j       = 0

fig_filename = csv_filename_base + "_" + chosen_header + "_" + str( chosen_i ) + "_" + str( chosen_j )

use_absolute_values     = True # Not yet implemented 20240603
plot_y_equals_zero_line = True


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
list_r       = []
list_j       = []
list_t       = []

chosen_r     = 0 # To be updated by the code.
chosen_t     = 0 # To be updated by the code.
		
#----- Read Profiles CSV -----
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
		f       = float( row[ i_f       ] )
		
		
		if not T_index in list_T_index:
			list_T_index.append( T_index )
			list_T.append( 0 ) # To be filled in later.
		
		if not i in list_i:
			list_i.append( i )
			list_r.append( 0 ) # To be filled in later.
		
		if not j in list_j:
			list_j.append( j )
			list_t.append( 0 ) # To be filled in later.
		
		if i == chosen_i and j == chosen_j:
			list_f.append( f )


#----- Read Time values CSV -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_1_time_values.csv" )

print( f"Reading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns
	i_T_index = headers.index( "T_index" )
	i_T       = headers.index( "T"       )
	
	for row in data:
		
		T_index = int  ( row[ i_T_index ] )
		T       = float( row[ i_T       ] )
		
		if T_index in list_T_index:
			list_T[ list_T_index.index( T_index ) ] = T


#----- Read gridpoints CSV -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_2_gridpoints.csv" )

print( f"Reading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns
	i_i_or_j = headers.index( "i_or_j" )
	i_r      = headers.index( "r" )
	i_t      = headers.index( "t" )
	
	for row in data:
		
		i_or_j = int  ( row[ i_i_or_j ] )
		
		r = row[ i_r ]
		if r != "":
			r = float( r )
		t = row[ i_t ]
		if t != "":
			t = float( t )
		
		if i_or_j in list_i:
			list_r[ list_i.index( i_or_j ) ] = r
		if i_or_j in list_j:
			list_t[ list_j.index( i_or_j ) ] = t

chosen_r = list_r[ list_i.index( chosen_i ) ]
chosen_t = list_t[ list_j.index( chosen_j ) ]

#----- Print lists of coordinates and times, and fail code if chosen T_index not present -----
print( "i:" )
print( list_i )
print( "j:" )
print( list_j )
print( "T_index:" )
print( list_T_index )
print( "\nT:" )
print( list_T )

if not chosen_i in list_i:
	raise ValueError( f"Chosen_i { chosen_i } not in Profiles CSV." )
if not chosen_j in list_j:
	raise ValueError( f"Chosen_j { chosen_j } not in Profiles CSV." )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( list_T, list_f )

ax.set_title( f"r[{ chosen_i }]={ round(chosen_r,3) },\tt[{ chosen_j }]={ round(chosen_t,3) }" )

ax.set_xlabel( r"$T$ (code units)" )

ax.set_ylabel( chosen_header )

if plot_y_equals_zero_line:
	ax.axhline( 0, color="grey", ls="--", lw=0.5 )

fig.tight_layout()

fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath )
print( f"\nSaved figure:\t{ fig_filepath }" )
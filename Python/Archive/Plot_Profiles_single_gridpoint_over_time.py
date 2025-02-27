'''
Plot_Profiles_single_gridpoint_over_time.py

Plot the evolution with time of a quantity evaluated at a single gridpoint.
Option to compare to a theoretical value.
Also read the history CSV to extract quantities such as Omega and R_LC.
'''

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20240412"
#csv_filename_base = "20240411_b_benchmark_normal"
csv_filename_base = "20240411_a_benchmark_divergenceless"

chosen_header = "E_r"
chosen_i      = 99
chosen_j      = 198

fig_filename = csv_filename_base + "_" + chosen_header + "_i_" + str( chosen_i ) + "_j_" + str( chosen_j )
fig_filename += "_exact_generic"

graph_type = "semilogy" # "normal" or "semilogy"

use_absolute_values = False

calculate_exact_value = True

use_custom_ylims = True
ylim_min = 1e-4
ylim_max = 0.2

use_custom_xlims = True
xlim_min = 0.1
xlim_max = 1.6

list_T_index = []
list_T       = []
list_f       = []
list_i       = []
list_r       = []
list_j       = []
list_t       = []
list_Omega   = []
list_B_r     = []   # B_r     as calculated by the code.
list_B_t     = []	# B_theta as calculated by the code.
list_f_exact = []


def exact_value( T_index ):
	i_T_index = list_T_index.index( T_index )
	
	if chosen_header == "E_r":
		return list_Omega[i_T_index] * chosen_r * np.sin( chosen_t ) * list_B_t[i_T_index]
		#return list_Omega[i_T_index] * chosen_r**(-2.0) * ( np.sin(chosen_t)**2 )
	if chosen_header == "E_t":
		return list_Omega[i_T_index] * chosen_r * np.sin( chosen_t ) * list_B_r[i_T_index] * -1.0
		#return list_Omega[i_T_index] * chosen_r**(-2.0) * -2.0 * np.cos(chosen_t)
	
	return 0


#----- Read CSV first time: Get time and coordinates -----

if graph_type == "semilogy":
	use_absolute_values = True

csv_filepath = os.path.join( csv_folder, csv_filename_base + "_1_profiles.csv" )

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
	i_T       = headers.index( "T"           )
	i_i       = headers.index( "i"           )
	i_r       = headers.index( "r"           )
	i_j       = headers.index( "j"           )
	i_t       = headers.index( "t"           )
	i_f       = headers.index( chosen_header )
	i_B_r     = headers.index( "B_r"         )
	i_B_t     = headers.index( "B_t"         )
	
	for row in data:
		
		T_index = int  ( row[ i_T_index ] )
		T       = float( row[ i_T       ] )
		i       = int  ( row[ i_i       ] )
		r       = float( row[ i_r       ] )
		j       = int  ( row[ i_j       ] )
		t       = float( row[ i_t       ] )
		f       = float( row[ i_f       ] )
		B_r     = float( row[ i_B_r     ] )
		B_t     = float( row[ i_B_t     ] )
		
		
		if not T_index in list_T_index:
			list_T_index.append( T_index )
			list_T      .append( T       )
		
		if not i in list_i:
			list_i.append( i )
			list_r.append( r )
		
		if not j in list_j:
			list_j.append( j )
			list_t.append( t )
		
		if i == chosen_i and j == chosen_j:
			if use_absolute_values:
				f = abs( f )
			list_f.append( f )
			list_B_r.append( B_r )
			list_B_t.append( B_t )



#----- Print list of times and fail code if chosen T_index not present -----
print( "T_index:" )
print( list_T_index )
print( "\nT:" )
print( list_T )
print( "\ni:" )
print( list_i )
print( "\nj:" )
print( list_j )

if not chosen_i in list_i:
	raise ValueError( f"Chosen_i { chosen_i } not in Profiles CSV." )
if not chosen_j in list_j:
	raise ValueError( f"Chosen_i { chosen_j } not in Profiles CSV." )

chosen_r = list_r[ list_i.index( chosen_j ) ]
chosen_t = list_t[ list_j.index( chosen_j ) ]

print( f"\nChosen i and r:\t{ chosen_i }\t{ chosen_r }" )
print( f"\nChosen j and t:\t{ chosen_j }\t{ chosen_t }" )


#----- Print extremal values to the screen for quick validation -----
T_index_min = min( list_T_index )
T_index_max = max( list_T_index )

T_min = min( list_T )
T_max = max( list_T )

f_min = min( list_f )
f_max = max( list_f )

print( f"\nmin( T_index ):\t{ T_index_min }" )
print(   f"max( T_index ):\t{ T_index_max }" )

print( f"\nmin( T ):\t{ T_min }" )
print(   f"max( T ):\t{ T_max }" )

print( f"\nmin( { chosen_header } ):\t{ f_min }" )
print(   f"max( { chosen_header } ):\t{ f_max }" )


#----- Read history CSV file to obtain R_LC values -----

csv_filepath = os.path.join( csv_folder, csv_filename_base + "_2_history.csv" )

print( f"\nReading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers_history = next( data )
	i_T_index_history = headers_history.index( "T_index" )
	i_Omega_history   = headers_history.index( "Omega"   )
	
	for row in data:
		
		T_index = int( row[ i_T_index_history ] )
		
		if T_index in list_T_index:
			Omega = float( row[ i_Omega_history ] )
			list_Omega.append( Omega )


#----- Calculate exact values -----
if calculate_exact_value:
	for T_index in list_T_index:
		f_exact = exact_value( T_index )
		if use_absolute_values:
			f_exact = abs( f_exact )
		list_f_exact.append( f_exact )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if graph_type == "normal":
	ax.plot( list_T, list_f, color="k", ls="-", label="Calculated" )
elif graph_type == "semilogy":
	ax.semilogy( list_T, list_f, color="k", ls="-" , label="Calculated" )
else:
	raise ValueError( f"graph_type = { graph_type } not understood. Please choose either 'normal' or 'semilogy'." )

if calculate_exact_value:
	if graph_type == "normal":
		ax.plot( list_T, list_f_exact, color="r", ls="--", label="Exact"      )
	elif graph_type == "semilogy":
		ax.semilogy( list_T, list_f_exact, color="r", ls="--", label="Exact"      )	

ax.set_xlabel( "Time (code units)" )

if chosen_header == "E_r":
	ax.set_ylabel( "E_r (code units)" )
if chosen_header == "E_t":
	ax.set_ylabel( "E_t (code units)" )

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

ax.legend( loc="best" )

fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath, bbox_inches="tight" )
print( f"\nSaved figure:\t{ fig_filepath }" )
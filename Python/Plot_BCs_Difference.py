'''
Plot_BCs_Difference.py

Plot the difference in boundary values (before and after application of boundary conditions) as a function of time.
Ability to plot for multiple evolutions on the same graph.
'''

import os
import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20240527"

csv_filename_bases = [ 
	"20240621_a_test"
	]

graph_labels = [ r"Old", r"Fixed outer BCs", "Half nt", "Half ellmax" ]

chosen_header = "E_inner_BCs_difference_p"

chosen_j = 30

fig_filename = "20240527_" + chosen_header + "_j_" + str( chosen_j )

graph_type = "normal" # "normal" or "semilogy"

set_title  = False
graph_title = r"$n_r=n_\theta=1000$ and $\ell_{\rm{max}}=20$; no $\phi$-evolution"

calculate_exact_value = False

use_absolute_values = False

use_custom_ylims = False
ylim_min = 0
ylim_max = 2

use_custom_xlims = False
xlim_min = 0
xlim_max = 0.11

colors = [ "black", "blue", "orange", "darkgreen", "grey", "yellow" ]

def exact_value( t ):
	
	# For now, enter these manually. Later, add code to read the log file and extract them.
	r_min = 1.0
	r_max = 10.024169921875
	r     = 1.367431640625
	theta = 0.261799387799149
	
	if chosen_header == "total_electric_energy":
		return ( r_min**-3 - r_max**-3 ) * Omega**2 * 32.0 * math.pi / 45.0
	
	if chosen_header == "E_p":
		return 0
		
		
	
	return 0

list_T_index_all = [ [] for lst in csv_filename_bases ]
list_T_all       = [ [] for lst in csv_filename_bases ] # BCs CSV only contains a sample of values. These might be different for different runs, so keep as n seprate lists.
list_T_index     = [ [] for lst in csv_filename_bases ]
list_T           = [ [] for lst in csv_filename_bases ]
list_f_calc      = [ [] for lst in csv_filename_bases ]


#----- Read Time values CSV -----
for n, filename_base in enumerate( csv_filename_bases ):
	
	csv_filepath = os.path.join( csv_folder, filename_base + "_1_time_values.csv" )
	
	print( f"Reading CSV file:\t{ csv_filepath }" )
	
	with open( csv_filepath, "r" ) as input_file:
		
		data = csv.reader( input_file )
		
		headers = next( data )
		
		# Get correct columns
		i_T_index = headers.index( "T_index" )
		i_T       = headers.index( "T"       )
		
		for row in data:
			
			list_T_index_all[n].append( int  ( row[ i_T_index ] ) )
			list_T_all      [n].append( float( row[ i_T       ] ) )



#----- Read BCs CSV -----
for n, filename_base in enumerate( csv_filename_bases ):
	
	csv_filepath = os.path.join( csv_folder, filename_base + "_7_BCs.csv" )
	
	print( f"Reading CSV file:\t{ csv_filepath }" )
	
	with open( csv_filepath, "r" ) as input_file:
		
		data = csv.reader( input_file )
		
		headers = next( data )
		
		if not chosen_header in headers:
			print( "Headers:" )
			print( headers )
			raise ValueError( f"Chosen header { chosen_header } not in CSV file headers." )
		
		# Get correct columns
		i_T_index = headers.index( "T_index"     )
		i_f_calc  = headers.index( chosen_header )
		i_j       = headers.index( "j"           ) # Unfortunate that i and j almost clash here, but consistent with the convention used in the rest of the code.
		
		for row in data:
			f_calc  = float( row[ i_f_calc  ] )
			T_index = int  ( row[ i_T_index ] )
			T       = list_T_all[n][ list_T_index_all[n].index( T_index ) ]
			if int( row[ i_j ] ) == chosen_j:
				list_f_calc [n].append( f_calc  )
				list_T_index[n].append( T_index )
				list_T      [n].append( T       )


if use_absolute_values:
	for n in range(len( list_f_calc ) ):
		for i in range(len( list_f_calc[n] ) ):
			list_f_calc[n][i] = abs( list_f_calc[n][i] )

#----- Calculate exact values -----
# To do 20240412: Combine list_T_1 and list_T_2 and sort from low to high, so that all values are included.
if calculate_exact_value:
	
	'''
	#--- Guarantee that all time values are included ---
	for n in range(len( list_T )):
		for i in range(len( list_T[n] )):
			if not list_T[n][i] in list_T_all:
				list_T_all.append( list_T[n][i] )
	
	list_T_all.sort()
	'''
	# Don't use all T_values because the Omega might be different for different runs. Instead, just use for the first run.
	
	#--- Calculate expected values for these times ---
	#for T in list_T_all:
	for T in list_T[0]:
		if use_absolute_values:
			list_f_exact.append( abs( exact_value( T ) ) )
		else:
			list_f_exact.append( exact_value( T ) )



#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if graph_type == "normal":
	for n in range(len( csv_filename_bases ) ):
		ax.plot( list_T[n], list_f_calc[n] , color=colors[n], ls="-" , label=graph_labels[n] )
elif graph_type == "semilogy":
	for n in range(len( csv_filename_bases ) ):
		ax.semilogy( list_T[n], list_f_calc[n] , color=colors[n], ls="-" , label=graph_labels[n] )
else:
	raise ValueError( f"graph_type = { graph_type } not understood. Please choose either 'normal' or 'semilogy'." )

if calculate_exact_value:
	if graph_type == "normal":
		ax.plot( list_T[0], list_f_exact, color="r", ls="--", label="Exact"      )
	elif graph_type == "semilogy":
		ax.semilogy( list_T[0], list_f_exact, color="r", ls="--", label="Exact"      )

ax.set_xlabel( "Time (code units)" )
if use_absolute_values:
	ax.set_ylabel( f"abs( {chosen_header} )" )
else:
	ax.set_ylabel( chosen_header )

if   chosen_header == "total_electric_energy":
	ax.set_ylabel( "total electric energy (code units)" )
elif chosen_header == "total_magnetic_energy":
	ax.set_ylabel( "total magnetic energy (code units)" )

if set_title:
	ax.set_title( graph_title )

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

if len( csv_filename_bases ) > 1 or calculate_exact_value:
	ax.legend( loc="best" )

fig.tight_layout()
fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath )
print( f"Saved figure:\t{ fig_filepath }"  )
'''
Plot_History.py

Plot the evolution of a quantity in the "History" CSV file with time.
If the quantity has a known exact value, option to calculate that and plot alongside.

New version allowing you to plot from multiple CSVs for direct comparison.

REASON ARCHIVED: CSV system changed so that time and coordinate values are output to different files, so now need to read two CSVs.
'''

import os
import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20240510"

csv_filename_bases = [ "20240510_a_fewer_gridpoints_ellmax_10",
					   "20240510_b_fewer_gridpoints_ellmax_15",
					   "20240510_c_fewer_gridpoints_ellmax_20",
					   "20240510_d_fewer_gridpoints_ellmax_25",
					   "20240510_e_fewer_gridpoints_ellmax_30" ]

graph_labels = [  "10", "15", "20", "25", "30" ]

chosen_header = "stdev_B_r"

fig_filename = "20240510_compare_" + chosen_header

graph_type = "semilogy" # "normal" or "semilogy"

set_title  = True
graph_title = r"Varying $\ell_{rm{max}}$"

calculate_exact_value = False

use_custom_ylims = False
ylim_min = 4
ylim_max = 10

use_custom_xlims = False
xlim_min = 0
xlim_max = 1.6

colors = [ "black", "blue", "orange", "darkgreen", "grey" ]

def exact_value( t ):
	
	if chosen_header == "total_magnetic_energy":
		r_min = 1.0
		r_max = 10.024169921875
		return ( r_min**-3 - r_max**-3 ) * 4.0* math.pi / 3.0
	
	return 0

list_T_index   = [ [] for lst in csv_filename_bases ]
list_T         = [ [] for lst in csv_filename_bases ]
list_f_calc    = [ [] for lst in csv_filename_bases ]
list_f_exact   = []


#----- Read CSVs -----
for n, filename_base in enumerate( csv_filename_bases ):
	
	csv_filepath = os.path.join( csv_folder, filename_base + "_5_history.csv" )
	
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
		i_T       = headers.index( "T"           )
		i_f_calc  = headers.index( chosen_header )
		
		for row in data:
			
			list_T_index[n].append( int  ( row[ i_T_index ] ) )
			list_T      [n].append( float( row[ i_T       ] ) )
			list_f_calc [n].append( float( row[ i_f_calc  ] ) )



#----- Calculate exact values -----
# To do 20240412: Combine list_T_1 and list_T_2 and sort from low to high, so that all values are included.
if calculate_exact_value:
	for T in list_T[0]:
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

ax.set_title( chosen_header.replace( "_", " " ) )

ax.set_xlabel( "Time (code units)" )

if chosen_header == "total_magnetic_energy":
	ax.set_ylabel( "total magnetc energy (code units)" )

if set_title:
	ax.set_title( graph_title )

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

if len( csv_filename_bases ) > 1:
	ax.legend( loc="best" )

fig.tight_layout()
fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath )
print( f"Saved figure:\t{ fig_filepath }"  )
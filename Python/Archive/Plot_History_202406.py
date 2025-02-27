'''
Plot_History.py

Plot the evolution of a quantity in the "History" CSV file with time.
If the quantity has a known exact value, option to calculate that and plot alongside.

New version allowing you to plot from two CSVs for direct comparison.
To use one, set compare_two_files = False.
'''

import os
import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20241016"
csv_filename_1_base = "20241014_a_test"
csv_filename_2_base = "20240411_b_benchmark_normal"

chosen_header = "stdev_B_r"

fig_filename = "20241014_a_test_" + chosen_header

graph_type = "semilogy" # "normal" or "semilogy"

graph_label_1 = "divergenceless"
graph_label_2 = "normal"

calculate_exact_value = False
compare_two_files     = False

use_custom_ylims = False
ylim_min = 1e-7
ylim_max = 1e-1

use_custom_xlims = False
xlim_min = 0
xlim_max = 1.6

def exact_value( t ):
	
	if chosen_header == "total_magnetic_energy":
		r_min = 1.0
		r_max = 10.024169921875
		return ( r_min**-3 - r_max**-3 ) * 4.0* math.pi / 3.0
	
	return 0

list_T_index_1 = []
list_T_index_2 = []
list_T_1       = []
list_T_2       = []
list_f_calc_1  = []
list_f_calc_2  = []
list_f_exact   = []


#----- Read first CSV -----
csv_filepath_1 = os.path.join( csv_folder, csv_filename_1_base + "_2_history.csv" )

print( f"Reading CSV file:\t{ csv_filepath_1 }" )

with open( csv_filepath_1, "r" ) as input_file_1:
	
	data = csv.reader( input_file_1 )
	
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
		
		list_T_index_1.append( int  ( row[ i_T_index ] ) )
		list_T_1      .append( float( row[ i_T       ] ) )
		list_f_calc_1 .append( float( row[ i_f_calc  ] ) )



#----- Read second CSV -----
if compare_two_files:
	csv_filepath_2 = os.path.join( csv_folder, csv_filename_2_base + "_2_history.csv" )
	
	print( f"Reading CSV file:\t{ csv_filepath_2 }" )
	
	with open( csv_filepath_2, "r" ) as input_file_2:
		
		data = csv.reader( input_file_2 )
		
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
			
			list_T_index_2.append( int  ( row[ i_T_index ] ) )
			list_T_2      .append( float( row[ i_T       ] ) )
			list_f_calc_2 .append( float( row[ i_f_calc  ] ) )



#----- Calculate exact values -----
# To do 20240412: Combine list_T_1 and list_T_2 and sort from low to high, so that all values are included.
if calculate_exact_value:
	for T in list_T_1:
		list_f_exact.append( exact_value( T_1 ) )



#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

if graph_type == "normal":
	ax.plot( list_T_1, list_f_calc_1 , color="k", ls="-" , label=graph_label_1 )
elif graph_type == "semilogy":
	ax.semilogy( list_T_1, list_f_calc_1 , color="k", ls="-" , label=graph_label_1 )
else:
	raise ValueError( f"graph_type = { graph_type } not understood. Please choose either 'normal' or 'semilogy'." )

if compare_two_files:
	if graph_type == "normal":
		ax.plot( list_T_2, list_f_calc_2 , color="b", ls="-" , label=graph_label_2 )
	elif graph_type == "semilogy":
		ax.semilogy( list_T_2, list_f_calc_2 , color="b", ls="-" , label=graph_label_2 )

if calculate_exact_value:
	if graph_type == "normal":
		ax.plot( list_T_1, list_f_exact, color="r", ls="--", label="Exact"      )
	elif graph_type == "semilogy":
		ax.semilogy( list_T_1, list_f_exact, color="r", ls="--", label="Exact"      )		

ax.set_title( chosen_header.replace( "_", " " ) )

ax.set_xlabel( "Time (code units)" )

if chosen_header == "total_magnetic_energy":
	ax.set_ylabel( "total magnetc energy (code units)" )

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

ax.legend( loc="best" )

fig.tight_layout()
fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath )
print( f"Saved figure:\t{ fig_filepath }"  )
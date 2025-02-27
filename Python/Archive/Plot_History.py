'''
Plot_History.py

Plot the evolution of a quantity in the "History" CSV file with time.
If the quantity has a known exact value, option to calculate that and plot alongside.
'''

import os
import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../CSV"
fig_folder = "../Figures/20240412"
csv_filename_base = "20240411_a_benchmark_divergenceless"

chosen_header = "total_magnetic_energy"

calculate_exact_value = False

use_custom_ylims = True
ylim_min = 6
ylim_max = 7

use_custom_xlims = True
xlim_min = 0
xlim_max = 1.6

def exact_value( t ):
	
	if chosen_header == "total_magnetic_energy":
		r_min = 1.0
		r_max = 10.024169921875
		return ( r_min**-3 - r_max**-3 ) * 4.0* math.pi / 3.0
	
	return 0

list_T_index     = []
list_T           = []
list_T_SI        = []
list_value_calc  = []
list_value_exact = []


#----- Read CSV -----
csv_filepath = os.path.join( csv_folder, csv_filename_base + "_2_history.csv" )

print( f"Reading CSV file:\t{ csv_filepath }" )

with open( csv_filepath, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	if not chosen_header in headers:
		print( "Headers:" )
		print( headers )
		raise ValueError( f"Chosen header { chosen_coeff } not in CSV file headers." )
	
	# Get correct columns
	i_T_index    = headers.index( "T_index"     )
	i_T          = headers.index( "T"           )
	i_T_SI       = headers.index( "T_SI"        )
	i_value_calc = headers.index( chosen_header )
	
	for row in data:
		
		list_T_index   .append( int  ( row[ i_T_index ] )    )
		list_T         .append( float( row[ i_T       ] )    )
		list_T_SI      .append( float( row[ i_T_SI    ] )    )
		list_value_calc.append( float( row[ i_value_calc ] ) )



#----- Calculate exact values -----
if calculate_exact_value:
	for T in list_T:
		list_value_exact.append( exact_value( T ) )



#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( list_T, list_value_calc , color="k", ls="-" , label="Calculated" )

if calculate_exact_value:
	ax.plot( list_T, list_value_exact, color="b", ls="--", label="Exact"      )
	ax.legend( loc="best" )

ax.set_title( chosen_header.replace( "_", " " ) )

ax.set_xlabel( "Time (code units)" )

if chosen_header == "total_magnetic_energy":
	ax.set_ylabel( "total magnetc energy (code units)" )

ax.set_xlim( 0, list_T[-2]*1.05 )	# 20240323 The final T value is zero.

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

fig.tight_layout()
fig_filepath = os.path.join( fig_folder, csv_filename_base + "_" + chosen_header )
plt.savefig( fig_filepath )
print( f"Saved figure:\t{ fig_filepath }"  )
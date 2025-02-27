'''
Plot_History_Multiple.py

Plot the evolution of a quantity in the "History" CSV file with time.
If the quantity has a known exact value, option to calculate that and plot alongside.

New version allowing you to plot from multiple CSVs for direct comparison.
'''

import os
import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_folder = "../../../../../Documents/NS_code_offline_results/CSV" # "../CSV"
fig_folder = "../../../../../Documents/NS_code_offline_results/Figures/20250220" #"../Figures/20241017"

csv_filename_bases = [
	"20250205_a_detailed_run_nonrotating",
	#"20250205_b_detailed_run_nonrotating_ellmax1"
	#"20250205_c_detailed_run_rotating",
	"20250219_a_twisted",
	"20250220_a_twisted_axisymmetric"
]

graph_labels = [
	r"Nontwisted",
	r"Twisted (NH)",
	r"Twisted (Symm.)"
]

chosen_header = "E_p"

#fig_filename = csv_filename_bases[0] + "_" + chosen_header
#fig_filename = "Twisted_AS_" + chosen_header
fig_filename = "Twisted_Both_" + chosen_header

graph_type = "normal" # "normal" or "semilogy"

set_title  = False
graph_title = r"Which B-values to use during rotation ramp?"
graph_title = csv_filename_bases[0]

calculate_exact_value = True

use_absolute_values = False

use_custom_ylims = True
ylim_min = -0.06
ylim_max = 0.01

use_custom_xlims = True
xlim_min = 0
xlim_max = 9

add_axvlines = True
axvlines = [
	0.457763672
]

use_custom_legend_loc = False
custom_legend_loc = "upper right"

colors = [ "black", "blue", "orange", "darkgreen", "grey", "yellow" ]

def exact_value( t ):
	
	# For now, enter these manually. Later, add code to read the log file and extract them.
	r_min = 1.0
	r_max = 10.0020751953125
	r     = 4.67431640625
	theta = 0.644026493985908
	Omega = list_Omega[0][ list_T[0].index( t ) ] # The value at time t.
	
	total_electric_energy = ( r_max - r_min ) * Omega**2 * 32.0 * math.pi / 15.0
	total_magnetic_energy = ( r_min**-3 - r_max**-3 ) * 4.0* math.pi / 3.0
	
	if chosen_header == "total_electric_energy":
		return total_electric_energy
	
	if chosen_header == "total_magnetic_energy":
		return total_magnetic_energy
	
	if chosen_header == "total_energy":
		return total_electric_energy + total_magnetic_energy
	
	if chosen_header == "B_r":
		return 2.0 * math.cos(theta) / ( r**3 )
	
	if chosen_header == "B_t":
		return math.sin(theta) / ( r**3 )
	
	if chosen_header == "E_r":
		 return Omega * ( math.sin(theta)**2 ) / ( r**2 )
	
	if chosen_header == "E_t":
		return - 2.0 * Omega * math.sin(theta) * math.cos(theta) / ( r**2 )
	
	if chosen_header == "E_p":
		return 0
		
		
	
	return 0

list_T_index   = [ [] for lst in csv_filename_bases ]
list_T         = [ [] for lst in csv_filename_bases ]
list_Omega     = [ [] for lst in csv_filename_bases ] # Needed for some exact calculations.
list_f_calc    = [ [] for lst in csv_filename_bases ]
list_f_exact   = []
list_T_all     = []
list_Omega_all = []


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
		i_Omega   = headers.index( "Omega"   )
		
		for row in data:
			
			list_T_index[n].append( int  ( row[ i_T_index ] ) )
			list_T      [n].append( float( row[ i_T       ] ) )
			list_Omega  [n].append( float( row[ i_Omega   ] ) )
			list_f_calc [n].append( 0 ) # Will be added to when reading the next CSV.



#----- Read history CSV -----
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
		i_f_calc  = headers.index( chosen_header )
		
		for row in data:
			list_f_calc[n][ list_T_index[n].index( int( row[ i_T_index ] ) ) ] = float( row[ i_f_calc  ] )


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
if len( csv_filename_bases ) == 1:
	graph_labels[0] = "Numerical"

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
elif chosen_header == "total_energy":
	ax.set_ylabel( "total energy (code units)" )

if set_title:
	ax.set_title( graph_title )

if use_custom_ylims:
	ax.set_ylim( ylim_min, ylim_max )
if use_custom_xlims:
	ax.set_xlim( xlim_min, xlim_max )

if add_axvlines:
	for axvline in axvlines:
		plt.axvline( axvline, color="grey", ls="--", lw=0.5 )

if len( csv_filename_bases ) > 1 or calculate_exact_value:
	if use_custom_legend_loc:
		ax.legend( loc=custom_legend_loc )
	else:
		ax.legend( loc="best" )

fig.tight_layout()
fig_filepath = os.path.join( fig_folder, fig_filename )
plt.savefig( fig_filepath )
print( f"Saved figure:\t{ fig_filepath }"  )
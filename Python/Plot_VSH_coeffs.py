'''
Plot_VSH_coeffs.py

Plot the run of a VSH coefficient with r, for a given ell and time index, for either B or E.
Option to calculate expected value and plot alongside.
'''


import csv
import matplotlib.pyplot as plt
import math


#----- Define variables -----
csv_filename_base = "../CSV/20240320_c_check_radial_coordinate"

chosen_coeff = "B_r_L"

chosen_ell     = 1
chosen_T_index = 100

calculate_exact_coeff = True

def exact_coeff_value( r ):
	lamda_inverse = ( 2.0 * r - r_max - r_min ) / ( r_max - r_min )
	if chosen_coeff == "B_r_L" and chosen_ell == 1:
		return 4.0 * math.sqrt( math.pi / 3.0 ) * pow( lamda_inverse, -3 )
	return 0

list_T_index     = []
list_i           = []
list_r           = []
list_coeff_calc  = []
list_coeff_exact = []

n_points_r = 0
r_min      = 999
r_max      = 0


#----- Read CSV -----
csv_filename = csv_filename_base + "_3_VSH_coeffs.csv"

print( f"Reading CSV file:\t{ csv_filename }" )

with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	if not chosen_coeff in headers:
		print( "Headers:" )
		print( headers )
		raise ValueError( f"Chosen coefficient { chosen_coeff } not in CSV file headers." )
	
	# Get correct columns
	i_T_index    = headers.index( "T_index"             )
	i_i          = headers.index( "i" )
	i_r          = headers.index( "r" )
	i_ell        = headers.index( "ell" )
	i_coeff_calc = headers.index( chosen_coeff )
	
	for row in data:
		
		T_index    = int  ( row[ i_T_index ] )
		i          = int  ( row[ i_i ] )
		r          = float( row[ i_r ] )
		coeff_calc = float( row[ i_coeff_calc ] )
		
		
		if not T_index in list_T_index:
			list_T_index.append( T_index )
		
		if i+1 > n_points_r:
			n_points_r = i+1
		if r < r_min:
			r_min = r
		if r > r_max:
			r_max = r
		
		if T_index == chosen_T_index and int( row[ i_ell ] ) == chosen_ell:
			list_i         .append( i )
			list_r         .append( r )
			list_coeff_calc.append( coeff_calc )


if not chosen_T_index in list_T_index:
	print( "T_index list:" )
	print( list_T_index )
	raise ValueError( f"Chosen T_index { chosen_T_index } not in CSV file." )


#----- Calculate exact values -----
if calculate_exact_coeff:
	for r in list_r:
		list_coeff_exact.append( exact_coeff_value( r ) )



#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( list_r, list_coeff_calc , color="k", ls="-" , label="Calculated" )
if calculate_exact_coeff:
	ax.plot( list_r, list_coeff_exact, color="b", ls="--", label="Exact"      )
	ax.legend( loc="best" )
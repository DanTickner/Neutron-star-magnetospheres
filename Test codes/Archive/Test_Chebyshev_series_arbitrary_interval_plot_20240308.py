'''
Test_Chebyshev_series_arbitrary_interval_plot.py

Plot the calculated Chebyshev series of a function and its exact values on an arbitrary interval,
using results already outputted to a CSV by the C++ code
Test_Chebyshev_series_arbitrary_interval.cpp
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename      = "../CSV/20240308_Test_Chebyshev_series_arbitrary_interval_function_values.csv"
fig_filename_base = "../Figures/20240308/20240308_Test_Chebyshev_series_arbitrary_interval"

i_list           = []	# Number of values not known a priori. Not strictly necessary as can use range(len(stdev_r)) but an extra safety net in case i not ordered.
r_list           = []
f_exact_list     = []
f_series_list    = []
dfdx_exact_list  = []
dfdx_series_list = []




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })



#----- Read CSV for stdev as a function of r -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
	# but it's better than obtaining the wrong column without checking.
	i_i           = headers.index( "i"           )
	i_r           = headers.index( "r"           )
	i_f_exact     = headers.index( "f_exact"     )
	i_f_series    = headers.index( "f_series"    )
	i_dfdx_exact  = headers.index( "dfdx_exact"  )
	i_dfdx_series = headers.index( "dfdx_series" )
	
	for row in data:
		i_list          .append( int  ( row[ i_i           ] ) )
		r_list          .append( float( row[ i_r           ] ) )
		f_exact_list    .append( float( row[ i_f_exact     ] ) )
		f_series_list   .append( float( row[ i_f_series    ] ) )
		dfdx_exact_list .append( float( row[ i_dfdx_exact  ] ) )
		dfdx_series_list.append( float( row[ i_dfdx_series ] ) )



#----- Determine y-limits for graphs -----
f_min = min( [ min( f_exact_list ), min( f_series_list ) ] ) * 0.95
f_max = max( [ max( f_exact_list ), max( f_series_list ) ] ) * 1.05

dfdx_min = min( [ min( dfdx_exact_list ), min( dfdx_series_list ) ] ) * 0.95
dfdx_max = max( [ max( dfdx_exact_list ), max( dfdx_series_list ) ] ) * 1.05



'''#----- Plot exact and series values of function -----
fig_f = plt.figure( figsize=(8,5) )
ax_f  = fig_f.add_subplot( 111 )

ax_f.plot( r_list, f_exact_list , color="k", ls="-" , label="Exact"  )
ax_f.plot( r_list, f_series_list, color="b", ls="--", label="Series" )

ax_f.legend( loc="best" )
ax_f.set_ylabel( "Function value $f(r)$" )
ax_f.set_xlabel( r"Radial coordinate $r$" )

ax_f.set_xlim( min(r_list), max(r_list) )
ax_f.set_ylim( f_min, f_max )

fig_f.tight_layout()
plt.savefig( fig_filename_base + "_f" )
'''


#----- Plot exact and series values of first derivative -----
fig_dfdx = plt.figure( figsize=(8,5) )
ax_dfdx  = fig_dfdx.add_subplot( 111 )

ax_dfdx.plot( r_list, dfdx_exact_list , color="k", ls="-" , label="Exact"  )
ax_dfdx.plot( r_list, dfdx_series_list, color="b", ls="--", label="Series" )

ax_dfdx.legend( loc="best" )
ax_dfdx.set_ylabel( "First derivative of function $df/dr$" )
ax_dfdx.set_xlabel( r"Radial coordinate $r$" )

ax_dfdx.set_xlim( min(r_list), max(r_list) )
ax_dfdx.set_ylim( dfdx_min, dfdx_max )

fig_dfdx.tight_layout()
plt.savefig( fig_filename_base + "_dfdx" )
'''
Test_radial_derivatives_FD_plot.py

Plot the calculated and exact first and second radial derivatives of a function,
where the calculation method is by finite differencing.
'''

import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename = "../CSV/20240411_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3_function_values.csv"

fig_filename_base = "../Figures/20240411/20240411_Test_radial_derivatives_FD_five_points"
fig_filename_dr  = fig_filename_base + "_dr"
fig_filename_dr2 = fig_filename_base + "_dr2"

i_list                 = []
f_dr_exact_list        = []
f_dr_calc_list         = []
f_dr_abs_rel_err_list  = []	# Absolute relative error
f_dr2_exact_list       = []
f_dr2_calc_list        = []
f_dr2_abs_rel_err_list = []	# Absolute relative error

ylims_dr = [ 0, 10 ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Read CSV files for stdev as a function of r -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
	# but it's better than obtaining the wrong column without checking.
	i_i                 = headers.index( "i"                 )
	i_f_dr_exact        = headers.index( "f_dr_exact"        )
	i_f_dr_calc         = headers.index( "f_dr_FD"           )
	i_f_dr_abs_rel_err  = headers.index( "abs_rel_err_f_dr"  )
	i_f_dr2_exact       = headers.index( "f_dr2_exact"       )
	i_f_dr2_calc        = headers.index( "f_dr2_FD"          )
	i_f_dr2_abs_rel_err = headers.index( "abs_rel_err_f_dr2" )
	
	for row in data:
		i_list                .append( int  ( row[ i_i                 ] ) )
		f_dr_exact_list       .append( float( row[ i_f_dr_exact        ] ) )
		f_dr_calc_list        .append( float( row[ i_f_dr_calc         ] ) )
		f_dr_abs_rel_err_list .append( float( row[ i_f_dr_abs_rel_err  ] ) )
		f_dr2_exact_list      .append( float( row[ i_f_dr2_exact       ] ) )
		f_dr2_calc_list       .append( float( row[ i_f_dr2_calc        ] ) )
		f_dr2_abs_rel_err_list.append( float( row[ i_f_dr2_abs_rel_err ] ) )




#----- Plot abs rel err of first derivative as a function of r -----
fig_dr = plt.figure( figsize=(8,5) )
ax_dr  = fig_dr.add_subplot( 111 )

ax_dr.semilogy( i_list, f_dr_abs_rel_err_list, color="k", ls="-" )

ax_dr.set_ylabel( "Absolute relative error" )
ax_dr.set_xlabel( r"$r$-index $i$" )

ax_dr.set_xlim( min(i_list)-1, max(i_list)+1 )
#ax_dr.set_ylim( ylims_dr[0], ylims_dr[1] )

fig_dr.tight_layout()
plt.savefig( fig_filename_dr )
print( f"Figure saved:\t{ fig_filename_dr }" )




#----- Plot abs rel err of second derivative as a function of r -----
fig_dr2 = plt.figure( figsize=(8,5) )
ax_dr2  = fig_dr2.add_subplot( 111 )

ax_dr2.semilogy( i_list, f_dr2_abs_rel_err_list, color="k", ls="-" )

ax_dr2.set_ylabel( "Absolute relative error" )
ax_dr2.set_xlabel( r"$r$-index $i$" )

ax_dr2.set_xlim( min(i_list)-1, max(i_list)+1 )
#ax_dr2.set_ylim( ylims_dr2[0], ylims_dr[1] )

fig_dr2.tight_layout()
plt.savefig( fig_filename_dr2 )
print( f"Figure saved:\t{ fig_filename_dr2 }" )
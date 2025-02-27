'''
Test_Chebyshev_series_arbitrary_interval.py

Plot the calculated and exact first and second radial derivatives of a function,
where the calculation method is by finite differencing.
'''

import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename_base = "../CSV/20240405_Test_Chebyshev_series_arbitrary_interval_nmax_NNN_one_over_r3_function_values.csv"

n_list = [ 8, 10, 12 ]

fig_filename_base = "../Figures/20240405/20240405_Test_Chebyshev_series_arbitrary_interval"
fig_filename_f   = fig_filename_base + "_f"
fig_filename_dr  = fig_filename_base + "_dr"
fig_filename_dr2 = fig_filename_base + "_dr2"

i_list                 = [ [] for n in n_list ]
f_exact_list           = [ [] for n in n_list ]
f_calc_list            = [ [] for n in n_list ]
f_dr_exact_list        = [ [] for n in n_list ]
f_dr_calc_list         = [ [] for n in n_list ]
f_dr_abs_rel_err_list  = [ [] for n in n_list ]	# Absolute relative error
f_dr2_exact_list       = [ [] for n in n_list ]
f_dr2_calc_list        = [ [] for n in n_list ]
f_dr2_abs_rel_err_list = [ [] for n in n_list ]	# Absolute relative error

ylims_dr = [ 0, 10 ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Read CSV files for stdev as a function of r -----
for n_i in range( len( n_list ) ):
	
	csv_filename = csv_filename_base.replace( "NNN", str( n_list[n_i] ) )
	
	with open( csv_filename, "r" ) as input_file:
		
		data = csv.reader( input_file )
		
		headers = next( data )
		
		# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
		# but it's better than obtaining the wrong column without checking.
		i_i           = headers.index( "i"           )
		i_f_exact     = headers.index( "f_dr_exact"  )
		i_f_calc      = headers.index( "f_dr_calc"   )
		i_f_dr_exact  = headers.index( "f_dr_exact"  )
		i_f_dr_calc   = headers.index( "f_dr_calc"   )
		#i_f_dr2_exact = headers.index( "f_dr2_exact" )
		#i_f_dr2_calc  = headers.index( "f_dr2_FD"    )
		
		for row in data:
			i_list          [n_i].append( int  ( row[ i_i           ] ) )
			f_exact_list    [n_i].append( float( row[ i_f_exact     ] ) )
			f_calc_list     [n_i].append( float( row[ i_f_calc      ] ) )
			f_dr_exact_list [n_i].append( float( row[ i_f_dr_exact  ] ) )
			f_dr_calc_list  [n_i].append( float( row[ i_f_dr_calc   ] ) )
			#f_dr2_exact_list[n_i].append( float( row[ i_f_dr2_exact ] ) )
			#f_dr2_calc_list [n_i].append( float( row[ i_f_dr2_calc  ] ) )
		



#----- Calculate absolute relative error -----
f_abs_rel_err_list     = np.zeros( shape = ( len( n_list ), len( i_list[0] ) ) )
f_dr_abs_rel_err_list  = np.zeros( shape = ( len( n_list ), len( i_list[0] ) ) )
#f_dr2_abs_rel_err_list = np.zeros( shape = ( len( n_list ), len( i_list[0] ) ) )

for n_i in range( len( n_list ) ):
	for i in range( len( i_list[n_i] ) ):
		
		if f_exact_list[n_i][i] != 0:
			f_abs_rel_err_list[n_i][i]     = abs( 1.0 - f_calc_list    [n_i][i] / f_exact_list    [n_i][i] )
			
		if f_dr_exact_list[n_i][i] != 0:
			f_dr_abs_rel_err_list[n_i][i]  = abs( 1.0 - f_dr_calc_list [n_i][i] / f_dr_exact_list [n_i][i] )
		
		#if f_dr2_exact_list[n_i][i] != 0 :
		#	f_dr2_abs_rel_err_list[n_i][i] = abs( 1.0 - f_dr2_calc_list[n_i][i] / f_dr2_exact_list[n_i][i] )




#----- Plot abs rel err of function as a function of r -----
fig_f = plt.figure( figsize=(8,5) )
ax_f  = fig_f.add_subplot( 111 )

for n_i in range( len( n_list ) ):
	ax_f.semilogy( i_list[n_i], f_abs_rel_err_list[n_i], label=r"$n="+str(n_list[n_i])+"$" )

ax_f.legend( loc="best" )
ax_f.set_ylabel( "Relative error" )
ax_f.set_xlabel( r"$r$-index $i$" )

ax_f.set_xlim( min(i_list[0])-1, max(i_list[0])+1 )
#ax_f.set_ylim( ylims_f[0], ylims_f[1] )

fig_f.tight_layout()
plt.savefig( fig_filename_f )
print( f"Figure saved:\t{ fig_filename_f }" )




#----- Plot abs rel err of first derivative as a function of r -----
fig_dr = plt.figure( figsize=(8,5) )
ax_dr  = fig_dr.add_subplot( 111 )

for n_i in range( len( n_list ) ):
	ax_dr.semilogy( i_list[n_i], f_dr_abs_rel_err_list[n_i], label=r"$n="+str(n_list[n_i])+"$" )

ax_dr.legend( loc="best" )
ax_dr.set_ylabel( "Relative error" )
ax_dr.set_xlabel( r"$r$-index $i$" )

ax_dr.set_xlim( min(i_list[0])-1, max(i_list[0])+1 )
#ax_dr.set_ylim( ylims_dr[0], ylims_dr[1] )

fig_dr.tight_layout()
plt.savefig( fig_filename_dr )
print( f"Figure saved:\t{ fig_filename_dr }" )



'''
#----- Plot abs rel err of second derivative as a function of r -----
fig_dr2 = plt.figure( figsize=(8,5) )
ax_dr2  = fig_dr2.add_subplot( 111 )

for n_i in range( len( n_list ) ):
	ax_dr2.semilogy( i_list[n_i], f_dr2_abs_rel_err_list[n_i], label=r"$n="+str(n_list[n_i])+"$" )

ax_dr2.legend( loc="best" )
ax_dr2.set_ylabel( "Relative error" )
ax_dr2.set_xlabel( r"$r$-index $i$" )

ax_dr2.set_xlim( min(i_list[0])-1, max(i_list[0])+1 )
#ax_dr.set_ylim( ylims_dr[0], ylims_dr[1] )

fig_dr2.tight_layout()
plt.savefig( fig_filename_dr2 )
print( f"Figure saved:\t{ fig_filename_dr2 }" )
'''
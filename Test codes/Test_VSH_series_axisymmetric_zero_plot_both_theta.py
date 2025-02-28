'''
Test_VSH_series_axisymmetric_zero_plot_both_theta.py

Plot stdev of recalculated vector as a function of r, and then as a function of t.

In this code, we compare data generated from two different methods of producing
VSH series of divergenceless functions.
Chebyshev method removed - see archive for previous version which included Chebsyhev method.

Automatic ylim calculation removed and replaced by manual setting. See original code in archive to recover it.
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename_base_original = "../CSV/20240502_Test_VSH_series_axisymmetric_General"
csv_filename_base_FD       = "../CSV/20240502_Test_VSH_series_axisymmetric_Divergenceless_FD"

fig_filename_base = "../Figures/20240502/20240502_Test_VSH_series_axisymmetric_theta_component"
fig_filename_r = fig_filename_base + "_stdev_r"
fig_filename_t = fig_filename_base + "_stdev_t"

csv_filenames_r = []
csv_filenames_t = []

i_list    = [ [], [] ]	# Number of values not known a priori. Not strictly necessary as can use range(len(stdev_r)) but an extra safety net in case i not ordered.
stdev_r_t = [ [], [] ]
j_list    = [ [], [] ]
stdev_t_t = [ [], [] ]

ylims_r = [ 1e-9 , 1e-5 ]
ylims_t = [ 1e-10, 1e-5 ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })


#----- Build lists of CSV filenames -----
for filename in [ csv_filename_base_original, csv_filename_base_FD ]:
	csv_filenames_r.append( filename + "_values_stdev_r.csv" )
	csv_filenames_t.append( filename + "_values_stdev_t.csv" )




#----- Read CSV files for stdev as a function of r -----
for n, filename in enumerate( csv_filenames_r ):
	with open( filename, "r" ) as input_file_r:
		
		data_r = csv.reader( input_file_r )
		
		headers_r = next( data_r )
		
		# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
		# but it's better than obtaining the wrong column without checking.
		i_i         = headers_r.index( "i"       )
		i_stdev_r_t = headers_r.index( "stdev_t" )
		
		for row in data_r:
			i_list   [n].append( int  ( row[ i_i       ] ) )
			stdev_r_t[n].append( float( row[ i_stdev_r_t ] ) )




#----- Plot stdev as a function of r -----
fig_r = plt.figure( figsize=(8,5) )
ax_r  = fig_r.add_subplot( 111 )

ax_r.semilogy( i_list[0], stdev_r_t[0], color="k", ls="-" , label=r"General"           )
ax_r.semilogy( i_list[1], stdev_r_t[1], color="b", ls="--", label=r"Divergenceless FD" )

ax_r.legend( loc="best" )
ax_r.set_ylabel( "Standard deviation" )
ax_r.set_xlabel( r"$r$-index $i$" )

ax_r.set_xlim( min(i_list[0]), max(i_list[0]) )
ax_r.set_ylim( ylims_r[0], ylims_r[1] )

fig_r.tight_layout()
plt.savefig( fig_filename_r )
print( f"Figure saved:\t{ fig_filename_r }" )




#----- Read CSV files for stdev as a function of theta -----
for n, filename in enumerate( csv_filenames_t ):
	with open( filename, "r" ) as input_file_t:
		
		data_t = csv.reader( input_file_t )
		
		headers_t = next( data_t )
		
		# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
		# but it's better than obtaining the wrong column without checking.
		i_j         = headers_t.index( "j"       )
		i_stdev_t_t = headers_t.index( "stdev_t" )
		
		for row in data_t:
			j_list   [n].append( int  ( row[ i_j       ] ) )
			stdev_t_t[n].append( float( row[ i_stdev_t_t ] ) )




#----- Plot stdev as a function of theta -----
fig_t = plt.figure( figsize=(8,5) )
ax_t  = fig_t.add_subplot( 111 )

ax_t.semilogy( j_list[0], stdev_t_t[0], color="k", ls="-" , label=r"General"           )
ax_t.semilogy( j_list[1], stdev_t_t[1], color="b", ls="--", label=r"Divergenceless FD" )

ax_t.legend( loc="best" )
ax_t.set_ylabel( "Standard deviation" )
ax_t.set_xlabel( r"$\theta$-index $j$" )

ax_t.set_xlim( min(i_list[0]), max(i_list[0]) )
ax_t.set_ylim( ylims_t[0], ylims_t[1] )

fig_t.tight_layout()
plt.savefig( fig_filename_t )

print( f"Figure saved:\t{ fig_filename_t }" )
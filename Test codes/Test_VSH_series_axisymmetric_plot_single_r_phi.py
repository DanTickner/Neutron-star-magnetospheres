'''
Test_VSH_series_axisymmetric_plot_single_r_phi.py

Plot stdev of recalculated vector as a function of r, and then as a function of t.
Slightly modified 20240309 to skip the theta terms, so that those can all be plotted together instead.
To include the theta terms, un-comment the respective code lines.

Automatic ylim calculation removed and replaced by manual setting. See original code in archive to recover it.
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename_base = "../CSV/20240404_Test_VSH_series_axisymmetric_Original"
fig_filename_base = "../Figures/20240404/20240404_Test_VSH_series_axisymmetric_r_and_phi_components"

csv_filename_r = csv_filename_base + "_values_stdev_r.csv"
csv_filename_t = csv_filename_base + "_values_stdev_t.csv"

i_list    = []	# Number of values not known a priori. Not strictly necessary as can use range(len(stdev_r)) but an extra safety net in case i not ordered.
stdev_r_r = []	
stdev_r_t = []
stdev_r_p = []
j_list    = []
stdev_t_r = []	
stdev_t_t = []
stdev_t_p = []

ylims_r = [ 1e-9, 1e-4 ]
ylims_t = [ 1e-9, 1e-4 ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })



#----- Read CSV for stdev as a function of r -----
with open( csv_filename_r, "r" ) as input_file_r:
	
	data_r = csv.reader( input_file_r )
	
	headers_r = next( data_r )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
	# but it's better than obtaining the wrong column without checking.
	i_i         = headers_r.index( "i"       )
	i_stdev_r_r = headers_r.index( "stdev_r" )
	i_stdev_r_t = headers_r.index( "stdev_t" )
	i_stdev_r_p = headers_r.index( "stdev_p" )
	
	for row in data_r:
		i_list   .append( int  ( row[ i_i       ] ) )
		stdev_r_r.append( float( row[ i_stdev_r_r ] ) )
		stdev_r_t.append( float( row[ i_stdev_r_t ] ) )
		stdev_r_p.append( float( row[ i_stdev_r_p ] ) )




#----- Plot stdev as a function of r -----
fig_r = plt.figure( figsize=(8,5) )
ax_r  = fig_r.add_subplot( 111 )

ax_r.semilogy( i_list, stdev_r_r, color="k", ls="-" , label=r"$r$"      )
#ax_r.semilogy( i_list, stdev_r_t, color="b", ls="--", label=r"$\theta$" )
ax_r.semilogy( i_list, stdev_r_p, color="r", ls=":" , label=r"$\phi$"   )

ax_r.legend( loc="best" )
ax_r.set_ylabel( "Standard deviation" )
ax_r.set_xlabel( r"$r$-index $i$" )

ax_r.set_xlim( min(i_list), max(i_list) )
ax_r.set_ylim( ylims_r[0], ylims_r[1] )

fig_r.tight_layout()
plt.savefig( fig_filename_base + "_stdev_r" )




#----- Read CSV for stdev as a function of theta -----
with open( csv_filename_t, "r" ) as input_file_t:
	
	data_t = csv.reader( input_file_t )
	
	headers_t = next( data_t )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and j here,
	# but it's better than obtaining the wrong column without checking.
	i_j         = headers_t.index( "j"       )
	i_stdev_t_r = headers_t.index( "stdev_r" )
	i_stdev_t_t = headers_t.index( "stdev_t" )
	i_stdev_t_p = headers_t.index( "stdev_p" )
	
	for row in data_t:
		j_list   .append( int  ( row[ i_j       ] ) )
		stdev_t_r.append( float( row[ i_stdev_t_r ] ) )
		stdev_t_t.append( float( row[ i_stdev_t_t ] ) )
		stdev_t_p.append( float( row[ i_stdev_t_p ] ) )




#----- Plot stdev as a function of theta -----
fig_t = plt.figure( figsize=(8,5) )
ax_t  = fig_t.add_subplot( 111 )

ax_t.semilogy( j_list, stdev_t_r, color="k", ls="-" , label=r"$r$"      )
#ax_t.semilogy( j_list, stdev_t_t, color="b", ls="--", label=r"$\theta$" )
ax_t.semilogy( j_list, stdev_t_p, color="r", ls=":" , label=r"$\phi$"   )

ax_t.legend( loc="best" )
ax_t.set_ylabel( "Standard deviation" )
ax_t.set_xlabel( r"$\theta$-index $j$" )

ax_t.set_xlim( min(j_list), max(j_list) )
ax_t.set_ylim( ylims_t[0], ylims_t[1] )

fig_t.tight_layout()
plt.savefig( fig_filename_base + "_stdev_t" )
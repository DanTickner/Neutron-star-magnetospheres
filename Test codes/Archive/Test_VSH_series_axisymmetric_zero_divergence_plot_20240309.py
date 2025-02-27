'''
Test_VSH_series_axisymmetric_zero_divergence_plot.py

Plot stdev of recalculated vector as a function of r, and then as a function of t.

In this code, we compare data generated from three different methods of producing
VSH series of divergenceless functions.
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename_base_original  = "../CSV/20240309_Test_VSH_series_axisymmetric_Original"
csv_filename_base_FD        = "../CSV/20240309_Test_VSH_series_axisymmetric_Divergenceless_FD"
csv_filename_base_Chebyshev = "../CSV/20240309_Test_VSH_series_axisymmetric_Divergenceless_Chebyshev"

fig_filename_base = "../Figures/20240309/20240309_Test_VSH_series_axisymmetric_zero_divergence"

csv_filenames_r = []
csv_filenames_t = []

i_list    = [ [], [], [] ]	# Number of values not known a priori. Not strictly necessary as can use range(len(stdev_r)) but an extra safety net in case i not ordered.
stdev_r_t = [ [], [], [] ]
j_list    = [ [], [], [] ]
stdev_t_t = [ [], [], [] ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })


#----- Build lists of CSV filenames -----
for filename in [ csv_filename_base_original, csv_filename_base_FD, csv_filename_base_Chebyshev ]:
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




#----- Determine limits for graph for stdev as a function of r -----
xmax_r = max( [ max( lst ) for lst in i_list ] )
xmin_r = min( [ min( lst ) for lst in i_list ] )

stdev_r_without_zeros = []
for lst in stdev_r_t:
	for stdev in lst:
		if stdev != 0:
			stdev_r_without_zeros.append( stdev )

if stdev_r_without_zeros == []:
	ymax_r = 0.1
	ymin_r = 0.01
else:
	ymax_r = 10** int( np.log10( max( stdev_r_without_zeros ) ) )
	ymin_r = 10** int( np.log10( min( stdev_r_without_zeros ) ) - 1 )





#----- Plot stdev as a function of r -----
fig_r = plt.figure( figsize=(8,5) )
ax_r  = fig_r.add_subplot( 111 )

ax_r.semilogy( i_list[0], stdev_r_t[0], color="k", ls="-" , label=r"Original"  )
ax_r.semilogy( i_list[1], stdev_r_t[1], color="b", ls="--", label=r"FD"        )
ax_r.semilogy( i_list[2], stdev_r_t[2], color="r", ls=":" , label=r"Chebyshev" )

ax_r.legend( loc="best" )
ax_r.set_ylabel( "Standard deviation" )
ax_r.set_xlabel( r"$r$-index $i$" )

ax_r.set_xlim( xmin_r, xmax_r )
ax_r.set_ylim( ymin_r, ymax_r )

fig_r.tight_layout()
plt.savefig( fig_filename_base + "_r" )
print( f"Figure saved:\t{ fig_filename_base + '_r' }" )




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




#----- Determine limits for graph for stdev as a function of theta -----
xmax_t = max( [ max( lst ) for lst in j_list ] )
xmin_t = min( [ min( lst ) for lst in j_list ] )

stdev_t_without_zeros = []
for lst in stdev_t_t:
	for stdev in lst:
		if stdev != 0:
			stdev_t_without_zeros.append( stdev )

if stdev_t_without_zeros == []:
	ymax_t = 0.1
	ymin_t = 0.01
else:
	ymax_t = 10** int( np.log10( max( stdev_t_without_zeros ) ) )
	ymin_t = 10** int( np.log10( min( stdev_t_without_zeros ) ) - 1 )




#----- Plot stdev as a function of theta -----
fig_t = plt.figure( figsize=(8,5) )
ax_t  = fig_t.add_subplot( 111 )

ax_t.semilogy( j_list[0], stdev_t_t[0], color="k", ls="-" , label=r"Original"  )
ax_t.semilogy( j_list[1], stdev_t_t[1], color="b", ls="--", label=r"FD"        )
ax_t.semilogy( j_list[2], stdev_t_t[2], color="r", ls=":" , label=r"Chebyshev" )

ax_t.legend( loc="best" )
ax_t.set_ylabel( "Standard deviation" )
ax_t.set_xlabel( r"$\theta$-index $j$" )

ax_t.set_xlim( xmin_t, xmax_t )
ax_t.set_ylim( ymin_t, ymax_t )

fig_r.tight_layout()
plt.savefig( fig_filename_base + "_t" )

print( f"Figure saved:\t{ fig_filename_base + '_t' }" )
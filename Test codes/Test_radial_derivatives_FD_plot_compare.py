'''
Test_radial_derivatives_FD_plot_compare.py

Plot the calculated and exact first and second radial derivatives of a function,
where the calculation method is by finite differencing.

Slightly modified code to compare multiple CSVs and save having two separate figures in the report.
'''

import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filenames = [
	#"../CSV/20240502_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3_order_2_function_values.csv",
	#"../CSV/20240502_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3_order_4_function_values.csv",
	"../CSV/20240502_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3_order_6_function_values.csv",
	"../CSV/20240502_Test_radial_derivatives_FD_enhanced_endpoints_1_over_r3_order_6_with_8_near_endpoints_function_values.csv"
]

labels = [
	#r"$\mathcal{O}(h^2)$",
	#r"$\mathcal{O}(h^4)$",
	r"$\mathcal{O}(h^6)$",
	r"$\mathcal{O}(h^6-h^8)$ hybrid"
]

colors = [
	"k",
	"b",
	"r"
]

lss = [
	"-",
	"--",
	":"
]

fig_filename_base = "../Figures/20240502/20240502_Test_radial_derivatives_FD_hybrid_outer"
fig_filename_dr  = fig_filename_base + "_dr"
fig_filename_dr2 = fig_filename_base + "_dr2"

r_list         = [ [] for filename in csv_filenames ]	# Radial positions may be different for each CSV.
f_dr_ARE_list  = [ [] for filename in csv_filenames ]	# Absolute relative error
f_dr2_ARE_list = [ [] for filename in csv_filenames ]	# Absolute relative error

use_custom_xlim = True
custom_xlim_dr  = [ 9.44, 9.56 ]
custom_xlim_dr2 = [ 9.44, 9.56 ]

use_custom_ylim = True
custom_ylim_dr  = [ 1e-14, 1e-11 ]
custom_ylim_dr2 = [ 1e-13, 1e-10 ]




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Read CSV files for ARE as a function of r -----
for n, filename in enumerate( csv_filenames ):
	with open( filename, "r" ) as input_file:
		
		data = csv.reader( input_file )
		
		headers = next( data )
		
		# Get correct columns.
		i_r         = headers.index( "r"                 )
		i_f_dr_ARE  = headers.index( "abs_rel_err_f_dr"  )
		i_f_dr2_ARE = headers.index( "abs_rel_err_f_dr2" )
		
		for row in data:
			r_list        [n].append( float( row[ i_r         ] ) )
			f_dr_ARE_list [n].append( float( row[ i_f_dr_ARE  ] ) )
			f_dr2_ARE_list[n].append( float( row[ i_f_dr2_ARE ] ) )


# Extend r either side so the endpoints don't get obscured by the boundaries of the figure.
r_min         = min( [ min( lst ) for lst in r_list         ] ) * 0.98
r_max         = max( [ max( lst ) for lst in r_list         ] ) * 1.02
f_dr_ARE_min  = min( [ min( lst ) for lst in f_dr_ARE_list  ] )
f_dr_ARE_max  = max( [ max( lst ) for lst in f_dr_ARE_list  ] )
f_dr2_ARE_min = min( [ min( lst ) for lst in f_dr2_ARE_list ] )
f_dr2_ARE_max = max( [ max( lst ) for lst in f_dr2_ARE_list ] )

ylims_dr = [
	10**( np.floor( np.log10( f_dr_ARE_min ) ) ),
	10**( np.ceil ( np.log10( f_dr_ARE_max ) ) )
]

ylims_dr2 = [
	10**( np.floor( np.log10( f_dr2_ARE_min ) ) ),
	10**( np.ceil ( np.log10( f_dr2_ARE_max ) ) )
]

# Unfortunarely the default xticks won't include r=1 and r=10 and don't look right.
xticks_values = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
xticks_labels = [ rf"${ x }$" for x in xticks_values ]


#----- Plot abs rel err of first derivative as a function of r -----
fig_dr = plt.figure( figsize=(8,5) )
ax_dr  = fig_dr.add_subplot( 111 )

for n in range(len( csv_filenames )):
	ax_dr.semilogy( r_list[n], f_dr_ARE_list[n], label=labels[n], color=colors[n], ls=lss[n] )

ax_dr.legend( loc="best" )

ax_dr.set_ylabel( "Absolute relative error" )
ax_dr.set_xlabel( r"$r$" )

if use_custom_xlim:
	ax_dr.set_xlim( custom_xlim_dr[0], custom_xlim_dr[1] )
else:
	ax_dr.set_xlim( r_min, r_max )

if use_custom_ylim:
	ax_dr.set_ylim( custom_ylim_dr[0], custom_ylim_dr[1] )
else:
	ax_dr.set_ylim( ylims_dr[0], ylims_dr[1] )

#plt.xticks( xticks_values, xticks_labels )

fig_dr.tight_layout()
plt.savefig( fig_filename_dr )
print( f"Figure saved:\t{ fig_filename_dr }" )




#----- Plot abs rel err of second derivative as a function of r -----
fig_dr2 = plt.figure( figsize=(8,5) )
ax_dr2  = fig_dr2.add_subplot( 111 )

for n in range(len( csv_filenames )):
	ax_dr2.semilogy( r_list[n], f_dr2_ARE_list[n], label=labels[n], color=colors[n], ls=lss[n] )

ax_dr2.legend( loc="best" )

ax_dr2.set_ylabel( "Absolute relative error" )
ax_dr2.set_xlabel( r"$r$" )

if use_custom_xlim:
	ax_dr2.set_xlim( custom_xlim_dr2[0], custom_xlim_dr2[1] )
else:
	ax_dr2.set_xlim( r_min, r_max )

if use_custom_ylim:
	ax_dr2.set_ylim( custom_ylim_dr2[0], custom_ylim_dr2[1] )
else:
	ax_dr2.set_ylim( ylims_dr[0], ylims_dr[1] )

#plt.xticks( xticks_values, xticks_labels )

fig_dr2.tight_layout()
plt.savefig( fig_filename_dr2 )
print( f"Figure saved:\t{ fig_filename_dr2 }" )
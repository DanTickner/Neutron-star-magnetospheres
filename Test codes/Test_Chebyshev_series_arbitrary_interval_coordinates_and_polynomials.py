'''
Test_Chebyshev_series_arbitrary_interval_coordinates_and_polynomials.py

Read the output CSV containing values of transformed coordinates and Chebyshev polynomials.
Plot the polynomials as a function of r, and of R, to verify.
'''


import csv
import numpy as np
import matplotlib.pyplot as plt


#----- Define variables -----
csv_filename  = "../CSV/20240228_Test_Chebyshev_series_arbitrary_interval_coordinates.csv"

fig_filename = "../Figures/20240228/20240228_Test_Chebyshev_series_arbitrary_interval"

only_plot_chosen_n = True # Graph becomes quickly crowded. Only select all
n_to_plot = [ 0, 1, 2, 3, 4 ]

i_list              = []
r_list              = []
lambda_inverse_list = []
R_list              = []
T_n_list            = []




#----- Global matplotlib parameters -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size": 18 })




#----- Read CSV file -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	for header in headers:
		if header[:2] == "T_":
			T_n_list.append( [] )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and i here,
	# but it's better than obtaining the wrong column without checking.
	i_i              = headers.index( "i"              )
	i_r              = headers.index( "r"              )
	i_lambda_inverse = headers.index( "Lambda_inverse" )
	i_R              = headers.index( "R"              )
	i_T_0            = headers.index( "T_0"            )
	
	for row in data:
		i_list             .append( int  ( row[ i_i              ] ) )
		r_list             .append( float( row[ i_r              ] ) )
		lambda_inverse_list.append( float( row[ i_lambda_inverse ] ) )
		R_list             .append( float( row[ i_R              ] ) )
		T_n_list[0]        .append( float( row[ i_T_0            ] ) )
		
		for n in range( 1, len( T_n_list ) ):
			T_n_list[n].append( float( row[ i_T_0 + n ] ) )




'''
#----- Plot r, lambda_inverse and R as a function of i -----
fig_i = plt.figure( figsize=(8,5) )
ax_i  = fig_i.add_subplot( 111 )

ax_i.plot( i_list, r_list             , color="k", ls="-" , label=r"$r$"                                  )
ax_i.plot( i_list, lambda_inverse_list, color="b", ls="--", label=r"$\Lambda^{-1}(r)$"                    )
ax_i.plot( i_list, R_list             , color="r", ls=":" , label=r"$\cos^{-1}\big[\Lambda^{-1}(r)\big]$" )

ax_i.legend( loc="best" )

ax_i.set_xlabel( r"Radial index $i$" )
ax_i.set_ylabel( r"Value" )
'''




#----- Plot T_n as a function of r -----
fig_T_r = plt.figure( figsize=(8,5) )
ax_T_r  = fig_T_r.add_subplot( 111 )

for n in range( len( T_n_list ) ):
	if ( not(only_plot_chosen_n) or ( only_plot_chosen_n and n in n_to_plot ) ):
		ax_T_r.plot( r_list, T_n_list[n], label=f"$n={n}$" )

ax_T_r.legend( loc="best" )

ax_T_r.set_xlabel( r"Radial coordinate $r$" )
ax_T_r.set_ylabel( r"Chebyshev polynomial $T_n(r)$" )




#----- Plot T_n as a function of L=Lambda_inverse -----
fig_T_L = plt.figure( figsize=(8,5) )
ax_T_L  = fig_T_L.add_subplot( 111 )

for n in range( len( T_n_list ) ):
	if ( not(only_plot_chosen_n) or ( only_plot_chosen_n and n in n_to_plot ) ):
		ax_T_L.plot( lambda_inverse_list, T_n_list[n], label=f"$n={n}$" )

ax_T_L.legend( loc="best" )

ax_T_L.set_xlabel( r"Mapped radial coordinate $\Lambda^{-1}(r)$" )
ax_T_L.set_ylabel( r"Chebyshev polynomial $T_n\big[\Lambda^{-1}(r)\big]$" )




'''
#----- Plot T_n as a function of R -----
fig_T_R = plt.figure( figsize=(8,5) )
ax_T_R  = fig_T_R.add_subplot( 111 )

for n in range( len( T_n_list ) ):
	if ( not(only_plot_chosen_n) or ( only_plot_chosen_n and n in n_to_plot ) ):
		ax_T_R.plot( R_list, T_n_list[n], label=f"$n={n}$" )

ax_T_R.legend( loc="best" )

ax_T_R.set_xlabel( r"Chebyshev series integration variable $R=\cos^{-1}\big[\Lambda^{-1}(r)\big]$" )
ax_T_R.set_ylabel( r"Chebyshev polynomial $T_n(R)$" )
'''
'''
20241015_Test_twist_smoothing_function_plot.py

Read the smoothed twist values calculated by
20241015_Test_twist_smoothing_function.cpp
from a saved CSV file, and plot them in a heatmap.
'''

import numpy as np
import csv
import matplotlib.pyplot as plt


#----- Define variables -----

csv_filepath = "../CSV/20250217_Test_twist_smoothing_function_values.csv"
fig_filepath = "../Figures/20250217/20250217_Test_twist_smoothing_function_values_exponential"

r     = []
t     = []
twist = []

n_points_r = 0
n_points_t = 0

twist_max_value = 0

plot_lines_defining_region = True




#----- Read CSV first time: get number of datapoints -----

with open( csv_filepath, "r" ) as csv_file:
	
	data = csv.reader( csv_file )
	
	region_points = next( data )
	r_1 = float( region_points[ region_points.index( "r_1" ) + 1 ] )
	r_2 = float( region_points[ region_points.index( "r_2" ) + 1 ] )
	t_1 = float( region_points[ region_points.index( "t_1" ) + 1 ] )
	t_2 = float( region_points[ region_points.index( "t_2" ) + 1 ] )
	print( f"r_1:\t{ r_1 }" )
	print( f"r_2:\t{ r_2 }" )
	print( f"t_1:\t{ t_1 }" )
	print( f"t_2:\t{ t_2 }" )
	
	headers = next( data )
	
	i_i           = headers.index( "i" )
	i_j           = headers.index( "j" )
	i_r           = headers.index( "r" )
	i_t           = headers.index( "t" )
	i_twist       = headers.index( "twist" )
	i_twist_ratio = headers.index( "twist_over_unsmoothed" );
	
	for row in data:
		
		i = int( row[ i_i ] )
		j = int( row[ i_j ] )
		
		if i+1 > len( r ):
			r.append( float( row[ i_r ] ) )
		
		if j+1 > len( t ):
			t.append( float( row[ i_t ] ) )



#----- Make arrays of the correct size -----
n_points_r = len( r )
n_points_t = len( t )

r_grid           = np.zeros( shape=( n_points_r, n_points_t ) )
t_grid           = np.zeros( shape=( n_points_r, n_points_t ) )
twist_grid       = np.zeros( shape=( n_points_r, n_points_t ) )
twist_ratio_grid = np.zeros( shape=( n_points_r, n_points_t ) )



#----- Read CSV second time: Populate twist grid -----

with open( csv_filepath, "r" ) as csv_file:
	
	data = csv.reader( csv_file )
	
	next( data ) # Skip row containing twisted region boundary coordinates.
	next( data ) # SKip headers
	
	for row in data:
		
		i = int( row[ i_i ] )
		j = int( row[ i_j ] )
		
		r_grid          [i][j] = float( row[i_r]           )
		t_grid          [i][j] = float( row[i_t]           )
		twist_grid      [i][j] = float( row[i_twist]       )
		twist_ratio_grid[i][j] = float( row[i_twist_ratio] )
		
		if twist_grid[i][j] > twist_max_value:
			twist_max_value = twist_grid[i][j]

print( f"TWist max:\t{ twist_max_value }" )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure( figsize=(7,5) )
ax  = fig.add_subplot( 111 )

cp = ax.contourf( r_grid, t_grid, twist_ratio_grid )

ax.set_xlabel( r"$r$ (NS radii)" )
ax.set_ylabel( r"$\theta$ (rad)" )

fig.colorbar( cp )

if plot_lines_defining_region:
	ax.plot( [ r_1, r_1 ] , [ t_1, t_2 ], color="grey", ls="--", lw=1 )
	ax.plot( [ r_2, r_2 ] , [ t_1, t_2 ], color="grey", ls="--", lw=1 )
	ax.plot( [ r_1, r_2 ] , [ t_1, t_1 ], color="grey", ls="--", lw=1 )
	ax.plot( [ r_1, r_2 ] , [ t_2, t_2 ], color="grey", ls="--", lw=1 )
	#ax.plot( [ x_p_minus, x_p_minus ] , [ -1     , f(x_p_minus) ], color="grey", ls="--", lw=0.5 )

plt.savefig( fig_filepath, bbox_inches="tight" )
print( f"\nSaved figure:\t{ fig_filepath }" )
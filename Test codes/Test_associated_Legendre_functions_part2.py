'''
Test_Associated_Legendre_Functions_Part2.py

Read a CSV file containing calculated values of associated Legendre functions for m=0 and m=1,
with ell=0 up to some ell_max which this code determines simply by the highest ell found,
and for discretely spaced coordinates theta_j.
Calculate the functions separately using scipy modules.
For each ell, create two plots for m=0 and two plots for m=1.
1) C++ value of P_ell^m as a function of theta, and Python value of P_ell^m as a function of theta, on the same axis.
2) Residual calculated by subtracting the C++ value from the Python value, as a function of theta. This should be as close to zero as possible (assuming that the Python value is 100% accurate, which may not be the case).
Also calculate standard deviations for m=0 and m=1, for each ell.
This code features a lot of overwriting values of theta, but it's safe because we can assume that we have already ran
Test_coordinates.cpp to check that the coordinates are correct.

https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_legendre.html
'''


import csv
import numpy as np
from scipy.special import eval_legendre, lpmv


#----- Define variables -----
csv_filename = "../CSV/20240223_Test_Associated_Legendre_Functions.csv"

j_max   = 0
ell_max = 0




#----- Read CSV first time to get n_points_t and ell_max -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	# Get correct columns. Unfortunate that we have a potentially confusing conflict between i and j here,
	# but it's better than obtaining the wrong column without checking.
	i_j         = headers.index( "j"       )
	i_t_j       = headers.index( "t_j"     )
	i_ell       = headers.index( "ell"     )
	i_P_ell_0_j = headers.index( "P_ell_0" )
	i_P_ell_1_j = headers.index( "P_ell_1" )
	
	for row in data:
		j   = int( row[ i_j   ] )
		ell = int( row[ i_ell ] )
		
		if j > j_max:
			j_max = j
		
		if ell > ell_max:
			ell_max = ell



#----- Define arrays for theta and the functions, now that we know their sizes -----
n_points_t   = j_max   + 1
n_values_ell = ell_max + 1

t  = np.zeros( n_points_t ) # theta
ct = np.zeros( n_points_t ) # cos(theta)

P_ell_0_cpp = np.zeros( shape=( n_points_t, n_values_ell ) )
P_ell_1_cpp = np.zeros( shape=( n_points_t, n_values_ell ) )

P_ell_0_python = np.zeros( shape=( n_points_t, n_values_ell ) )
P_ell_1_python = np.zeros( shape=( n_points_t, n_values_ell ) )




#----- Read CSV second time to populate lists -----
with open( csv_filename, "r" ) as input_file:
	
	data = csv.reader( input_file )
	
	headers = next( data )
	
	for row in data:
		j         = int  ( row[ i_j         ] )
		ell       = int  ( row[ i_ell       ] )
		t_j       = float( row[ i_t_j       ] )
		P_ell_0_j = float( row[ i_P_ell_0_j ] )
		P_ell_1_j = float( row[ i_P_ell_1_j ] )
		
		P_ell_0_cpp[j][ell] = P_ell_0_j
		P_ell_1_cpp[j][ell] = P_ell_1_j
		
		if ell == 0:
			t [j] = t_j
			ct[j] = np.cos( t_j )





#----- Calculate the expected values of the associated Legendre functions -----
for j in range( n_points_t ):
	for ell in range( n_values_ell ):
		P_ell_0_python[j][ell] = eval_legendre( ell, ct[j] )
		P_ell_1_python[j][ell] = lpmv( 1, ell, ct[j] )




#----- Calculate stdev of functions for each ell -----
stdev_0 = np.zeros( n_values_ell )
stdev_1 = np.zeros( n_values_ell )

for j in range( n_points_t ):
	for ell in range( n_values_ell ):
		stdev_0[ell] += ( P_ell_0_cpp[j][ell] - P_ell_0_python[j][ell] )**2
		stdev_1[ell] += ( P_ell_1_cpp[j][ell] - P_ell_1_python[j][ell] )**2

print( "Standard deviation for associated Legendre functions" )
print( "ell\tm=0\tm=1" )

for ell in range( n_values_ell ):
	
	stdev_0[ell] = np.sqrt( stdev_0[ell] / n_points_t )
	stdev_1[ell] = np.sqrt( stdev_1[ell] / n_points_t )
	
	#print( f"{ ell }\t{ stdev_0[ell] }\t{ stdev_1[ell] }" )# For reading on the screen
	print( f"{ ell } & { stdev_0[ell] } & { stdev_1[ell] }\\\\" )# For copying into LaTeX table
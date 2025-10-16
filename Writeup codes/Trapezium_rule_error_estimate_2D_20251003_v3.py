'''
Trapezium_rule_error_estimate_2D_20251003_v3

Error estimate of my 2D trapezium rule.
Expect order( 1/N1^3 + 1/N2^3 )

Version 3: Read from CSV file instead of copying into lists.

Fit to f(N) = ( A1 * N1**k1 ) + ( A2 * N2**k2 )
with A1, A2, k1, k2 to be determined.
Expect k1, k2 around -3.
'''

import csv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

use_initial_guess = True
initial_guess = [ 0, -3, 0, -3 ]


#----- Input data from the C++ code executions in a CSV file -----
#csv_filename = "CSV/20251004_Test_trapezium_rule_2D_Cartesian.csv"
csv_filename = "CSV/20251004_Test_trapezium_rule_2D_Polar.csv"


N_points_1 = []
N_points_2 = []
abs_err    = []

with open( csv_filename ) as csv_file:
	
	data = csv.reader( csv_file )
	
	next( data ) # Skip first row with outputted parameters.
	next( data ) # Skip headers.
	
	for row in data:
		
		N_points_1.append( int  ( row[0] ) )
		N_points_2.append( int  ( row[1] ) )
		abs_err   .append( float( row[3] ) )


#----- Manually perform the power-law fit -----
def func( N_tuple, A1, k1, A2, k2 ):
	# The fitting function is probably the only thing to change between Cartesian and polar.
	N1, N2 = N_tuple
	return A1 * N1**k1 + A2 * N2**k2

if use_initial_guess:
	popt, pcov = curve_fit(
		func,
		( N_points_1, N_points_2 ),
		abs_err,
		p0=initial_guess
	)
else:
	popt, pcov = curve_fit(
		func,
		( N_points_1, N_points_2 ),
		abs_err
	)

perr = np.sqrt(np.diag( pcov )) # One-standard-deviation errors

abs_err_fit = []
abs_err_fit_rel_diff = [] # Relative difference compared to the CSV values of the absolute error.
for i in range(len( N_points_1 )):
	err  = popt[0] * N_points_1[i]**popt[1]
	err += popt[2] * N_points_2[i]**popt[3]
	abs_err_fit.append( err )
	abs_err_fit_rel_diff.append( abs( ( err - abs_err[i] ) / abs_err[i] ) )

#----- Output results -----
#for i in range( len( N_points_1 ) ):
#	print( f"{ N_points_1[i] }\t{ N_points_2[i] }\t{ abs_err[i] }\t{ abs_err_fit[i] }\t{ abs_err_fit_rel_diff }" )

print( "\nOptimised fitting parameters" )
print( popt )
print( "Associated uncertainties" )
print( perr )
'''
Trapezium_rule_error_estimate_1D_20251003_v2

Error estimate of my 1D trapezium rule.
Expect order( 1/N^3 )

Version 2: Rewrite from the 2D version, with a CSV reader instead of manually inputting data.

Fit to f(N) = ( A1 * N**k )
with A, k to be determined.
Expect k around -3.
'''

import csv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


#----- Input data from the C++ code executions in a CSV file -----
csv_filename = "CSV/20251004_Test_trapezium_rule_1D_Cartesian.csv"


N_points = []
abs_err  = []

with open( csv_filename ) as csv_file:
	
	data = csv.reader( csv_file )
	
	next( data ) # Skip first row with outputted parameters.
	next( data ) # Skip headers.
	
	for row in data:
		
		N_points.append( int  ( row[0] ) )
		abs_err .append( float( row[2] ) )


#----- Manually perform the power-law fit -----
def func( N, A, k ):
	# The fitting function is probably the only thing to change between Cartesian and polar.
	return A * N**k

popt, pcov = curve_fit(
	func,
	N_points,
	abs_err
)

perr = np.sqrt(np.diag( pcov )) # One-standard-deviation errors

abs_err_fit = []
abs_err_fit_rel_diff = [] # Relative difference compared to the CSV values of the absolute error.
for i in range(len( N_points )):
	err  = popt[0] * N_points[i]**popt[1]
	abs_err_fit.append( err )
	abs_err_fit_rel_diff.append( abs( ( err - abs_err[i] ) / abs_err[i] ) )

#----- Output results -----
for i in range( len( N_points ) ):
	print( f"{ N_points[i] }\t{ abs_err[i] }\t{ abs_err_fit[i] }\t{ abs_err_fit_rel_diff }" )

print( "\nOptimised fitting parameters" )
print( popt )
print( "Associated uncertainties" )
print( perr )
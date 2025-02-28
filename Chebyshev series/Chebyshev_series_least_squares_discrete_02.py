'''
Chebyshev_series_least_squares_discrete_02.py

Fit a set of discrete datapoints to a set of Chebyshev polynomials by least-squares fitting.
Useful for prototyping.

Modified from
20230613 VSH decomposition/Fit_function_to_Legendre_rational_functions.py
and
Chebyshev_series_least_squares_01.py

V02: Do everything explicitly as a way of prototyping for C++ code.

'''

import numpy as np


def T_n( x, n, a, b ):
	# Chebyshev polynompial of order n evalauted at x, scaled for the interval [a, b].
	Lambda_inverse = ( 2.0 - x - a - b ) / ( b - a )
	return np.cos( n * Lambda_inverse )


def gaussian_elimination( A, B ):
	
	#----- Elimination phase -----
	for j in range(len( B )):
		for i in range( j+1, len(B) ):
			
			if A[i][j] != 0:
				Lambda = A[i][j] / A[j][j]
				
				for k in range( j+1, len(B) ):
					A[i][k] -= Lambda * A[j][k]
					
				B[i] -= Lambda * B[j]
	
	#----- Back-substitution phase -----
	for i in range( len(B)-1, -1, -1 ):
		dotproduct = sum( [ A[i][j] * B[j] for j in range( i+1, len(B) ) ] )
		B[i] = ( B[i] - dotproduct ) / A[i][i]
	
	return B




#------ Data -----
n_max = 3		# Maximum order of Chebyshev polynomials to fit to.
x_min = 1		# Range over which the fit needs to be valid.
x_max = 2		# Range over which the fit needs to be valid.

#x_data = [ 1.234, 12.432, 45.32, 55.6 ]
x_data = [ 1.234, 1.432, 1.532, 1.866 ]
y_data = [ 3.0*x**3 - 1.5*x**2 for x in x_data ]



#----- Build matrix A -----
A = np.zeros( shape=(n_max+1,n_max+1) )

for n in range( n_max+1 ):
	for m in range( n_max+1 ):
		
		for i in range( len( x_data ) ):
			A[n][m] += T_n( x_data[i], n, x_min, x_max ) * T_n( x_data[i], m, x_min, x_max )


#----- Build vector B -----
B = np.zeros( n_max+1 )


for n in range( n_max+1 ):
	
	for i in range( len( x_data ) ):
		B[n] += y_data[i] * T_n( x_data[i], n, x_min, x_max )


#----- Perform Gaussian elimination -----
print( B )

gaussian_elimination( A, B )

print( B )


#----- Calculate fitted y-values -----
y_fit = np.zeros( len( x_data ) )

for i in range( len( x_data ) ):
	
	for n in range( n_max+1 ):
		y_fit[i] += B[n] * T_n( x_data[i], n, x_min, x_max )


#----- Output correct and fitted y-values -----

for i in range( len( x_data ) ):
	print( f"{ x_data[i] }\t{ y_data[i] }\t{ y_fit[i] }" )





'''
#----- Test Gaussian elimination function -----
// file:///C:/Users/Dan/OneDrive/PhD/TA%20work/CMP-4005Y%20Mathematics%20for%20Computing%20B/Exercises/06-matrix-exercises1-2022-09-12(2).pdf
// Q3c. Correct result is [ 1, 2, 1 ].

matrix_A = [ [ 1, 3, 1 ], [ 2, 1, 3 ], [ 1, 1, -1 ] ]
vector_B = [ 8, 7, 2 ]

gaussian_elimination( matrix_A, vector_B )
print( vector_B )

print( f"{ 1 * vector_B[0] + 3 * vector_B[1] + 1  * vector_B[2] }\t{ 8 }" )
print( f"{ 2 * vector_B[0] + 1 * vector_B[1] + 3  * vector_B[2] }\t{ 7 }" )
print( f"{ 1 * vector_B[0] + 1 * vector_B[1] + -1 * vector_B[2] }\t{ 2 }" )
'''



#----- Now try extrapolating outward a bit -----
x_extrapolate = 1.865

y_extrapolate_exact = 3.0*x_extrapolate**3 - 1.5*x_extrapolate**2

y_extrapolate_fit = 0
for n in range( n_max+1 ):
	y_extrapolate_fit += B[n] * T_n( x_extrapolate, n, x_min, x_max )

print( f"\nExtrapolate to x = { x_extrapolate }" )
print( f"y_exact = { y_extrapolate_exact }" )
print( f"y_fit   = { y_extrapolate_fit   }" )
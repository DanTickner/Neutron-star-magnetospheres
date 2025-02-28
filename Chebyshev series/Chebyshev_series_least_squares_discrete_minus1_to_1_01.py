'''
Chebyshev_series_least_squares_discrete_minus1_to_1_01.py

Fit a set of discrete datapoints to a set of Chebyshev polynomials by least-squares fitting.
Useful for prototyping.
Only valid on the interval [-1,1] until I get it working.

Modified from
20230613 VSH decomposition/Fit_function_to_Legendre_rational_functions.py
and
Chebyshev_series_least_squares_01.py
and
Chebyshev_series_least_squares_discrete_02
The latter includes a test to prove that gaussian_elimination( A, B ) works.

V02: Do everything explicitly as a way of prototyping for C++ code.

'''

import numpy as np


def T_n( x, n ):
	# Chebyshev polynompial of order n evalauted at x, scaled for the interval [a, b]
	return np.cos( n * np.arccos( x ) )


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

x_data = [ -0.834, -0.432, 0.532, 0.866 ]
y_data = [ 3.0*x**3 - 1.5*x**2 for x in x_data ]



#----- Build matrix A -----
A = np.zeros( shape=(n_max+1,n_max+1) )

for n in range( n_max+1 ):
	for m in range( n_max+1 ):
		
		for i in range( len( x_data ) ):
			A[n][m] += T_n( x_data[i], n ) * T_n( x_data[i], m )


#----- Build vector B -----
B = np.zeros( n_max+1 )

for n in range( n_max+1 ):
	
	for i in range( len( x_data ) ):
		B[n] += y_data[i] * T_n( x_data[i], n )


#----- Perform Gaussian elimination -----
print( B )

gaussian_elimination( A, B )

B_expected = [ -0.75, 2.25, -0.75, 0.75 ]

print( "\nMatrix of coefficients, followed by expected matrix of coefficients:")
print( B )
print( B_expected )


#----- Calculate fitted y-values -----
y_fit = np.zeros( len( x_data ) )

for i in range( len( x_data ) ):
	
	for n in range( n_max+1 ):
		y_fit[i] += B[n] * T_n( x_data[i], n )


print()
print( B[0] * T_n(x_data[0],0) )
print( B[1] * T_n(x_data[0],1) )
print( B[2] * T_n(x_data[0],2) )
print( B[3] * T_n(x_data[0],3) )

#----- Output correct and fitted y-values -----

print( "\nx_data[i]\ty_data[i]\ty_fit[i]" )
for i in range( len( x_data ) ):
	print( f"{ x_data[i] }\t{ y_data[i] }\t{ y_fit[i] }" )


#----- Now try interpolating outward a bit -----
# Defining the Cheb poly's as cos(ncos^-1(x)) means we CANNOT extrapolate because the functions are not defined outside [-1,1].
# Instead, we might be able to use the recursion relation.
x_extrapolate = 0.965

y_extrapolate_exact = 3.0*x_extrapolate**3 - 1.5*x_extrapolate**2

y_extrapolate_fit = 0
for n in range( n_max+1 ):
	y_extrapolate_fit += B[n] * T_n( x_extrapolate, n )

print( f"\nExtrapolate to x = { x_extrapolate }" )
print( f"y_exact = { y_extrapolate_exact }" )
print( f"y_fit   = { y_extrapolate_fit   }" )




#----- Generate Chebyshev polynomials by recursion relation -----
# An issue with this is that the polynomials grow/plummet unboundedly beyond [-1,1].
# This is because any power of x (power greater than 1) does this.
# With "not-nice" data, i.e. not engineered to follow some polynomial anyway, we will probably get significant error by extrapolating.

x_extrapolate_2 = 1.965

T_n_recursion_list = np.zeros( n_max+1 )

T_n_recursion_list[0] = 1
T_n_recursion_list[1] = x_extrapolate_2

for n in range( 2, n_max+1 ):
	T_n_recursion_list[n] = 2.0 * x_extrapolate_2 * T_n_recursion_list[n-1] - T_n_recursion_list[n-2]


print()
print( T_n_recursion_list )

y_extrapolate_2_exact = 3.0*x_extrapolate_2**3 - 1.5*x_extrapolate_2**2

y_extrapolate_2_fit = 0
for n in range( n_max+1 ):
	y_extrapolate_2_fit += B[n] * T_n_recursion_list[n]

print( f"\nExtrapolate to x = { x_extrapolate_2 }" )
print( f"y_exact = { y_extrapolate_2_exact }" )
print( f"y_fit   = { y_extrapolate_2_fit   }" )
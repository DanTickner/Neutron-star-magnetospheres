'''
Chebyshev_series_least_squares_discrete_01.py

Fit a set of discrete datapoints to a set of Chebyshev polynomials by least-squares fitting.
Useful for prototyping.

Modified from
20230613 VSH decomposition/Fit_function_to_Legendre_rational_functions.py
and
Chebyshev_series_least_squares_01.py

'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#----- Define functions -----
def fit_function( x, a, b, c, d, e ):
	T0 = 1
	T1 = x
	T2 = 2*x*x - 1
	T3 = 4*x**3 - 3*x
	T4 = 8*x**4 - 8*x**2 + 1
	return a*T0 + b*T1 + c*T2 + d*T3 + e*T4

def y_function( x ):
	return 1.0 / x

def T_n( x, n, a, b ):
	# Chebyshev polynompial of order n evalauted at x, scaled for the interval [a, b].
	Lambda_inverse = ( 2.0 - x - a - b ) / ( b - a )
	print( Lambda_inverse )
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
	for i in range( len(B)-1, -1 ):
		dotproduct = sum( [ A[i][j] * B[j] for j in range( i+1, len(B) ) ] )
		B[i] = ( B[i] - dotproduct ) / A[i][i]
	
	return B



x_values = [ 79.5255, 79.6046, 79.6837, 79.7628, 79.8418 ]
y_values = [ -1.47344E-07, -1.46760E-07, -1.46178E-07, -1.45599E-07, -1.45023E-07 ]

n_max = 3	# Highest fitting order


print( T_n(3,2,2,5) )


#----- Define variables -----
# Change the multiplication factors to simulate extrapolation.
#x_min = min( x_values ) * 1.0
#x_max = max( x_values ) * 1.3
n_points = 100


#----- Build matrix A -----
A = np.zeros( shape=(n_max+1,n_max+1) )
x_min = 1
x_max = 80

for n in range( n_max+1 ):
	A[n][n] = sum( [ T_n(x_values[i],n,x_min,x_max)**2 for i in range(len( x_values )) ] )
	for m in range( n ):
		A[n][m] = sum( [ T_n(x_values[i],n,x_min,x_max) * T_n(x_values[i],m,x_min,x_max) for i in range(len( x_values )) ] )
		A[m][n] = A[n][m]

# Divide by two because the a_0 coefficient is halved. Divide A[0][0] by 4 in total.
for n in range( n_max+1 ):
	A[0][n] *= 0.5
	A[n][0] *= 0.5

B = np.zeros( n_max+1 )
for n in range( n_max+1 ):
	B[n] = sum( [ y_values[i] * T_n(x_values[i],n,x_min,x_max) for i in range(len( x_values )) ] )


C = gaussian_elimination( A, B )

#----- Perform the fit -----
#popt, pcov = curve_fit( fit_function, x_list, y_data_list )
#y_fit_list = fit_function( x_list, *popt )
# *popt is a shortcut to apply all elements of popt as arguments to the function.
# See the Examples section of the webpage.

#print( "Fitted coefficients" )
#print( popt )


#----- Plot graph -----
fig = plt.figure()
ax  = fig.add_subplot( 111 )

#ax.plot( x_list, y_fit_list , label="Fit"  )
ax.scatter( x_values, y_values, label="Data" )

ax.legend( loc="best" )
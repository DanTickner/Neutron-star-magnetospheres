'''
Chebyshev_series_least_squares_extrapolate_01.py

Fit f(x) to a set of Chebyshev polynomials by least-squares fitting.
Useful if the integrals for the coefficients do not converge.
Extrapolate the fit

Modified from
20230613 VSH decomposition/Fit_function_to_Legendre_rational_functions.py

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



#----- Define variables -----
x_min = 2
x_max = 3
x_max_extrapolate = 4
n_points = 100


#----- Build lists of data -----
x_list = np.linspace( x_min, x_max, n_points )
y_data_list = y_function( x_list )

x_list_extrapolate = np.linspace( x_min, x_max_extrapolate, n_points )
y_data_list_extrapolate = y_function( x_list_extrapolate )


#----- Perform the fit -----
popt, pcov = curve_fit( fit_function, x_list, y_data_list )
y_fit_list = fit_function( x_list, *popt )
# *popt is a shortcut to apply all elements of popt as arguments to the function.
# See the Examples section of the webpage.

y_fit_list_extrapolate = fit_function( x_list_extrapolate, *popt )

print( "Fitted coefficients" )
print( popt )


#----- Plot graph -----
fig = plt.figure()
ax  = fig.add_subplot( 111 )

#ax.plot( x_list, y_data_list, label="Data" )
#ax.plot( x_list, y_fit_list , label="Fit"  )

ax.plot( x_list_extrapolate, y_data_list_extrapolate, label="Data extrapolated" )
ax.plot( x_list_extrapolate, y_fit_list_extrapolate , label="Fit extrapolated"  )

ax.legend( loc="best" )
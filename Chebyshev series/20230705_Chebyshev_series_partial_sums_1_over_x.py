'''

20230705_Chebyshev_series_partial_sums_1_over_x.py

We have expressions for the Chebyshev series coefficients of f(x)=1/x.
Evaluate the partial sums of Chebyshev series coefficients
(i.e. truncations of the series at a given n_max)
and see if the approximation gets closer and closer to 1/x.
'''

import numpy as np
import matplotlib.pyplot as plt

#----- Define functions -----

def a_n( n ):
	if n % 4 == 1:
		return 2
	if n % 4 == 3:
		return -2
	return 0

def T_n( x, n ):
	return np.cos( n * np.arccos( x ) )


#----- Define variables -----
x = 0.321	# Value at which to evaluate the series. x in [-1,1].
n_max = 10	# Largest value of n at which to truncate the series.


#----- Build lists -----
a_n_list = []
T_n_list = []

for n in range( n_max + 1 ):
	a_n_list.append( a_n(n)   )
	T_n_list.append( T_n(x,n) )


#----- Calculate approximate values -----
f_exact = 1.0 / x
f_guess_list = []
f_error_list = []

for n in range( n_max + 1 ):
	
	f_guess = 0
	for k in range( n + 1 ):
		f_guess += a_n_list[k] * T_n_list[k]
	
	f_guess_list.append( f_guess )
	f_error_list.append( f_exact - f_guess )



#----- Plot graph -----
plt.rcParams.update({ "font.size":18 })
plt.rc( "text", usetex=True )

fig = plt.figure()
ax  = fig.add_subplot( 111 )

ax.plot( range( n_max+1 ), f_error_list )
ax.axhline( 0, ls="--", color="grey", lw=0.5 )
ax.set_xlabel( r"$n_{\rm{max}}$" )
ax.set_ylabel( r"$f(x)-f_{\rm{Chebyshev}}(x,n)$" )
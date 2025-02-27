'''
20250226_Chebyshev_arbitrary_interval_test.py

Check that I have the right expressions for the Chebyshev polynomials on an arbitrary interval.

'''
import numpy as np
import matplotlib.pyplot as plt

def Lambda( x, A, B ):
	return 0.5 * (B-A) * x + 0.5 * (A+B)

def Lambda_inverse( x, A, B ):
	return 2.0/(B-A) * X - (B+A)/(B-A)


A = 2
B = 7

n_points = 1000

X_list = np.linspace( A, B, n_points )

T_2_list = []

for X in X_list:
	T_2_list.append( 2.0 * ( Lambda_inverse(X,A,B) )**2 - 1.0 )


#----- Plot graph -----
plt.rc( "text", usetex=True )
plt.rcParams.update({ "font.size":18 })

fig = plt.figure( figsize=(8,5) )
ax  = fig.add_subplot( 111 )

ax.plot( X_list, T_2_list )
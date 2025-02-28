'''
Check mathematical expressions/Check_P_ell_m_integral_RHS_at_0_20230803.py
Plot the expression for the integral

int_0^pi P_ell^m[cos(theta)] sin(theta) dtheta = 
f(m) [ (-1)^ell+(-1)^M ] 2^{M-2} M Gamma(ell/2) Gamma( [ell+M+1]/2 ) / ( Gamma( [ell-m+2]/2 ) Gamma( [ell+3]/2 )
M = abs(m)
f(m) = 1 if m >= 0 and (-1)^{m+1} Gamma(ell+m+1) / Gamma(ell-m+1) if m < 0.
There is an issue at ell=m=0 so take the limit.
Plot for ell in [0,1] and m in [-1,1].
Don't consider the left limit because ell>=0.
'''

from scipy.special import gamma
import numpy as np

def function( ell, m ):
	
	
	M = abs(m)
	
	f = 1
	if m < 0:
		f = (-1)**(m+1) * gamma( ell+m+1 ) / gamma( ell-m+1 )
	
	#return f * ( (-1)**ell + (-1)**M ) * 2**(M-2) * M * gamma(0.5*ell) * gamma( 0.5*(ell+M+1) ) / ( gamma( 0.5*(ell-M+2) ) * gamma( 0.5*(ell+3) ) )
	
	return f * ( (-1)**ell + (-1)**M ) * 2**(M-2) * gamma( 0.5*(ell+M+1) ) / ( gamma( 0.5*(ell-M+2) ) * gamma( 0.5*(ell+3) ) )
	
	
	#return abs(m) * gamma( 0.5*ell )



n_points_ell = 100
n_points_m   = 100

values = np.zeros( shape=(n_points_ell,n_points_m), dtype=complex )
# Need complex because Gamma(-x) in C for x not an integer.

delta_ell = ( 0   - (-1.0) ) / ( n_points_ell - 1.0 )
delta_m   = ( 1.0 - (-1.0) ) / ( n_points_m   - 1.0 )

ell_list = [ 0  + i * delta_ell for i in range( n_points_ell ) ]
m_list   = [ -1 + i * delta_m   for i in range( n_points_m   ) ]

for i in range( n_points_ell ):
	for j in range( n_points_m ):
		values[i][j] = function( ell_list[i], m_list[j] )
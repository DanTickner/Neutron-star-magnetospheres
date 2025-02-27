'''
Decimal_search_NS_neutron_number_density.py

Find the value of x for which f(x) = 0, given a defined function f(x).
Start from an initial guess x_0 and increase accuracy up to a specified number of decimal places.

Application of the function to find the equilibrium neutron number density in a neutron star.
This uses the equation as given in (7.1) of Ryan-Norton 2010, which implicitly requires a mass density.
Once working, derive an equation for a mass fraction so that the mass density explicitly appears.
'''

import numpy as np # np.sign(), np.pi


#----- Define variables -----
m_p = 1.672621926E-27
m_n = 1.674927498E-27
m_e = 9.109383714E-31
h   = 6.62607015E-34
c   = 299792458.0

rho          = 1.3E17 # Mass density in kg m^-3. Free parameter.
#x = rho / ( m_n + m_p + m_e )	# First value at which to evaluate the function.
x = 5.003e43
tens_initial = int( np.log10( x / 10 ) )	# A in 10^A to use first in the search. Es.g. A=-1 means we initially jump by 0.1 each time.
n_sig_figs   = 7	# Desired number of significant figures in result.

print( f"Initial x: { x }" )
print( f"Initial x-jump: { tens_initial }\n" )



#----- Define function -----
def f(x):
	
	n_p = ( rho - m_n * x ) / ( m_p + m_e )
	
	'''
	term_1 = ( 3.0 * n_p / ( 8.0 * np.pi ) )**(1.0/3.0) * h * c
	term_2 = ( 3.0 * n_p / ( 8.0 * np.pi ) )**(2.0/3.0) * h*h / ( 2.0 * m_p )
	term_3 = ( 3.0 * x   / ( 8.0 * np.pi ) )**(2.0/3.0) * h*h / ( 2.0 * m_n ) * -1.0
	term_4 = ( m_n - m_p ) * c*c * -1.0
	return term_1 + term_2 + term_3 + term_4
	'''
	
	'''
	ret = n_p**(2.0/3.0) / m_p - x**(2.0/3.0) / m_n
	ret *= ( 3.0 / ( 8.0 * np.pi ) )**(2.0/3.0) * h*h * 0.5
	ret += ( 3.0 * n_p / ( 8.0 * np.pi ) )**(1.0/3.0) * h * c
	ret -= ( m_n - m_p ) * c*c
	return ret
	'''

	ret = x**(2.0/3.0) / m_n - n_p**(2.0/3.0) / m_p
	ret *= ( 3.0 / ( 8.0 * np.pi ) )**(2.0/3.0) * h*h * 0.5
	ret -= ( 3.0 * n_p / ( 8.0 * np.pi ) )**(1.0/3.0) * h * c
	ret += ( m_n - m_p ) * c*c
	return ret


#----- Set initial values -----
dx           = 10**tens_initial
f_i          = f( x - dx )		# Start one point below, so that the loop begins at f(x).
f_prev       = f( x - dx - dx )	# To store the previous value so that sign changes can be checked. Choose initial value that is likely to have same sign as f.

#----- Perform decimal search -----
for n in range( n_sig_figs - 1 ):
	
	for i in range( 11 ):
		# Go to 11 so that the jump 0.9 to 1.0 is included.
		x += dx
		f_i = f( x )
		
		if ( np.sign( f_i ) != np.sign( f_prev ) ) or ( abs( f_i.imag ) > 0 ):
			x -= dx
			dx *= 0.1
			break
		
		f_prev = f_i


#----- Output final values -----
print( f"Lower estimate: x = { x       }\tf(x) = { f(x      ) }" )
print( f"Upper estimate: x = { x+10*dx }\tf(x) = { f(x+10*dx) }\n" )

x_avg = ( x + (x+10*dx) ) * 0.5
n_p = n_p = ( rho - m_n * x_avg ) / ( m_p + m_e )

print( f"Average: x = { x_avg }" )
print( f"Neutron-to-proton ratio: { x_avg / n_p }\t:\t1" )
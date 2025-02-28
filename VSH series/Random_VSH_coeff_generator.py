'''
Random_VSH_coeff_generator.py
Generate a set of random complex VSH series coefficients in [0,1] with three decimal places.
'''


import random


'''
n_coeffs = 10


for n in range( n_coeffs ):
	print( f"{ round( random.random(),3 ) }\t{ round( random.random(),3 ) }" )
	
'''



ell_max = 3

for ell in range( ell_max+1 ):
	for m in range( -ell, ell+1 ):
		
		string = f"if( ( ell == { ell } ) and ( m == "
		
		if m >= 0:
			string += " "
		
		string += f"{ m } ) ){{ return std::complex<double> {{ "
		# Use double curly braces as an escape character for curly braces with f-strings.
		
		string_a = str( round( random.random(),3 ) )
		
		if len(string_a) < 5:
			string_a += "0" * ( 5 - len(string_a) )
		
		string += f"{ string_a }, "
			
		string_b = str( round( random.random(),3 ) )
		
		if len(string_b) < 5:
			string_b += "0" * ( 5 - len(string_b) )
			
		string += f"{ string_b } }}; }}"
		
		print( string )
	print()
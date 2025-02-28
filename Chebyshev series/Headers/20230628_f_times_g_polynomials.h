/*
20230628_f_times_g_polynomials.h

Functions f and g for which we wish to determine a Chebyshev series of fg in terms of the series of f and of g.

f(x) = 3x^4 + 4x^3 - 2x^2 -  x + 1
g(x) =  x^4 + 5x^3 + 3x^2 - 2x + 3

Coeffs found by matching powers of x.
Note that a_{n_max} and a_{n_max-1} can be found immediately because they are the only coeffs multiplying x^{n_max} and x^{n_max-1}. The other a_n are uncovered two-at-a-time because the Tn are even, then odd.
See photo 2023068 22:00.

*/

double f( double x ){
	return 3.0*pow(x,4) + 4.0*pow(x,3) - 2.0*pow(x,2) - x + 1.0;
}




double g( double x ){
	return pow(x,4) + 5.0*pow(x,3) + 3.0*pow(x,2) - 2.0*x + 3;
}




double fg( double x ){
	return f(x) * g(x);
}




double fg_explicit( double x ){
	return 3.0*pow(x,8) + 19.0*pow(x,7) + 27.0*pow(x,6) - 5.0*pow(x,5) - 9.0*pow(x,4) + 18.0*pow(x,3) - pow(x,2) - 5.0*x + 3.0;
}






double f_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 9.0/8;
		case 1: return 2.0;
		case 2: return 1.0/2;
		case 3: return 1.0;
		case 4: return 3.0/8;
		
		default: return 0;
		
	}
	
}




double g_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 39.0/8;
		case 1: return 7.0/4;
		case 2: return 2.0;
		case 3: return 5.0/4;
		case 4: return 1.0/8;
		
		default: return 0;
		
	}
	
}




double fg_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 1073.0/128;
		case 1: return 1009.0/64;
		case 2: return 287.0/32;
		case 3: return 587.0/64;
		case 4: return 147.0/32;
		case 5: return 113.0/64;
		case 6: return 33.0/32;
		case 7: return 19.0/64;
		case 8: return 3.0/128;
		
		default: return 0;
		
	}
	
}
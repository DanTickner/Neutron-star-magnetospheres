/*
20230703_f_times_g_other_functions.h

Functions f and g for which we wish to determine a Chebyshev series of fg in terms of the series of f and of g.

f(x) = sin(3x)
g(x) = e^(-x/2)

Haven't bothered to find the exact Chebyshev coefficients of f, g or fg.
This code is more a test that it can be applied to unknown, not-exactly-decomposable-in-finitely-many-terms functions, and still get accurate results.

*/

double f( double x ){
	return sin( 3.0 * x );
}




double g( double x ){
	return exp( -0.5 * x );
}




double fg( double x ){
	return f(x) * g(x);
}




double fg_explicit( double x ){
	return sin( 3.0 * x ) * exp( -0.5 * x );
}






double f_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}




double g_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}




double fg_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}
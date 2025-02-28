/*
20230705_x_times_f_other_functions.h

Functions f and g for which we wish to determine a Chebyshev series of fg in terms of the series of f and of g.

f(x) = x
g(x) = e^(1.3*x) + 3/2

Haven't bothered to find the exact Chebyshev coefficients of f, g or fg.
This code is more a test that it can be applied to unknown, not-exactly-decomposable-in-finitely-many-terms functions, and still get accurate results.

*/

double f( double x ){
	return x;
}




double g( double x ){
	return exp( 1.3 * x ) + 1.5;
}




double fg( double x ){
	return f(x) * g(x);
}




double fg_explicit( double x ){
	return exp( 1.3 * x ) + 1.5;
}






double f_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 0;
		
		default: return 0;
		
	}
	
}




double g_n_guess( int n ){
	
	switch( n ){
		
		case 0: return 2.96928;
		case 1: return 1.59466;
		case 2: return 0.485235;
		case 3: return 0.101629;
		case 4: return 0.0161777;
		case 5: return 0.00207417;
		case 6: return 0.000222481;
		case 7: return 2.05051e-5;
		
		default: return 0;
		
	}
	
}




double fg_n_guess( int n ){
	
	if( n == 0 ){
		return 0.5 * g_n_guess( 1 );
	}
	if( n == 1 ){
		return g_n_guess( 0 ) + 0.5 * g_n_guess( 2 );
	}
	if( n > 1 ){
		return 0.5 * g_n_guess( n-1 ) + 0.5 * g_n_guess( n+1 );
	}
	
	return 0;
	
}